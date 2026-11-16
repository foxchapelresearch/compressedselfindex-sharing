// compressed_self_index.cpp -- a chunked FM-index over (optionally FSST-compressed) text.
//
// Build:  make   (needs libsais; FSST and CRoaring are git submodules)
//
// HOW IT WORKS
//
// 1. FM-index. Sort all suffixes of a text T. The Burrows-Wheeler transform (BWT) lists, for
//    each suffix in sorted order ("row"), the byte that precedes it in T. All suffixes that
//    start with a pattern P occupy a contiguous range of rows, and the range for cP follows
//    from the range for P with two rank queries (rank(c, i) = number of c's in bwt[0..i)).
//    So P is found by feeding its bytes right to left ("backward search").
//
// 2. PFO20 tiles. The BWT is stored in 4096-byte tiles, one page each, so that everything a
//    query step needs for a row is on the same page:
//        640 bytes   256 rank counters of 20 bits: occurrences of each byte value before this
//                    tile, relative to the start of its superpage (256 tiles)
//        384 bytes   one bit per row: is the suffix array sampled here?
//       3072 bytes   the BWT bytes themselves
//    A small table holds the absolute counts at the start of every superpage. rank is one
//    table lookup, one 20-bit counter and a scan of at most 3072 bytes.
//
// 3. Locate and extract. LF(row) is the row of the suffix starting one byte earlier. Every
//    sample_rate-th text position is sampled: its row is marked in the tile, and sa_samples[k]
//    is the row of position k * sample_rate. To extract text we start from the sample after
//    the range and walk LF backwards, reading one BWT byte per step. To locate a row we walk
//    LF to a marked row; the *locate cache* (a separate, optional file: a rank directory over
//    the mark bits plus the sample numbers in row order) then tells us which sample it is.
//    Without the cache that last step is a linear scan over all samples.
//
// 4. FSST. Optionally T is first compressed with FSST: up to 255 symbols of 1..8 bytes are
//    replaced by one-byte codes, any other byte b becomes the pair (ESC, b). The FM-index is
//    built over the code stream. A pattern then has many possible spellings, depending on
//    where symbol boundaries fall. The query planner enumerates them as "plans" (a code
//    sequence plus the offset of the pattern inside its first symbol), each plan is searched
//    for, and each hit is mapped back to a text position (via stream/text offset pairs stored
//    every 4096 codes) and verified by decoding. Pruning policies cut down the number of plans
//    by assuming the encoder parses greedily.
//
// 5. Case folding. Optionally T is lowercased before indexing and the positions of capitals
//    are kept in a Roaring bitmap. Case-sensitive queries check each hit against the bitmap;
//    extraction restores case. Without folding, a case-insensitive query searches for every
//    upper/lower spelling of the pattern (at most 12 letters).
//
// 6. Chunks. The text is cut into chunks that overlap by 4096 bytes and are indexed
//    independently. A match that lies entirely inside the overlap at the start of a chunk was
//    already reported by the previous chunk and is dropped.
//
// The index is one file that is mmap'ed and used in place. It is written in host byte order.
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <memory>
#include <optional>
#include <print>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#include <libsais.h>
#include <roaring/roaring.h>

#include "fsst/fsst.h"

#ifdef _WIN32
#include <windows.h>
#else
#include <fcntl.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#endif

namespace csi {

using u8 = std::uint8_t;
using u32 = std::uint32_t;
using u64 = std::uint64_t;
using i32 = std::int32_t;
using i64 = std::int64_t;
using Bytes = std::span<const u8>;

constexpr u32  tile_bytes = 4096, tile_header_bytes = 640, tile_marks_bytes = 384, tile_payload_bytes = 3072;
constexpr u32  superpage_tiles = 256;
constexpr u32  overlap_bytes = 4096;        // chunk overlap; also the longest pattern
constexpr u32  locate_block_bits = 4096;    // rows per entry of the locate cache's rank directory
constexpr u32  fsst_pos_sample_rate = 4096; // codes between stream/text offset pairs
constexpr u32  decode_window = 256;         // stream bytes fetched at a time when decoding
constexpr u64  section_align = 4096;
constexpr u32  flag_fsst = 1, flag_fold = 2;
constexpr char magic[8] = "CSIDX02";

struct FileHeader {
    char magic[8];
    u32  flags, sa_sample_rate;
    u64  text_size, chunk_bytes, nchunks, chunk_table;
    u8   fsst_dict[2072]; // the symbol table as serialized by fsst_export()

    bool fsst() const { return flags & flag_fsst; }
    bool fold() const { return flags & flag_fold; }
};

struct ChunkRecord {
    i64 file_offset; // where this chunk starts in the original text
    i64 raw_len;     // text bytes covered
    i64 n;           // indexed bytes (== raw_len unless FSST) == number of BWT rows
    i64 primary;     // row of the suffix that starts at position 0
    i64 sa_sample_rate, sa_sample_count;
    u64 num_tiles, pos_sample_count, uppercase_bytes;
    u64 base_table, tiles, sa_samples, pos_samples, uppercase; // section file offsets
    i64 c_table[257];                                          // c_table[c] = bytes smaller than c in the chunk
};

struct BuildOptions {
    bool fsst = false, fold_case = false;
    u32  sample_rate = 32;
    u64  chunk_bytes = 256ull << 20;
};

enum class Pruning { off, unsafe, safe };

constexpr u8   lower(u8 c) { return c >= 'A' && c <= 'Z' ? c + 32 : c; }
constexpr u8   upper(u8 c) { return c >= 'a' && c <= 'z' ? c - 32 : c; }
constexpr bool is_upper(u8 c) { return c >= 'A' && c <= 'Z'; }
constexpr bool is_letter(u8 c) { return lower(c) >= 'a' && lower(c) <= 'z'; }

class MappedFile {
public:
    explicit MappedFile(const std::string &path)
    {
#ifdef _WIN32
        HANDLE file = CreateFileA(path.c_str(), GENERIC_READ, FILE_SHARE_READ, nullptr, OPEN_EXISTING, 0, nullptr);
        if (file != INVALID_HANDLE_VALUE) {
            LARGE_INTEGER size;
            HANDLE        view = CreateFileMappingA(file, nullptr, PAGE_READONLY, 0, 0, nullptr);
            if (view && GetFileSizeEx(file, &size))
                if (auto *p = static_cast<const u8 *>(MapViewOfFile(view, FILE_MAP_READ, 0, 0, 0))) bytes_ = {p, size_t(size.QuadPart)};
            if (view) CloseHandle(view);
            CloseHandle(file);
        }
#else
        int         fd = open(path.c_str(), O_RDONLY);
        struct stat st;
        if (fd >= 0 && fstat(fd, &st) == 0 && st.st_size > 0)
            if (void *p = mmap(nullptr, size_t(st.st_size), PROT_READ, MAP_PRIVATE, fd, 0); p != MAP_FAILED)
                bytes_ = {static_cast<const u8 *>(p), size_t(st.st_size)};
        if (fd >= 0) close(fd);
#endif
        if (bytes_.empty()) throw std::runtime_error("cannot read " + path + " (missing or empty)");
    }
    MappedFile(const MappedFile &) = delete;
    ~MappedFile()
    {
#ifdef _WIN32
        UnmapViewOfFile(bytes_.data());
#else
        munmap(const_cast<u8 *>(bytes_.data()), bytes_.size());
#endif
    }
    Bytes bytes() const { return bytes_; }
    template <class T> const T *at(u64 offset) const { return reinterpret_cast<const T *>(bytes_.data() + offset); }

private:
    Bytes bytes_;
};

class Writer {
public:
    explicit Writer(const std::string &path) : out_(path, std::ios::binary)
    {
        if (!out_) throw std::runtime_error("cannot create " + path);
    }
    // Appends a section, padded to a page, and returns its file offset.
    template <class T> u64 put(std::span<const T> data)
    {
        const u64 at = offset_, bytes = data.size_bytes(), pad = -bytes % section_align;
        out_.write(reinterpret_cast<const char *>(data.data()), std::streamsize(bytes));
        out_.write(std::string(pad, '\0').data(), std::streamsize(pad));
        offset_ += bytes + pad;
        return at;
    }
    template <class T> u64 put(const std::vector<T> &data) { return put(std::span<const T>(data)); }
    template <class T> void finish(const T &header)
    {
        out_.seekp(0).write(reinterpret_cast<const char *>(&header), sizeof header).flush();
        if (!out_) throw std::runtime_error("error writing file");
    }

private:
    std::ofstream out_;
    u64           offset_ = 0;
};

// ---------------------------------------------------------------------------------------------
// PFO20: the BWT in page-sized tiles, each with its own rank counters and sample marks
// ---------------------------------------------------------------------------------------------

struct Pfo20 {
    const u32 *base_table = nullptr; // [superpages][256] occurrences before each superpage
    const u8  *tiles = nullptr;      // [num_tiles][tile_bytes]
    u64        n = 0;

    static u32 delta20(const u8 *tile, u8 c)
    {
        u64 word;
        std::memcpy(&word, tile + c * 20 / 8, 8);
        return u32(word >> (c * 20 % 8)) & 0xFFFFF;
    }
    // Number of occurrences of c in bwt[0..i).
    i64 rank(u8 c, i64 i) const
    {
        if (i <= 0) return 0;
        const u64 pos = std::min(u64(i), n);
        const u64 tile = (pos - 1) / tile_payload_bytes, within = pos - tile * tile_payload_bytes; // 1..3072
        const u8 *t = tiles + tile * tile_bytes, *payload = t + tile_header_bytes + tile_marks_bytes;
        return base_table[tile / superpage_tiles * 256 + c] + delta20(t, c) + std::count(payload, payload + within, c);
    }
    u8 access(i64 row) const
    {
        return tiles[u64(row) / tile_payload_bytes * tile_bytes + tile_header_bytes + tile_marks_bytes + u64(row) % tile_payload_bytes];
    }
    u8 mark_byte(u64 index) const // byte `index` of the conceptual mark bitvector
    {
        return tiles[index / tile_marks_bytes * tile_bytes + tile_header_bytes + index % tile_marks_bytes];
    }
    bool is_marked(i64 row) const { return mark_byte(u64(row) >> 3) >> (row & 7) & 1; }

    // Marked rows in [first, last]; first must be a multiple of 8.
    u32 marks_in(u64 first, u64 last) const
    {
        u32 count = 0;
        for (u64 byte = first >> 3; byte < last >> 3; byte++) count += u32(std::popcount(mark_byte(byte)));
        return count + u32(std::popcount(mark_byte(last >> 3) & ((2u << (last & 7)) - 1)));
    }

    struct Built {
        std::vector<u32> base_table;
        std::vector<u8>  tiles;
    };
    static Built build(Bytes bwt, std::span<const i32> marked_rows)
    {
        const u64            num_tiles = (bwt.size() + tile_payload_bytes - 1) / tile_payload_bytes;
        Built                out{{}, std::vector<u8>(num_tiles * tile_bytes)};
        std::array<u32, 256> global{}, superpage_start{};
        for (u64 tile = 0; tile < num_tiles; tile++) {
            if (tile % superpage_tiles == 0) {
                out.base_table.append_range(global);
                superpage_start = global;
            }
            u64 header[tile_header_bytes / 8] = {};
            for (u32 c = 0; c < 256; c++) {
                const u64 delta = global[c] - superpage_start[c], bit = c * 20;
                header[bit / 64] |= delta << (bit % 64);
                if (bit % 64 > 44) header[bit / 64 + 1] |= delta >> (64 - bit % 64);
            }
            const Bytes payload = bwt.subspan(tile * tile_payload_bytes, std::min<u64>(tile_payload_bytes, bwt.size() - tile * tile_payload_bytes));
            for (u8 c : payload) global[c]++;
            u8 *t = out.tiles.data() + tile * tile_bytes;
            std::memcpy(t, header, tile_header_bytes);
            std::ranges::copy(payload, t + tile_header_bytes + tile_marks_bytes);
        }
        for (u64 row : marked_rows)
            out.tiles[row / tile_payload_bytes * tile_bytes + tile_header_bytes + row % tile_payload_bytes / 8] |= u8(1u << (row % 8));
        return out;
    }
};

// ---------------------------------------------------------------------------------------------
// FSST symbol table and the query planner
// ---------------------------------------------------------------------------------------------

struct Fsst {
    u32                   nsym = 0;
    u8                    sym_len[256] = {};
    u8                    sym[256][8] = {};
    std::vector<u8>       codes_by_first[256]; // codes whose symbol starts with a given byte
    std::array<bool, 256> has_literal_code{};  // is there a one-byte symbol for this byte?
    u8                    max_sym_len = 1;

    explicit Fsst(const u8 *serialized)
    {
        fsst_decoder_t decoder;
        if (!fsst_import(&decoder, serialized)) throw std::runtime_error("bad fsst symbol table");
        for (int i = 0; i < 8; i++) nsym += serialized[9 + i]; // number of symbols of each length
        for (u32 c = 0; c < nsym; c++) {
            sym_len[c] = decoder.len[c];
            std::memcpy(sym[c], &decoder.symbol[c], 8);
            codes_by_first[sym[c][0]].push_back(u8(c));
            max_sym_len = std::max(max_sym_len, sym_len[c]);
            if (sym_len[c] == 1) has_literal_code[sym[c][0]] = true;
        }
    }
    Bytes symbol(u8 code) const { return {sym[code], sym_len[code]}; }
};

// One spelling of a pattern in the code stream: the pattern starts start_offset bytes into
// the text that codes[0] decodes to, and may end inside the symbol of the last code.
struct Plan {
    std::vector<u8> codes;
    u8              start_offset;
};

class Planner {
public:
    Planner(const Fsst &fsst, Bytes pattern) : f_(fsst), pat_(pattern), buf_(pattern.size() * 2) {}

    std::vector<Plan> run(Pruning policy)
    {
        const bool allow_boundary_split = policy == Pruning::safe && pat_.size() <= 511;

        // What may the unit covering pattern[0] be when it starts exactly at pattern[0]?
        //   off:    anything that matches.
        //   unsafe: only what a greedy encoder emits for the pattern itself.
        //   safe:   what a greedy encoder emits for any prefix of the pattern (the encoder may
        //           have been looking at a shorter remainder than we assume).
        std::array<bool, 256> start_allow{};
        bool                  start_allow_escape = policy == Pruning::off;
        if (policy == Pruning::off) start_allow.fill(true);
        for (size_t d = 1, max_d = std::min<size_t>(pat_.size(), f_.max_sym_len); policy == Pruning::safe && d <= max_d; d++) {
            const Emission e = greedy(pat_.first(d));
            (e.code_len == 2 ? start_allow_escape : start_allow[e.codes[0]]) = true;
        }
        const Emission boundary = greedy(pat_);
        if (policy == Pruning::unsafe) (boundary.code_len == 2 ? start_allow_escape : start_allow[boundary.codes[0]]) = true;

        auto continue_from = [&](size_t consumed, u32 len, u8 offset) {
            policy == Pruning::off ? dfs(consumed, len, offset) : pruned_continue(consumed, len, offset, allow_boundary_split);
        };

        // First unit: any symbol, entered at any offset, whose tail matches the pattern's head.
        for (u32 code = 0; code < f_.nsym; code++) {
            for (u8 offset = 0; offset < f_.sym_len[code]; offset++) {
                const size_t consumed = matches(u8(code), pat_, offset);
                if (!consumed) continue;
                buf_[0] = u8(code);
                if (consumed == pat_.size())
                    emit(1, offset);
                else if (offset != 0 || (start_allow[code] && (policy != Pruning::unsafe || consumed == boundary.raw_len)))
                    continue_from(consumed, 1, offset);
            }
        }
        if (start_allow_escape && !f_.has_literal_code[pat_[0]]) {
            buf_[0] = FSST_ESC, buf_[1] = pat_[0];
            pat_.size() == 1 ? emit(2, 0) : continue_from(1, 2, 0);
        }
        return std::move(out_);
    }

private:
    struct Emission {
        u8 codes[2], code_len, raw_len;
    };

    void emit(u32 len, u8 start_offset) { out_.push_back({{buf_.begin(), buf_.begin() + len}, start_offset}); }

    // How many bytes of `remaining` does the symbol's text from `offset` on match? Zero unless
    // it matches all the way to the end of the symbol or of `remaining`, whichever is first.
    size_t matches(u8 code, Bytes remaining, u8 offset) const
    {
        const Bytes  tail = f_.symbol(code).subspan(offset);
        const size_t len = std::min(tail.size(), remaining.size());
        return std::equal(tail.begin(), tail.begin() + i64(len), remaining.begin()) ? len : 0;
    }

    // What a greedy encoder writes next for `remaining`: the longest symbol, else an escape.
    Emission greedy(Bytes remaining) const
    {
        Emission best{{FSST_ESC, remaining[0]}, 2, 1};
        for (u8 code : f_.codes_by_first[remaining[0]])
            if (const u8 len = f_.sym_len[code]; len <= remaining.size() && matches(code, remaining, 0) == len && (best.code_len == 2 || len > best.raw_len))
                best = {{code, 0}, 1, len};
        return best;
    }

    // Every way of spelling pattern[pos..) given buf_[0..len).
    void dfs(size_t pos, u32 len, u8 start_offset)
    {
        const Bytes remaining = pat_.subspan(pos);
        if (!f_.has_literal_code[remaining[0]]) {
            buf_[len] = FSST_ESC, buf_[len + 1] = remaining[0];
            remaining.size() == 1 ? emit(len + 2, start_offset) : dfs(pos + 1, len + 2, start_offset);
        }
        for (u8 code : f_.codes_by_first[remaining[0]]) {
            const size_t consumed = matches(code, remaining, 0);
            if (!consumed) continue;
            buf_[len] = code;
            consumed == remaining.size() ? emit(len + 1, start_offset) : dfs(pos + consumed, len + 1, start_offset);
        }
    }

    // The pattern may end inside a symbol that is longer than what is left of it.
    void terminal_partial_symbols(Bytes remaining, u32 len, u8 start_offset)
    {
        for (u8 code : f_.codes_by_first[remaining[0]]) {
            if (f_.sym_len[code] <= remaining.size() || !matches(code, remaining, 0)) continue;
            buf_[len] = code;
            emit(len + 1, start_offset);
        }
    }

    void greedy_segment(Bytes segment, u32 &len)
    {
        for (size_t pos = 0; pos < segment.size();) {
            const Emission e = greedy(segment.subspan(pos));
            std::copy_n(e.codes, e.code_len, buf_.begin() + len);
            len += e.code_len, pos += e.raw_len;
        }
    }

    // Spell pattern[pos..) the way a greedy encoder would. With allow_boundary_split, also try
    // cutting each multi-byte symbol short once: the text after the pattern can make the real
    // encoder choose a different (longer) symbol further on, which shifts the boundaries here.
    void pruned_continue(size_t pos, u32 len, u8 start_offset, bool allow_boundary_split)
    {
        while (pos < pat_.size()) {
            const Bytes remaining = pat_.subspan(pos);
            terminal_partial_symbols(remaining, len, start_offset);
            const Emission e = greedy(remaining);
            for (u8 split = 1; allow_boundary_split && e.code_len == 1 && split < e.raw_len; split++) {
                u32 split_len = len;
                greedy_segment(pat_.subspan(pos, split), split_len);
                pruned_continue(pos + split, split_len, start_offset, false);
            }
            std::copy_n(e.codes, e.code_len, buf_.begin() + len);
            len += e.code_len, pos += e.raw_len;
        }
        emit(len, start_offset);
    }

    const Fsst       &f_;
    Bytes             pat_;
    std::vector<u8>   buf_;
    std::vector<Plan> out_;
};

// ---------------------------------------------------------------------------------------------
// Building
// ---------------------------------------------------------------------------------------------

using Encoder = std::unique_ptr<fsst_encoder_t, decltype(&fsst_destroy)>;

// Trains a symbol table on 4096-byte windows of the first 64 MiB and stores it in the header.
Encoder train(FileHeader &h, Bytes text)
{
    constexpr size_t window = 4096;
    std::vector<u8>  sample(std::from_range, text.first(std::min<size_t>(text.size(), 64 << 20)));
    if (h.fold()) std::ranges::transform(sample, sample.begin(), lower);
    std::vector<size_t>     lengths;
    std::vector<const u8 *> windows;
    for (size_t i = 0; i < sample.size(); i += window) {
        windows.push_back(sample.data() + i);
        lengths.push_back(std::min(window, sample.size() - i));
    }
    Encoder encoder(fsst_create(windows.size(), lengths.data(), windows.data(), 0), fsst_destroy);
    fsst_export(encoder.get(), h.fsst_dict);
    return encoder;
}

std::vector<u8> encode(fsst_encoder_t *encoder, Bytes text)
{
    std::vector<u8> out(2 * text.size() + 7);
    size_t          in_len = text.size(), out_len;
    const u8       *in = text.data();
    u8             *start;
    if (fsst_compress(encoder, 1, &in_len, &in, out.size(), out.data(), &out_len, &start) != 1) throw std::runtime_error("fsst failed");
    out.resize(out_len);
    return out;
}

ChunkRecord build_chunk(const FileHeader &h, const Fsst *fsst, fsst_encoder_t *encoder, Bytes raw, i64 file_offset, Writer &w)
{
    ChunkRecord rec{};
    rec.file_offset = file_offset;
    rec.raw_len = i64(raw.size());

    Bytes           text = raw; // the bytes that actually get indexed
    std::vector<u8> folded, codes;
    if (h.fold()) {
        std::unique_ptr<roaring_bitmap_t, decltype(&roaring_bitmap_free)> capitals(roaring_bitmap_create(), roaring_bitmap_free);
        for (size_t i = 0; i < raw.size(); i++)
            if (is_upper(raw[i])) roaring_bitmap_add(capitals.get(), u32(i));
        roaring_bitmap_run_optimize(capitals.get());
        std::vector<char> frozen(roaring_bitmap_frozen_size_in_bytes(capitals.get()));
        roaring_bitmap_frozen_serialize(capitals.get(), frozen.data());
        rec.uppercase_bytes = frozen.size();
        rec.uppercase = w.put(frozen);
        folded.assign_range(raw | std::views::transform(lower));
        text = folded;
    }
    if (h.fsst()) {
        text = codes = encode(encoder, text);
        // Every fsst_pos_sample_rate codes, remember (stream offset, text offset); plus the end.
        std::vector<u64> stream_offsets, text_offsets;
        u64              text_pos = 0, ncodes = 0;
        for (size_t i = 0; i < codes.size(); ncodes++) {
            if (ncodes % fsst_pos_sample_rate == 0) stream_offsets.push_back(i), text_offsets.push_back(text_pos);
            text_pos += codes[i] == FSST_ESC ? 1 : fsst->sym_len[codes[i]];
            i += codes[i] == FSST_ESC ? 2 : 1;
        }
        stream_offsets.push_back(codes.size()), text_offsets.push_back(text_pos);
        rec.pos_sample_count = stream_offsets.size();
        stream_offsets.append_range(text_offsets);
        rec.pos_samples = w.put(stream_offsets);
    }
    const i64 n = i64(text.size());
    if (n >= INT32_MAX) throw std::runtime_error("chunk too large; use a smaller chunk size");

    // The sample rate is rounded up to a power of two; tiny chunks sample everything.
    const i64 rate = n <= 64 ? 1 : i64(std::bit_ceil(std::min<u64>(h.sa_sample_rate, u64(n))));
    rec.n = n;
    rec.sa_sample_rate = rate;
    rec.sa_sample_count = (n + rate - 1) / rate;

    // BWT over the cyclic text: the row of the suffix at position 0 ("primary") holds the last
    // byte of the text instead of an end marker. sa_samples[k] = row of position k * rate.
    const size_t     rows = size_t(n), samples = size_t(rec.sa_sample_count);
    std::vector<i32> suffix(rows), sa_samples(samples);
    if (libsais(text.data(), suffix.data(), i32(n), 0, nullptr) != 0) throw std::runtime_error("suffix sorting failed");
    std::vector<u8> bwt(rows);
    for (i64 row = 0; row < n; row++) {
        const i64 pos = suffix[size_t(row)];
        if (pos == 0) rec.primary = row;
        bwt[size_t(row)] = text[size_t(pos == 0 ? n - 1 : pos - 1)];
        if (pos % rate == 0) sa_samples[size_t(pos / rate)] = i32(row);
    }
    for (u8 c : bwt) rec.c_table[c + 1]++;
    for (int c = 0; c < 256; c++) rec.c_table[c + 1] += rec.c_table[c];

    const auto pfo = Pfo20::build(bwt, sa_samples);
    rec.num_tiles = pfo.tiles.size() / tile_bytes;
    rec.base_table = w.put(pfo.base_table);
    rec.tiles = w.put(pfo.tiles);
    rec.sa_samples = w.put(sa_samples);
    return rec;
}

void build(const std::string &text_path, const std::string &index_path, BuildOptions opt)
{
    if (opt.sample_rate == 0 || opt.sample_rate > 1024) throw std::runtime_error("sample rate must be 1..1024");
    if (opt.chunk_bytes <= overlap_bytes || opt.chunk_bytes > 1ull << 30) throw std::runtime_error("chunk size must be 4097 bytes .. 1 GiB");
    const MappedFile file(text_path);
    const Bytes      text = file.bytes();
    const u64        stride = opt.chunk_bytes - overlap_bytes;

    FileHeader h{};
    std::ranges::copy(magic, h.magic);
    h.flags = (opt.fsst ? flag_fsst : 0) | (opt.fold_case ? flag_fold : 0);
    h.sa_sample_rate = opt.sample_rate;
    h.text_size = text.size();
    h.chunk_bytes = opt.chunk_bytes;
    h.nchunks = 1 + (text.size() - 1) / stride;
    const Encoder       encoder = opt.fsst ? train(h, text) : Encoder(nullptr, fsst_destroy);
    std::optional<Fsst> fsst;
    if (opt.fsst) fsst.emplace(h.fsst_dict);

    Writer w(index_path);
    w.put(std::span<const FileHeader>(&h, 1));
    std::vector<ChunkRecord> table;
    for (u64 i = 0; i < h.nchunks; i++) {
        const Bytes chunk = text.subspan(i * stride, std::min(opt.chunk_bytes, text.size() - i * stride));
        table.push_back(build_chunk(h, fsst ? &*fsst : nullptr, encoder.get(), chunk, i64(i * stride), w));
    }
    h.chunk_table = w.put(table);
    w.finish(h);
}

// ---------------------------------------------------------------------------------------------
// Querying
// ---------------------------------------------------------------------------------------------

struct SearchOptions {
    bool    ignore_case = false;
    u64     max_results = 0; // 0 = unlimited
    Pruning pruning = Pruning::off;
};

class Index {
public:
    explicit Index(const std::string &path) : path_(path), file_(path), h_(*file_.at<FileHeader>(0))
    {
        const u64 size = file_.bytes().size();
        if (size < sizeof h_ || !std::ranges::equal(h_.magic, magic) || h_.chunk_table > size ||
            h_.nchunks > (size - h_.chunk_table) / sizeof(ChunkRecord))
            throw std::runtime_error(path + " is not an index file");
        if (h_.fsst()) fsst_.emplace(h_.fsst_dict);
        for (u64 i = 0; i < h_.nchunks; i++) {
            const ChunkRecord &rec = file_.at<ChunkRecord>(h_.chunk_table)[i];
            Chunk              ch(rec, fsst_ ? &*fsst_ : nullptr);
            ch.pfo = {file_.at<u32>(rec.base_table), file_.at<u8>(rec.tiles), u64(rec.n)};
            ch.sentinel_char = ch.pfo.access(rec.primary);
            ch.sa_samples = {file_.at<i32>(rec.sa_samples), size_t(rec.sa_sample_count)};
            ch.stream_offsets = {file_.at<u64>(rec.pos_samples), rec.pos_sample_count};
            ch.text_offsets = {ch.stream_offsets.data() + rec.pos_sample_count, rec.pos_sample_count};
            if (h_.fold()) ch.capitals = roaring_bitmap_frozen_view(file_.at<char>(rec.uppercase), rec.uppercase_bytes);
            chunks_.push_back(ch);
        }
        try {
            load_locate_cache();
        } catch (const std::runtime_error &) { // no usable cache: locate falls back to scanning the samples
            cache_.reset();
            for (Chunk &ch : chunks_) ch.locate_rank_dir = ch.locate_sample_ids = {};
        }
    }
    Index(const Index &) = delete;
    ~Index()
    {
        for (Chunk &ch : chunks_) roaring_bitmap_free(ch.capitals); // frees only the view
    }

    u64  text_size() const { return h_.text_size; }
    bool has_locate_cache() const { return bool(cache_); }

    // Positions of all occurrences of the pattern, ascending.
    std::vector<u64> search(std::string_view pattern, const SearchOptions &options) const;

    // text[position .. position+length), clipped to the text.
    std::vector<u8> extract(u64 position, u64 length) const
    {
        std::vector<u8> out;
        const u64       stride = h_.chunk_bytes - overlap_bytes;
        for (length = std::min(length, h_.text_size - std::min(position, h_.text_size)); length;) {
            const Chunk &ch = chunks_[std::min(position / stride, h_.nchunks - 1)];
            const i64    a = i64(position) - ch.rec.file_offset, len = std::min(i64(length), ch.rec.raw_len - a);
            out.append_range(ch.read_text(a, a + len));
            position += u64(len), length -= u64(len);
        }
        return out;
    }

    // Writes the locate cache next to the index: for every chunk, the number of marked rows
    // before each block of locate_block_bits rows, then the sample numbers in row order.
    void build_locate_cache() const
    {
        Writer w(path_ + ".locate");
        w.put(std::span<const u64>(std::array<u64, 2>{locate_block_bits, h_.nchunks}));
        for (const Chunk &ch : chunks_) {
            std::vector<u32> ids(std::from_range, std::views::iota(0u, u32(ch.sa_samples.size())));
            std::ranges::sort(ids, {}, [&](u32 id) { return ch.sa_samples[id]; });
            std::vector<u32> rank_dir;
            for (u64 block = 0, j = 0; block <= u64(ch.rec.n) / locate_block_bits; block++) {
                while (j < ids.size() && u64(ch.sa_samples[ids[j]]) < block * locate_block_bits) j++;
                rank_dir.push_back(u32(j));
            }
            w.put(rank_dir);
            w.put(ids);
        }
    }

    void print_info() const;

private:
    struct Chunk {
        const ChunkRecord      &rec;
        const Fsst             *fsst;
        Pfo20                   pfo;
        u8                      sentinel_char = 0; // the byte stored at the primary row (the text's last byte)
        std::span<const i32>    sa_samples;        // [k] = row of the suffix at position k * sa_sample_rate
        std::span<const u64>    stream_offsets, text_offsets; // FSST position samples
        const roaring_bitmap_t *capitals = nullptr;
        std::span<const u32>    locate_rank_dir, locate_sample_ids; // from the locate cache, if present

        Chunk(const ChunkRecord &record, const Fsst *symbols) : rec(record), fsst(symbols) {}

        // The row of the suffix one byte earlier in the (cyclic) text. Because the primary row
        // stores a real byte where the end marker would be, that byte needs special care.
        i64 lf(i64 row) const
        {
            const u8 c = pfo.access(row);
            if (row == rec.primary) return rec.c_table[c];
            return rec.c_table[c] + pfo.rank(c, row) + (c == sentinel_char && row < rec.primary);
        }

        // Rows [sp, ep) whose suffixes start with the pattern. The text's last byte is stored
        // at the primary row but really precedes the (absent) empty suffix, which sorts before
        // row 0: so it counts as an occurrence before every row, and not at the primary row.
        // The initial range includes the empty suffix, hence no correction on the first step.
        std::pair<i64, i64> backward_search(Bytes pattern) const
        {
            i64 sp = 0, ep = rec.n;
            for (size_t i = pattern.size(); i-- > 0 && sp < ep;) {
                const u8  c = pattern[i];
                const i64 old_sp = sp, old_ep = ep;
                sp = rec.c_table[c] + pfo.rank(c, old_sp) + (c == sentinel_char && old_sp <= rec.primary && i + 1 < pattern.size());
                ep = rec.c_table[c] + pfo.rank(c, old_ep) + (c == sentinel_char && old_ep <= rec.primary);
            }
            return {sp, ep};
        }

        // The position in the indexed bytes of the suffix at `row`, or -1.
        i64 locate(i64 row) const
        {
            for (i64 steps = 0; steps < rec.sa_sample_rate; steps++, row = lf(row)) {
                if (!pfo.is_marked(row)) continue;
                i64 sample = -1;
                if (!locate_rank_dir.empty()) {
                    const u64 block = u64(row) / locate_block_bits;
                    sample = locate_sample_ids[locate_rank_dir[block] + pfo.marks_in(block * locate_block_bits, u64(row)) - 1];
                } else if (auto it = std::ranges::find(sa_samples, row); it != sa_samples.end()) {
                    sample = it - sa_samples.begin();
                }
                return sample < 0 ? -1 : (sample * rec.sa_sample_rate + steps) % rec.n;
            }
            return -1;
        }

        // Indexed bytes [start, end): for each sample that follows part of the range, walk LF
        // backwards from it for up to sa_sample_rate steps.
        std::vector<u8> read_indexed(i64 start, i64 end) const
        {
            start = std::clamp<i64>(start, 0, rec.n), end = std::clamp(end, start, rec.n);
            std::vector<u8> out(size_t(end - start));
            const i64       rate = rec.sa_sample_rate;
            for (i64 k = start / rate + 1; k <= std::min((end + rate - 1) / rate, rec.sa_sample_count); k++) {
                i64 row = k < rec.sa_sample_count ? sa_samples[size_t(k)] : rec.primary;
                i64 pos = (k < rec.sa_sample_count ? k * rate : rec.n) - 1;
                for (i64 steps = 0; steps < rate && pos >= start; steps++, pos--, row = lf(row))
                    if (pos < end) out[size_t(pos - start)] = pfo.access(row);
            }
            return out;
        }

        // FSST: the text offset that stream offset `code_pos` decodes to.
        i64 stream_to_text(i64 code_pos) const
        {
            if (code_pos <= 0) return 0;
            if (code_pos >= rec.n) return rec.raw_len;
            const size_t idx = size_t(std::ranges::upper_bound(stream_offsets, u64(code_pos)) - stream_offsets.begin()) - 1;
            u64          text_pos = text_offsets[idx];
            const auto   codes = read_indexed(i64(stream_offsets[idx]), code_pos);
            for (size_t i = 0; i < codes.size(); i++) {
                if (codes[i] == FSST_ESC && i + 1 >= codes.size()) break;
                text_pos += codes[i] == FSST_ESC ? 1 : fsst->sym_len[codes[i]];
                i += codes[i] == FSST_ESC;
            }
            return i64(text_pos);
        }

        // FSST: text bytes [a, b), decoded window by window starting from the nearest position
        // sample at or before a.
        std::vector<u8> decode_text(i64 a, i64 b) const
        {
            a = std::clamp<i64>(a, 0, rec.raw_len), b = std::clamp(b, a, rec.raw_len);
            std::vector<u8> out;
            const size_t    idx = size_t(std::ranges::upper_bound(text_offsets, u64(a)) - text_offsets.begin()) - 1;
            i64             code_pos = i64(stream_offsets[idx]), text_pos = i64(text_offsets[idx]);
            while (text_pos < b && code_pos < rec.n) {
                const auto codes = read_indexed(code_pos, code_pos + decode_window);
                size_t     i = 0;
                while (i < codes.size() && text_pos < b) {
                    const bool escape = codes[i] == FSST_ESC;
                    if (escape && i + 1 >= codes.size()) break; // its literal is in the next window
                    const Bytes unit = escape ? Bytes(&codes[i + 1], 1) : fsst->symbol(codes[i]);
                    for (u8 byte : unit)
                        if (text_pos++ >= a && text_pos <= b) out.push_back(byte);
                    i += escape ? 2 : 1;
                }
                if (i == 0) break;
                code_pos += i64(i);
            }
            return out;
        }

        // Original text bytes [a, b) of the chunk, with capitals restored.
        std::vector<u8> read_text(i64 a, i64 b) const
        {
            a = std::max<i64>(a, 0);
            auto out = fsst ? decode_text(a, b) : read_indexed(a, b);
            for (size_t i = 0; capitals && i < out.size(); i++)
                if (roaring_bitmap_contains(capitals, u32(a) + u32(i))) out[i] = upper(out[i]);
            return out;
        }
    };

    void load_locate_cache()
    {
        cache_ = std::make_unique<MappedFile>(path_ + ".locate");
        if (cache_->at<u64>(0)[0] != locate_block_bits || cache_->at<u64>(0)[1] != h_.nchunks) throw std::runtime_error("stale locate cache");
        const auto page_align = [](u64 bytes) { return (bytes + section_align - 1) / section_align * section_align; };
        u64        offset = section_align;
        for (Chunk &ch : chunks_) {
            const u64 blocks = u64(ch.rec.n) / locate_block_bits + 1;
            ch.locate_rank_dir = {cache_->at<u32>(offset), blocks};
            offset += page_align(blocks * 4);
            ch.locate_sample_ids = {cache_->at<u32>(offset), size_t(ch.rec.sa_sample_count)};
            offset += page_align(ch.locate_sample_ids.size_bytes());
        }
        if (offset != cache_->bytes().size()) throw std::runtime_error("stale locate cache");
    }

    std::string                 path_;
    MappedFile                  file_;
    const FileHeader           &h_;
    std::optional<Fsst>         fsst_;
    std::vector<Chunk>          chunks_;
    std::unique_ptr<MappedFile> cache_;
};

std::vector<u64> Index::search(std::string_view pattern, const SearchOptions &options) const
{
    if (pattern.empty() || pattern.size() > overlap_bytes) throw std::runtime_error("pattern must be 1..4096 bytes");
    const Bytes original(reinterpret_cast<const u8 *>(pattern.data()), pattern.size());
    const i64   m = i64(original.size());

    // The spellings to search for. On a folded index: the lowercased pattern. On an unfolded
    // index with a case-insensitive query: every upper/lower combination of its letters.
    std::vector<std::vector<u8>> variants(1, {std::from_range, original});
    if (h_.fold() || options.ignore_case) std::ranges::transform(variants[0], variants[0].begin(), lower);
    if (options.ignore_case && !h_.fold()) {
        std::vector<size_t> letters;
        for (size_t i = 0; i < original.size(); i++)
            if (is_letter(original[i])) letters.push_back(i);
        if (letters.size() > 12) throw std::runtime_error("case-insensitive search on a case-sensitive index allows at most 12 letters");
        variants.resize(size_t(1) << letters.size(), variants[0]);
        for (size_t v = 0; v < variants.size(); v++)
            for (size_t bit = 0; bit < letters.size(); bit++)
                if (v >> bit & 1) variants[v][letters[bit]] = upper(variants[v][letters[bit]]);
    }
    const bool verify_case = h_.fold() && !options.ignore_case;

    std::vector<std::vector<Plan>> plans; // FSST: per variant
    for (const auto &variant : variants)
        if (fsst_) plans.push_back(Planner(*fsst_, variant).run(options.pruning));

    std::vector<u64> matches;
    bool             stop = false;
    for (size_t chunk_index = 0; chunk_index < chunks_.size() && !stop; chunk_index++) {
        const Chunk &ch = chunks_[chunk_index];
        const i64    first_match = i64(matches.size());

        auto push_match = [&](i64 text_pos) {
            // A match entirely inside the leading overlap belongs to the previous chunk.
            if (chunk_index > 0 && text_pos + m <= i64(overlap_bytes)) return;
            for (i64 i = 0; verify_case && i < m; i++)
                if (is_letter(original[size_t(i)]) && is_upper(original[size_t(i)]) != roaring_bitmap_contains(ch.capitals, u32(text_pos + i))) return;
            matches.push_back(u64(ch.rec.file_offset + text_pos));
        };

        for (size_t v = 0; v < variants.size(); v++) {
            if (!fsst_) {
                const auto [sp, ep] = ch.backward_search(variants[v]);
                for (i64 row = sp; row < ep; row++)
                    if (const i64 pos = ch.locate(row); pos >= 0) push_match(pos);
                continue;
            }
            for (const Plan &plan : plans[v]) {
                const auto [sp, ep] = ch.backward_search(plan.codes);
                for (i64 row = sp; row < ep; row++) {
                    const i64 code_pos = ch.locate(row);
                    if (code_pos < 0) continue;
                    // The hit may be a coincidence (e.g. it starts on the literal byte of an
                    // escape), so decode the text there and compare.
                    const i64 text_pos = ch.stream_to_text(code_pos) + plan.start_offset;
                    if (text_pos + m > ch.rec.raw_len || !std::ranges::equal(ch.decode_text(text_pos, text_pos + m), variants[v])) continue;
                    push_match(text_pos);
                }
            }
        }
        std::ranges::sort(matches.begin() + first_match, matches.end());
        matches.erase(std::unique(matches.begin() + first_match, matches.end()), matches.end());
        if (options.max_results && matches.size() >= options.max_results) matches.resize(options.max_results), stop = true;
    }
    return matches;
}

void Index::print_info() const
{
    u64 indexed = 0, tiles = 0, base = 0, samples = 0, pos_samples = 0, capitals = 0;
    for (const Chunk &ch : chunks_) {
        indexed += u64(ch.rec.n);
        tiles += ch.rec.num_tiles * tile_bytes;
        base += (ch.rec.num_tiles + superpage_tiles - 1) / superpage_tiles * 1024;
        samples += u64(ch.rec.sa_sample_count) * 4;
        pos_samples += ch.rec.pos_sample_count * 16;
        capitals += ch.rec.uppercase_bytes;
    }
    const u64    cache = cache_ ? cache_->bytes().size() : 0;
    const double pct = 100.0 / double(h_.text_size);
    std::println("text size         {:12} bytes", h_.text_size);
    std::println("index size        {:12} bytes  {:5.1f}% of text", file_.bytes().size(), double(file_.bytes().size()) * pct);
    std::println("codec             {}{}", h_.fsst() ? "fsst" : "raw", h_.fold() ? ", case-folded" : "");
    std::println("chunks            {:12}  ({} bytes each)", h_.nchunks, h_.chunk_bytes);
    std::println("sample rate       {:12}", chunks_[0].rec.sa_sample_rate);
    std::println("indexed bytes     {:12}  {:5.1f}%", indexed, double(indexed) * pct);
    for (auto [name, bytes] : std::initializer_list<std::pair<const char *, u64>>{
             {"pfo20 tiles", tiles}, {"superpage table", base}, {"sa samples", samples}, {"fsst positions", pos_samples}, {"capitals", capitals}})
        std::println("  {:<15} {:12}  {:5.1f}%", name, bytes, double(bytes) * pct);
    std::println("locate cache      {:12}  {:5.1f}%{}", cache, double(cache) * pct, cache ? "" : "  (not built; run `csi cache`)");
}

// Builds every index flavour for a text file and compares search and extract with a direct
// scan of the text. Meant for files of around a megabyte.
int selftest(const std::string &text_path, int queries)
{
    const MappedFile  file(text_path);
    const Bytes       text = file.bytes();
    const std::string index_path = text_path + ".selftest.csi";
    u64               state = 88172645463325252ull;
    auto              rng = [&] { return state ^= state << 13, state ^= state >> 7, state ^= state << 17; };
    int               failures = 0;
    for (int flavour = 0; flavour < 8; flavour++) {
        const BuildOptions opt{bool(flavour & 1), bool(flavour & 2), flavour & 4 ? 8u : 32u, flavour & 4 ? text.size() / 5 + 5000 : 256ull << 20};
        std::remove((index_path + ".locate").c_str());
        build(text_path, index_path, opt);
        if (flavour & 4) Index(index_path).build_locate_cache();
        const Index index(index_path);
        const int   before = failures;
        if (!std::ranges::equal(index.extract(0, text.size()), text)) failures++, std::println("  FAIL full extract");
        for (int i = 0; i < queries; i++) {
            const u64 pos = rng() % text.size(), len = 1 + rng() % 300;
            if (!std::ranges::equal(index.extract(pos, len), text.subspan(pos, std::min(len, text.size() - pos)))) failures++, std::println("  FAIL extract {}+{}", pos, len);
        }
        for (int i = 0; i < queries; i++) {
            const size_t m = 2 + rng() % (i % 4 == 0 ? 3 : 14);
            if (m >= text.size()) continue;
            const size_t at = rng() % (text.size() - m);
            std::string  pattern(reinterpret_cast<const char *>(text.data()) + at, m);
            if (i % 5 == 1) pattern[rng() % m] ^= 32;         // flip a letter's case
            if (i % 7 == 2) pattern[rng() % m] = char(rng()); // probably no match
            for (bool ignore_case : {false, true}) {
                std::vector<u64> want, got;
                for (size_t p = 0; p + m <= text.size(); p++)
                    if (std::ranges::equal(text.subspan(p, m), pattern, [&](u8 a, char b) { return ignore_case ? lower(a) == lower(u8(b)) : a == u8(b); })) want.push_back(p);
                SearchOptions options;
                options.ignore_case = ignore_case;
                try {
                    got = index.search(pattern, options);
                } catch (const std::runtime_error &) { // more than 12 letters
                    continue;
                }
                if (got != want) failures++, std::println("  FAIL {}\"{}\": want {}, got {}", ignore_case ? "-i " : "", pattern, want.size(), got.size());
            }
        }
        std::println("{:<4} {:<6} {:<15} {}", opt.fsst ? "fsst" : "raw", opt.fold_case ? "folded" : "", flavour & 4 ? "chunked + cache" : "", failures == before ? "ok" : "FAILED");
    }
    std::remove(index_path.c_str());
    std::remove((index_path + ".locate").c_str());
    return failures != 0;
}

} // namespace csi

// ---------------------------------------------------------------------------------------------
// Command line
// ---------------------------------------------------------------------------------------------

static int usage()
{
    std::println(stderr, "usage: csi build [--fsst] [--fold-case] [--sample-rate N] [--chunk-bytes N] TEXT INDEX\n"
                         "       csi cache    INDEX                           (build the locate cache, INDEX.locate)\n"
                         "       csi count    [-i] [--pruning off|unsafe|safe] INDEX PATTERN\n"
                         "       csi search   [-i] [--pruning ...] [-n MAX_MATCHES] [-c CONTEXT_BYTES] INDEX PATTERN\n"
                         "       csi extract  INDEX [POSITION [LENGTH]]       (whole text by default)\n"
                         "       csi info     INDEX\n"
                         "       csi selftest TEXT [QUERIES]");
    return 2;
}

int main(int argc, char **argv)
try {
    using csi::u64;
    const std::string        command = argc > 1 ? argv[1] : "";
    std::vector<std::string> args;
    csi::BuildOptions        build;
    csi::SearchOptions       search;
    u64                      context = 40;
    for (int i = 2; i < argc; i++) {
        const std::string_view a = argv[i];
        const auto             value = [&] { return std::strtoull(argv[++i], nullptr, 0); };
        const bool             has_value = i + 1 < argc;
        if (a == "--fsst") build.fsst = true;
        else if (a == "--fold-case") build.fold_case = true;
        else if (a == "--sample-rate" && has_value) build.sample_rate = csi::u32(value());
        else if (a == "--chunk-bytes" && has_value) build.chunk_bytes = value();
        else if (a == "-i") search.ignore_case = true;
        else if (a == "-n" && has_value) search.max_results = value();
        else if (a == "-c" && has_value) context = value();
        else if (a == "--pruning" && has_value) {
            const std::string_view policy = argv[++i];
            search.pruning = policy == "safe" ? csi::Pruning::safe : policy == "unsafe" ? csi::Pruning::unsafe : csi::Pruning::off;
        } else args.emplace_back(a);
    }

    if (command == "build" && args.size() == 2) {
        csi::build(args[0], args[1], build);
        return 0;
    }
    if (command == "selftest" && !args.empty()) return csi::selftest(args[0], args.size() > 1 ? std::stoi(args[1]) : 100);
    if (args.empty() || !std::ranges::contains(std::array<std::string_view, 5>{"cache", "count", "search", "extract", "info"}, command)) return usage();
    const csi::Index index(args[0]);

    if (command == "info" && args.size() == 1) {
        index.print_info();
    } else if (command == "cache" && args.size() == 1) {
        index.build_locate_cache();
    } else if (command == "count" && args.size() == 2) {
        std::println("{}", index.search(args[1], search).size());
    } else if (command == "search" && args.size() == 2) {
        for (u64 pos : index.search(args[1], search)) {
            const u64 from = pos - std::min(pos, context);
            auto      line = index.extract(from, pos - from + args[1].size() + context);
            std::ranges::replace_if(line, [](csi::u8 c) { return c < 32 || c == 127; }, ' ');
            std::println("{}\t{}", pos, std::string_view(reinterpret_cast<const char *>(line.data()), line.size()));
        }
    } else if (command == "extract" && args.size() <= 3) {
        u64 pos = args.size() > 1 ? std::stoull(args[1]) : 0, left = args.size() > 2 ? std::stoull(args[2]) : UINT64_MAX;
        for (std::vector<csi::u8> piece; left && !(piece = index.extract(pos, std::min<u64>(left, 1 << 16))).empty();) {
            std::fwrite(piece.data(), 1, piece.size(), stdout);
            pos += piece.size(), left -= piece.size();
        }
    } else {
        return usage();
    }
    return 0;
} catch (const std::exception &error) {
    std::println(stderr, "csi: {}", error.what());
    return 1;
}
