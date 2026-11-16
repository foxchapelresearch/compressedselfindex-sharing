# Needs a recent clang (C++26), libsais (brew install llvm libsais) and the FSST and CRoaring
# submodules (git submodule update --init).
LLVM    ?= $(shell brew --prefix llvm 2>/dev/null)
SAIS    ?= $(shell brew --prefix libsais 2>/dev/null)
CC      := $(if $(LLVM),$(LLVM)/bin/clang,clang)
CXX     := $(if $(LLVM),$(LLVM)/bin/clang++,clang++)
FLAGS   := -O2 -g -Wall -Wextra -Iroaring/include $(if $(SAIS),-I$(SAIS)/include)
LIBS    := $(if $(SAIS),-L$(SAIS)/lib) -lsais
FSST    := build/fsst/libfsst.o build/fsst/fsst_avx512.o
ROARING := $(patsubst %.c,build/%.o,$(wildcard roaring/src/*.c roaring/src/containers/*.c roaring/src/art/*.c))

csi: compressed_self_index.cpp $(FSST) $(ROARING)
	$(CXX) -std=c++26 $(FLAGS) $< $(FSST) $(ROARING) $(LIBS) -o $@

build/fsst/%.o: fsst/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) -std=c++17 -O3 -w -c $< -o $@

build/roaring/%.o: roaring/%.c
	@mkdir -p $(dir $@)
	$(CC) -std=c11 -O3 -w -Iroaring/include -c $< -o $@

clean:
	rm -rf csi build *.dSYM

.PHONY: clean
