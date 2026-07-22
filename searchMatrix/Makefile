# Build searchmatrix from src/.
#
# THIS FILE EXISTS BECAUSE THE DEPLOYED BINARY WENT STALE. Commit 8a853a7
# (2026-07-15, "test-before-build in searchM hot loop, 2.9x") was committed and
# compiled, but `build/searchmatrix` was never refreshed -- it stayed at the
# 2026-07-11 build. Every sweep afterwards ran the pre-optimisation binary while
# the optimised one sat unused beside it. Measured cost: 38.4s vs 10.3s on one
# chunk, a 3.7x tax paid silently for a week.
#
# -O2 IS NOT OPTIONAL. The historical hand-written command (README, RUNBOOK)
# carried no -O flag, so the binary was built at -O0. That is a further 3.8x
# (10.3s -> 2.7s). Hit sets are byte-identical at every optimisation level, so
# this is pure throughput. The two factors compound to ~14x.
#
#   make            build if any source is newer than the binary
#   make stale      exit 1 if the binary is missing or older than src/ (no build)
#   make rebuild    force a clean rebuild
#   make clean      remove the binary only; *.bak snapshots are preserved
#
# `make stale` is the guard: submit_inverted.sh calls it before dispatching a
# sweep, so a stale binary fails at submit time instead of quietly producing a
# week of slow -- or wrong -- results.

CXX      ?= g++
# -DSILENT suppresses the per-hit chatter the original code writes to stdout.
CXXFLAGS ?= -O2 -DSILENT
SRC      := src/searchMatrix.cpp
HDRS     := $(wildcard src/*.h)
BIN      := build/searchmatrix

.PHONY: all stale rebuild clean

all: $(BIN)

# Depend on every header: searchMatrix.cpp is one translation unit that #includes
# them all, so a change to searchControl.h must trigger a rebuild.
$(BIN): $(SRC) $(HDRS)
	@mkdir -p build
	$(CXX) $(CXXFLAGS) $(SRC) -o $@
	@echo "built $@ with: $(CXX) $(CXXFLAGS)"

stale:
	@if [ ! -x $(BIN) ]; then \
	    echo "STALE: $(BIN) does not exist -- run 'make' in $(CURDIR)" >&2; exit 1; \
	fi; \
	newer=$$(find $(SRC) $(HDRS) -newer $(BIN) 2>/dev/null); \
	if [ -n "$$newer" ]; then \
	    echo "STALE: $(BIN) is older than:" >&2; \
	    echo "$$newer" | sed 's/^/  /' >&2; \
	    echo "run 'make' in $(CURDIR) before sweeping" >&2; exit 1; \
	fi; \
	echo "$(BIN) is up to date with src/"

rebuild:
	@rm -f $(BIN)
	@$(MAKE) --no-print-directory $(BIN)

# Deliberately does NOT touch build/*.bak -- those are the provenance trail for
# which binary produced which sweep.
clean:
	rm -f $(BIN)
