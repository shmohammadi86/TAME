#===============================================================================
# TAME -- Triangular AlignMEnt
# https://github.com/shmohammadi86/TAME
#
# Dependencies: Armadillo (>= 6), BLAS, LAPACK, and an OpenMP-capable C++17
# compiler. See README.md for platform-specific install commands.
#===============================================================================

CXX       ?= g++
PREFIX    ?= /usr/local

OPTFLAGS  ?= -O2
OPENMP    ?= -fopenmp
WARNINGS   = -Wall -Wno-unused-variable -Wno-unused-but-set-variable \
             -Wno-sign-compare -Wno-write-strings -Wno-unused-result

CPPFLAGS  += -I./include
CXXFLAGS  ?= -std=c++17 $(OPTFLAGS) $(WARNINGS) $(OPENMP)
LDLIBS    ?= -larmadillo -llapack -lblas
LDFLAGS   += $(OPENMP)

PROGRAM    = tri-match

SRC        = src/main.cc src/graph.cc src/io.cc src/triangle.cc src/tensor.cc
BMATCH_SRC = src/bMatch/mtxReader.cpp src/bMatch/bSuitor.cpp \
             src/bMatch/bSuitorD.cpp src/bMatch/PG.cpp \
             src/bMatch/localDom.cpp src/bMatch/PGDP.cpp \
             src/bMatch/Node.cpp

OBJ        = $(SRC:.cc=.o) $(BMATCH_SRC:.cpp=.o)

.PHONY: all check check-badinput install clean distclean

all: $(PROGRAM)

$(PROGRAM): $(OBJ)
	$(CXX) $(LDFLAGS) -o $@ $(OBJ) $(LDLIBS)

%.o: %.cc
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c -o $@ $<

# Smoke test: a small synthetic pair, so this finishes in seconds. Use
# ./sparse_run.sh for a full alignment of the bundled interactomes.
check: $(PROGRAM)
	@rm -rf output && mkdir -p output
	./$(PROGRAM) -t smat \
	    -G tests/data/G_small.smat \
	    -H tests/data/H_small.smat \
	    -S tests/data/GH_small.smat \
	    -x seqsim -Y 2 -a 0.85 -b 0.1 --iter 2 --post_iter 1 -o output
	@ls output/* >/dev/null 2>&1 || { echo "FAIL: no output produced"; exit 1; }
	@echo "OK: TAME smoke test passed"

# A headerless similarity file must be rejected, not read out of bounds.
check-badinput: $(PROGRAM)
	@rm -rf output && mkdir -p output
	@printf '1 2 0.5\n3 4 0.25\n' > /tmp/tame_bad.evals
	@./$(PROGRAM) -t smat -G tests/data/G_small.smat -H tests/data/H_small.smat \
	    -S /tmp/tame_bad.evals -x uniform -Y 2 --iter 1 --post_iter 0 -o output \
	    > /tmp/tame_bad.log 2>&1; \
	rc=$$?; \
	if [ $$rc -ge 128 ]; then \
	    echo "FAIL: crashed on malformed input (signal $$((rc - 128)))"; exit 1; \
	fi; \
	grep -q 'readSeqSim_SMAT::' /tmp/tame_bad.log \
	    && echo "OK: malformed similarity file rejected, no crash" \
	    || { echo "FAIL: malformed similarity file silently accepted"; exit 1; }

install: $(PROGRAM)
	install -d $(DESTDIR)$(PREFIX)/bin
	install -m 0755 $(PROGRAM) $(DESTDIR)$(PREFIX)/bin/$(PROGRAM)

clean:
	rm -f $(PROGRAM) $(OBJ) *~ src/*~ include/*~

distclean: clean
	rm -rf output
