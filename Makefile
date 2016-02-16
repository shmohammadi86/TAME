CXX=g++
INCLUDE=-I./include -I./armadillo/include
LIB_FLAGS=-lblas -llapack #-L./armadillo -larmadillo
CXXFLAGS=$(INCLUDE) -march=native -g -O3 -funroll-loops -msse2  -Wall -Wno-write-strings
SRC=src/main.cc src/graph.cc src/io.cc src/triangle.cc src/tensor.cc
OBJ=$(SRC:.cc=.o)
PROGRAM=tri-match

all : $(PROGRAM)

src/%.o: src/%.cc src/%.h
	$(CXX) $(CXXFLAGS) -c -o $@  $<

$(PROGRAM): $(OBJ) 
	$(CXX) $(CXXFLAGS)  -o $@ $(OBJ) $(LIB_FLAGS)	

.PHONY: clean ar

clean:
	rm -f $(PROGRAM)  src/*.o src/*~ include/*~ *~

ar:
	make clean
	tar -czvf ../$(PROGRAM)"(`date`)".tar.gz *
