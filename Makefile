CXX=g++
INCLUDE=-I./include -I./armadillo/include
LIB_FLAGS=-lblas -llapack #-L./armadillo -larmadillo
CXXFLAGS=$(INCLUDE) -g3 -march=native -O4 -std=c++11 -fomit-frame-pointer -funroll-loops -fforce-addr -fexpensive-optimizations -msse2 -fopenmp
SRC=src/main.cc src/graph.cc src/io.cc src/triangle.cc src/tensor.cc
OBJ=$(SRC:.cc=.o)
PROGRAM=tri-match
bMatch_OBJECTS = \
	src/bMatch/mtxReader.o \
	src/bMatch/bSuitor.o \
	src/bMatch/bSuitorD.o \
	src/bMatch/PG.o \
	src/bMatch/localDom.o \
	src/bMatch/PGDP.o \
	src/bMatch/Node.cpp \

	
all : $(PROGRAM) message

src/%.o: src/%.cc src/%.h
	$(CXX) $(CXXFLAGS) -c -o $@  $<

src/bMatch/%.o: src/bMatch/%.cc
	$(CXX) $(CXXFLAGS) -c -o $@  $<


$(PROGRAM): $(OBJ) $(bMatch_OBJECTS)
	$(CXX) $(CXXFLAGS)  -o $@ $(OBJ) $(bMatch_OBJECTS) $(LIB_FLAGS)	

.PHONY: clean ar

clean:
	rm -f $(PROGRAM)  src/*.o src/*~ include/*~ *~

ar:
	make clean
	tar -czvf ../$(PROGRAM)"(`date`)".tar.gz *

message:
	echo "Executable: $(PROGRAM) has been created"	
