GCC = g++
OPT = -O3
LIBS = -lconfig++ -lfftw3
CULIBS = -lcurand

OBJ = \
obj/main.o \
obj/NeighbourList.o \
obj/mathfuncs.o \
obj/config.o \
obj/geom.o \
obj/error.o \
obj/spinwaves.o \
obj/util.o

CULIBS = \
obj/cumalloc.o \
obj/cuthermal.o \
obj/cuheun.o \
obj/cuintegrate.o

obj/%.o: src/%.cpp
	$(GCC) $(OPT) -c -o $@ $< 

obj/%.o: src/%.cu
	$(NVCC) -DCUDA -G -c -o $@ $<


ASD: $(OBJ)
	$(GCC) $(OPT) -o $@ $^ $(LIBS)

clean:
	@rm -f obj/*.o

