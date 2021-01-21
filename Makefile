GCC = g++
OPT = -O3
LIBS = -lconfig++ -lfftw3

OBJ = \
obj/main.o \
obj/NeighbourList.o \
obj/mathfuncs.o \
obj/config.o \
obj/geom.o \
obj/error.o \
obj/spinwaves.o


obj/%.o: src/%.cpp
	$(GCC) $(OPT) -c -o $@ $< 

ASD: $(OBJ)
	$(GCC) $(OPT) -o $@ $^ $(LIBS)

clean:
	@rm -f obj/*.o

