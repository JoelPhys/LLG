GCC = g++
OPT = -O3
LIBS = -lfftw3

OBJ = \
../obj/lswt.o \
../obj/mathfuncs.o \


../obj/%.o: ../src/%.cpp
	$(GCC) $(OPT) -c -o $@ $< 

lswtexe: $(OBJ)
	$(GCC) $(OPT) -o $@ $^ $(LIBS)

clean:
	@rm -f obj/*.o

