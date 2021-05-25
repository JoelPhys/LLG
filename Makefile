GCC = g++
OPT = -O3
LIBS = -lconfig++ -lfftw3
CULIBS = -lcurand
NVCC = nvcc

# OBJ = \
# obj/main.o \
# obj/NeighbourList.o \
# obj/mathfuncs.o \
# obj/config.o \
# obj/geom.o \
# obj/error.o \
# obj/spinwaves.o \
# obj/util.o

# CULIBS = \
# obj/cumalloc.o \
# obj/cuthermal.o \
# obj/cuheun.o \
# obj/cuintegrate.o

# obj/%.o: src/%.cpp
# 	$(GCC) $(OPT) -c -o $@ $< 

# obj/%.o: src/%.cu
# 	$(NVCC) -DCUDA -G -c -o $@ $<


# ASD: $(OBJ)
# 	$(GCC) $(OPT) -o $@ $^ $(LIBS)

# clean:
# 	@rm -f obj/*.o

# GCC = g++
# NVCC = nvcc
# OPT = -O3
# LIBS = /cm/shared/apps/fftw/openmpi/gcc/64/3.3.4/lib/libfftw3.a /home/b6033256/libs/libconfig-1.5/lib/libconfig++.a
# CULIBS = -lcurand

OBJ = \
obj/main.o \
obj/NeighbourList.o \
obj/mathfuncs.o \
obj/config.o \
obj/fields.o \
obj/geom.o \
obj/error.o \
obj/spinwaves.o \
obj/util.o

OBJNVCC = \
obj/cumalloc.o \
obj/cuthermal.o \
obj/cufields.o \
obj/cuheun.o \
obj/cufuncs.o

obj/%.o: src/%.cpp
	$(GCC) $(OPT) -DCUDA -c -o $@ $< #-I/cm/shared/apps/fftw/openmpi/gcc/64/3.3.4/include/ -I/home/b6033256/libs/libconfig-1.5/include/

obj/%.o: src/%.cu
	$(NVCC) -O3 -DCUDA -c -o $@ $< #-I/cm/shared/apps/fftw/openmpi/gcc/64/3.3.4/include/


#obj/%.o: src/%.cpp
#       $(NVCC) -DCUDA -G -c -o $@ $< -I/cm/shared/apps/fftw/openmpi/gcc/64/3.3.4/include/ -I/home/b6033256/libs/libconfig-1.5/include/


#ASD: $(OBJ)
#       $(GCC) $(OPT) -o $@ $^ $(LIBS)

ASDcu: $(OBJ) $(OBJNVCC)
	$(NVCC) -DCUDA -G -o $@ $^ $(LIBS) $(CULIBS)


clean:
	@rm -f obj/*.o	

