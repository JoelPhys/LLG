include MakeInclude/Make.inc.local

#########################################################################################
# Options are SERIAL, OMP, OMPGPU or GPU
#some compiler definitions to sort out
CPUCOMP='"GNU C++ Compiler $(shell g++ --version | head -n 1 | cut -b 5-)"'
GPUCOMP='"NVIDIA C++ Compiler $(shell nvcc --version | tail -n 2 | head -n 1)"'
GITINFO='"$(shell git rev-parse HEAD)"'
GITDIRTY='"$(shell git status -s | grep -v ? | grep -e 'src\/' -e 'inc\/' | wc -l)"'
HOSTNAME='"$(shell hostname)"'
#########################################################################################


OBJ = \
obj/spins.o \
obj/neighbourlist.o \
obj/mathfuncs.o \
obj/config.o \
obj/fields.o \
obj/geom.o \
obj/error.o \
obj/spinwaves.o \
obj/util.o \
obj/heun.o \
obj/main.o

OBJNVCC = \
obj/cumalloc.o \
obj/cuthermal.o \
obj/cufields.o \
obj/cuheun.o \
obj/cufuncs.o

obj/%.o: src/%.cpp
	$(GCC) $(OPT) -DCUDA -c -o $@ $< -DCPUCOMP=$(CPUCOMP) -DGPUCOMP=$(GPUCOMP) -DHOSTNAME=$(HOSTNAME) -DGIT_SHA1=$(GITINFO) -DGITDIRTY=$(GITDIRTY)

obj/%.o: src/%.cu
	$(NVCC) -O3 -DCUDA -G -c -o $@ $< -DGPUCOMP=$(GPUCOMP) -DHOSTNAME=$(HOSTNAME) -DGIT_SHA1=$(GITINFO) -DGITDIRTY=$(GITDIRTY)


# ASD: $(OBJ)
# 	$(GCC) $(OPT) -o $@ $^ $(LIBS)

# clean:
# 	@rm -f obj/*.o

# GCC = g++
# NVCC = nvcc
# OPT = -O3
# LIBS = /cm/shared/apps/fftw/openmpi/gcc/64/3.3.4/lib/libfftw3.a /home/b6033256/libs/libconfig-1.5/lib/libconfig++.a
# CULIBS = -lcurand

# OBJ = \
# obj/spins.o \
# obj/neighbourlist.o \
# obj/mathfuncs.o \
# obj/config.o \
# obj/fields.o \
# obj/geom.o \
# obj/error.o \
# obj/spinwaves.o \
# obj/util.o \
# obj/heun.o \
# obj/main.o 

# OBJNVCC = \
# obj/cumalloc.o \
# obj/cuthermal.o \
# obj/cufields.o \
# obj/cuheun.o \
# obj/cufuncs.o

# obj/%.o: src/%.cpp
# 	$(GCC) $(OPT) -DCUDA -c -o $@ $< #-I/cm/shared/apps/fftw/openmpi/gcc/64/3.3.4/include/ -I/home/b6033256/libs/libconfig-1.5/include/

# obj/%.o: src/%.cu
# 	$(NVCC) -O3 -DCUDA -c -o $@ $< #-I/cm/shared/apps/fftw/openmpi/gcc/64/3.3.4/include/


#obj/%.o: src/%.cpp
#       $(NVCC) -DCUDA -G -c -o $@ $< -I/cm/shared/apps/fftw/openmpi/gcc/64/3.3.4/include/ -I/home/b6033256/libs/libconfig-1.5/include/


# ASD: $(OBJ)
# 	$(GCC) $(OPT) -o $@ $^ $(LIBS)

ASDcu: $(OBJ) $(OBJNVCC)
	$(NVCC) -DCUDA -O3 -G -o $@ $^ $(CONFIGLIB) $(FFTW3LIB) $(CULIBS)


clean:
	@rm -f obj/*.o	

