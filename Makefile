include hostargs/Make.inc.local

#########################################################################################
# Options are SERIAL, OMP, OMPGPU or GPU
#some compiler definitions to sort out
CPUCOMP='"GNU C++ Compiler $(shell g++ --version | head -n 1 | cut -b 5-)"'
GPUCOMP='"NVIDIA C++ Compiler $(shell nvcc --version | tail -n 2 | head -n 1)"'
GITINFO='"$(shell git rev-parse HEAD)"'
GITDIRTY='"$(shell git status -s | grep -v ? | grep -e 'src\/' -e 'inc\/' | wc -l)"'
HOSTNAME='"$(shell hostname)"'
#########################################################################################


OBJECTS = \
obj/spins.o \
obj/neighbourlist.o \
obj/mathfuncs.o \
obj/config.o \
obj/fields.o \
obj/geom.o \
obj/defects.o \
obj/error.o \
obj/spinwaves.o \
obj/util.o \
obj/main.o

OBJNVCC = \
obj/cumalloc.o \
obj/cuthermal.o \
obj/cufields.o \
obj/cuheun.o \
obj/cufuncs.o

CUDAOBJ=$(OBJECTS:.o=_cuda.o)

all: ASDcu ASD

ASDcu: $(CUDAOBJ) $(OBJNVCC)
	$(NVCC) -DCUDA $(OPT) -o $@ $^ $(CONFIGLIB) $(FFTW3LIB) $(CULIBS)

$(CUDAOBJ): obj/%_cuda.o: src/%.cpp
	$(GCC) $(OPT) -DCUDA -c -o $@ $< $(CONFIGINC) $(FFTW3INC) -DCPUCOMP=$(CPUCOMP) -DGPUCOMP=$(GPUCOMP) -DHOSTNAME=$(HOSTNAME) -DGIT_SHA1=$(GITINFO) -DGITDIRTY=$(GITDIRTY)

obj/%.o: src/%.cu
	$(NVCC) $(OPT) -DCUDA -c -o $@ $< $(CONFIGINC) $(FFTW3INC) -DGPUCOMP=$(GPUCOMP) -DHOSTNAME=$(HOSTNAME) -DGIT_SHA1=$(GITINFO) -DGITDIRTY=$(GITDIRTY)


ASD: $(OBJECTS)
	$(GCC) $(OPT) -o $@ $^ $(CONFIGLIB) $(FFTW3LIB) 

$(OBJECTS): obj/%.o: src/%.cpp
	$(GCC) $(OPT) -c -o $@ $< $(CONFIGINC) $(FFTW3INC) -DCPUCOMP=$(CPUCOMP) -DGPUCOMP=$(GPUCOMP) -DHOSTNAME=$(HOSTNAME) -DGIT_SHA1=$(GITINFO) -DGITDIRTY=$(GITDIRTY)


clean:
	@rm -f obj/*.o	

