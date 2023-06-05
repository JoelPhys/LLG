include $(HOSTARGS)

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
obj/thermal.o \
obj/spinwaves.o \
obj/util.o \
obj/heun.o \
obj/main.o

NVCCOBJ = \
obj/cumalloc.o \
obj/cuthermal.o \
obj/cufields.o \
obj/cuheun.o \
obj/cufuncs.o

OBJMP = \
obj/spins.o \
obj/neighbourlist.o \
obj/mathfuncs.o \
obj/config.o \
obj/fields.o \
obj/geom.o \
obj/defects.o \
obj/error.o \
obj/thermal.o \
obj/spinwaves.o \
obj/util.o \
obj/mpheun.o \
obj/main.o


CUDAOBJ=$(OBJECTS:.o=_cuda.o)
MPOBJ=$(OBJMP:.o=_mp.o)

all: ASDcu ASD ASDmp

$(CUDAOBJ): obj/%_cuda.o: src/%.cpp
	$(GCC) $(OPT) -DCUDA -c -o $@ $< $(CONFIGINC) $(FFTW3INC) -DCPUCOMP=$(CPUCOMP) -DGPUCOMP=$(GPUCOMP) -DHOSTNAME=$(HOSTNAME) -DGIT_SHA1=$(GITINFO) -DGITDIRTY=$(GITDIRTY)

$(NVCCOBJ): obj/%.o: src/%.cu
	$(NVCC) $(OPT) -DCUDA -c -o $@ $< $(CONFIGINC) $(FFTW3INC) -DGPUCOMP=$(GPUCOMP) -DHOSTNAME=$(HOSTNAME) -DGIT_SHA1=$(GITINFO) -DGITDIRTY=$(GITDIRTY)

$(MPOBJ): obj/%_mp.o: src/%.cpp
	$(GCC) $(OPT) -fopenmp -DMP -c -o $@ $< $(CONFIGINC) $(FFTW3INC) -DCPUCOMP=$(CPUCOMP) -DGPUCOMP=$(GPUCOMP) -DHOSTNAME=$(HOSTNAME) -DGIT_SHA1=$(GITINFO) -DGITDIRTY=$(GITDIRTY)

$(OBJECTS): obj/%.o: src/%.cpp
	$(GCC) $(OPT) -c -o $@ $< $(CONFIGINC) $(FFTW3INC) -DCPUCOMP=$(CPUCOMP) -DGPUCOMP=$(GPUCOMP) -DHOSTNAME=$(HOSTNAME) -DGIT_SHA1=$(GITINFO) -DGITDIRTY=$(GITDIRTY)

ASDmp: $(MPOBJ)
	$(GCC) -fopenmp -DMP $(OPT) -o $@ $^ $(CONFIGLIB) $(FFTW3LIB) $(CULIBS)

ASDcu: $(CUDAOBJ) $(NVCCOBJ)
	$(NVCC) -DCUDA $(OPT) -o $@ $^ $(CONFIGLIB) $(FFTW3LIB) $(CULIBS)

ASD: $(OBJECTS)
	$(GCC) $(OPT) -o $@ $^ $(CONFIGLIB) $(FFTW3LIB) 

clean:
	@rm -f obj/*.o	

