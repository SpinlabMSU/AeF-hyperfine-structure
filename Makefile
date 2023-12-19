## AeF-hyprfine-structure makefile
## Intended for use on linux only
AR:=gcc-ar
CXX:=g++
LD:=$(CXX)

HOST_COMPILER ?= g++
NVCC          := $(CUDA_PATH)/bin/nvcc -ccbin $(HOST_COMPILER)

MIN_GCC_VERSION = "13.0"
GCC_VERSION := "`gcc -dumpversion`"
IS_GCC_ABOVE_MIN_VERSION := $(shell expr "$(GCC_VERSION)" ">=" "$(MIN_GCC_VERSION)")
ifeq "$(IS_GCC_ABOVE_MIN_VERSION)" "1"
    # stuff that requires GCC_VERSION >= VERSION
    BUG_FLAGS:=
else
    BUG_FLAGS:=-freport-bug -save-temps 
endif
$(info $$BUG_FLAGS is [${BUG_FLAGS}])
INCLUDES:=-I./include -I./SpinlabHyperfineLib/include
CXX_VSN:=-std=gnu++23 -fmodules-ts -fopenmp
CXXFLAGS:=$(INCLUDES) $(CXX_VSN) -O4 -fPIC -flto $(BUG_FLAGS) -g -march=native -D_GNU_SOURCE
NCXXFLAGS:=$(INCLUDES) -std=gnu++20 -O4 -fPIC -flto $(BUG_FLAGS) -g -march=native
NVCCFLAGS:=$(INCLUDES) -O4 $(BUG_FLAGS) -g $(addprefix -Xcompiler ,$(NCXXFLAGS))
LDFLAGS=-L. -pthread -fopenmp -flto -static-libstdc++ -static-libgcc -g -march=native -O4
CUDA_LIBS:=-lcusolver -lcublas -lcublasLt -lcuda -lcudart_static
LDLIBS:=-l:./libSpinlabHyperfine.a -lgsl -lgslcblas $(CUDA_LIBS) -lz -lm

.PHONY: all clean libs AeF-hyperfine-structure.inl
all: libs aef_hyperfine_structure low_state_dumper stark_diagonalizer
libs: libSpinlabHyperfine.a libSpinlabHyperfine.so


SpinlabHyperfineLib/include/pch.h.gch: SpinlabHyperfineLib/include/pch.h
	$(CXX) -o $@ -x c++-header $(CXXFLAGS) -c $< 
LSPHF_OBJ:=$(patsubst %.cpp,%.o,$(wildcard SpinlabHyperfineLib/src/*.cpp)) $(patsubst %.cu,%.o,$(wildcard SpinlabHyperfineLib/src/*.cu))

libSpinlabHyperfine.a: $(LSPHF_OBJ)
	$(AR) rcs $@ $^
libSpinlabHyperfine.so: libSpinlabHyperfine.a $(LSPHF_OBJ)
	$(CXX) -o $@ -shared $(LDFLAGS) $(CXXFLAGS) $^ -lm -lgsl -lgslcblas $(CUDA_LIBS)
# SpinlabHyperfineLib/include/pch.h.gch
#_MAKEFILE_PROVIDES_DEVFLAG

CUDA_PATH ?= /usr/local/cuda

%.o: %.cu
	$(CXX) $(CXXFLAGS) -o $@ -x c++ -c $<

AeF-hyperfine-structure.inl:
	./prebuild.sh
AeF-hyperfine-structure.o: AeF-hyperfine-structure.cpp AeF-hyperfine-structure.inl
aef_hyperfine_structure: AeF-hyperfine-structure.o AeF-hyperfine-structure.inl libSpinlabHyperfine.so
	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) $< $(LDLIBS)
nodev_aef_hf: AeF-hyperfine-structure.cpp AeF-hyperfine-structure.inl libSpinlabHyperfine.so
	$(CXX) -o $@ $(CXXFLAGS) -D_MAKEFILE_PROVIDES_DEVFLAG $(LDFLAGS) $< $(LDLIBS)
deven_aef_hf: AeF-hyperfine-structure.cpp AeF-hyperfine-structure.inl libSpinlabHyperfine.so
	$(CXX) -o $@ $(CXXFLAGS) -D_MAKEFILE_PROVIDES_DEVFLAG -DENABLE_DEVONSHIRE $(LDFLAGS) $< $(LDLIBS)
low_state_dumper: LowStateDumper/LowStateDumper.o libSpinlabHyperfine.so
	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) $< $(LDLIBS)
stark_diagonalizer: StarkDiagonalizer/StarkDiagonalizer.o libSpinlabHyperfine.so AeF-hyperfine-structure.inl
	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) $< $(LDLIBS)

clean:
	$(RM) aef_hyperfine_structure AeF-hyperfine-structure.inl $(LSPHF_OBJ)\
	SpinlabHyperfineLib/include/pch.h.gch AeF-hyperfine-structure.o libSpinlabHyperfine.*
	StarkDiagonalizer/*.o NoStark_HyperfineTester/*.o LowStateDumper/*.o GenerateHamiltonianFiles/*.o\
	nodev_aef_hf deven_aef_hf low_state_dumper stark_diagonalizer
