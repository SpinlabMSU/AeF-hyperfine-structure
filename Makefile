## AeF-hyperfine-structure makefile
## Intended for use on linux only
AR:=gcc-ar
CXX:=g++
LD:=$(CXX)

HOST_COMPILER ?= g++

INCLUDES:=-I./include -I./SpinlabHyperfineLib/include
ROOT_CFLAGS:=$(shell root-config --cflags | sed 's/-std=.* //g') # note: want a more recent C++ version than root requires
ROOT_LDFLAGS:=$(shell root-config --ldflags)
ROOT_LDLIBS:=$(shell root-config --glibs)
CXX_VSN:=-std=gnu++23 -fmodules-ts -fopenmp
CXXFLAGS:=$(INCLUDES) $(CXX_VSN) -O4 -fPIC -flto $(BUG_FLAGS) -g -march=native -D_GNU_SOURCE
LDFLAGS=-L. -pthread -fopenmp -flto -static-libstdc++ -static-libgcc -g -march=native -O4
CUDA_LIBS:=-lcusolver -lcublas -lcublasLt -lcuda -lcudart_static
LDLIBS:=-l:./libSpinlabHyperfine.a -lgsl -lgslcblas $(CUDA_LIBS) -lz -lm

.PHONY: all clean libs AeF-hyperfine-structure.inl
all: libs aef_hyperfine_structure low_state_dumper stark_diagonalizer nodev_aef_hf deven_aef_hf\
operatorVisualizer perturbation_analyzer
libs: libSpinlabHyperfine.a libSpinlabHyperfine.so


SpinlabHyperfineLib/include/pch.h.gch: SpinlabHyperfineLib/include/pch.h
	$(CXX) -o $@ -x c++-header $(CXXFLAGS) -c $< 
LSPHF_OBJ:=$(patsubst %.cpp,%.o,$(wildcard SpinlabHyperfineLib/src/*.cpp))\
$(patsubst %.cu,%.o,$(wildcard SpinlabHyperfineLib/src/*.cu))\
$(patsubst %.cpp,%.o,$(wildcard SpinlabHyperfineLib/src/backends/*.cpp))\
$(patsubst %.cpp,%.o,$(wildcard SpinlabHyperfineLib/src/operators/*.cpp))

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
	$(CXX) -o $@ $(CXXFLAGS) -D_MAKEFILE_PROVIDES_DEVFLAG -DUSE_DEVONSHIRE $(LDFLAGS) $< $(LDLIBS)
nodev_cpu_hf: AeF-hyperfine-structure.cpp AeF-hyperfine-structure.inl libSpinlabHyperfine.so
	$(CXX) -o $@ $(CXXFLAGS) -D_MAKEFILE_PROVIDES_DEVFLAG -DDONT_USE_CUDA $(LDFLAGS) $< $(LDLIBS)
deven_cpu_hf: AeF-hyperfine-structure.cpp AeF-hyperfine-structure.inl libSpinlabHyperfine.so
	$(CXX) -o $@ $(CXXFLAGS) -D_MAKEFILE_PROVIDES_DEVFLAG -DUSE_DEVONSHIRE -DDONT_USE_CUDA $(LDFLAGS) $< $(LDLIBS)
low_state_dumper: LowStateDumper/LowStateDumper.o libSpinlabHyperfine.so
	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) $< $(LDLIBS)
stark_diagonalizer: StarkDiagonalizer/StarkDiagonalizer.o libSpinlabHyperfine.so AeF-hyperfine-structure.inl
	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) $< $(LDLIBS)

operatorVisualizer: operator_visualizer/operator_visualizer.cpp libSpinlabHyperfine.so AeF-hyperfine-structure.inl
	$(CXX) -o $@ $(CXXFLAGS) $(ROOT_CFLAGS) -Ioperator_visualizer/include $(LDFLAGS) $(ROOT_LDFLAGS) \
	$< $(ROOT_LDLIBS) $(LDLIBS)

perturbation_analyzer: PerturbationAnalyzer/PerturbationAnalyzer.o libSpinlabHyperfine.so
	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) $< $(LDLIBS)

clean:
	$(RM) aef_hyperfine_structure AeF-hyperfine-structure.inl $(LSPHF_OBJ) \
	SpinlabHyperfineLib/include/pch.h.gch AeF-hyperfine-structure.o libSpinlabHyperfine.* \
	StarkDiagonalizer/StarkDiagonalizer.o NoStark_HyperfineTester/NoStark_HyperfineTester.o \
	LowStateDumper/LowStateDumper.o GenerateHamiltonianFiles/GenerateHamiltonianFiles.o \
	nodev_aef_hf deven_aef_hf low_state_dumper stark_diagonalizer
