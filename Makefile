## AeF-hyprfine-structure makefile
## Intended for use on linux only
AR:=gcc-ar
CXX:=g++
LD:=$(CXX)

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
CXXFLAGS:=$(INCLUDES) $(CXX_VSN) -O4 -fPIC -flto $(BUG_FLAGS)
LDFLAGS=-L. -pthread -fopenmp -flto -static-libstdc++ -static-libgcc
LDLIBS:=-l:./libSpinlabHyperfine.a -lgsl -lgslcblas -lm -lz

.PHONY: all clean libs AeF-hyperfine-structure.inl
all: libs aef_hyperfine_structure low_state_dumper stark_diagonalizer
libs: libSpinlabHyperfine.a libSpinlabHyperfine.so


SpinlabHyperfineLib/include/pch.h.gch: SpinlabHyperfineLib/include/pch.h
	$(CXX) -o $@ -x c++-header $(CXXFLAGS) -c $< 
LSPHF_OBJ:=$(patsubst %.cpp,%.o,$(wildcard SpinlabHyperfineLib/src/*.cpp))

libSpinlabHyperfine.a: $(LSPHF_OBJ)
	$(AR) rcs $@ $^
libSpinlabHyperfine.so: $(LSPHF_OBJ)
	$(CXX) -o $@ -shared $(LDFLAGS) $(CXXFLAGS) $^ -lm -lgsl -lgslcblas
# SpinlabHyperfineLib/include/pch.h.gch



AeF-hyperfine-structure.inl:
	./prebuild.sh
AeF-hyperfine-structure.o: AeF-hyperfine-structure.cpp AeF-hyperfine-structure.inl
aef_hyperfine_structure: AeF-hyperfine-structure.o AeF-hyperfine-structure.inl libSpinlabHyperfine.so
	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) $< $(LDLIBS)
low_state_dumper: LowStateDumper/LowStateDumper.o libSpinlabHyperfine.so
	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) $< $(LDLIBS)
stark_diagonalizer: StarkDiagonalizer/StarkDiagonalizer.o libSpinlabHyperfine.so
	$(CXX) -o $@ $(CXXFLAGS) $(LDFLAGS) $< $(LDLIBS)

clean:
	$(RM) aef_hyperfine_structure AeF-hyperfine-structure.inl $(LSPHF_OBJ) SpinlabHyperfineLib/include/pch.h.gch AeF-hyperfine-structure.o libSpinlabHyperfine.*
