CXX=g++
LD:=$(CXX)
CXXFLAGS:=-I./include -I./SpinlabHyperfineLib/include -std=gnu++23 -O4 -fopenmp -fPIC -flto
LDFLAGS=-L. -pthread -fopenmp -flto -static-libstdc++ -static-libgcc
LDLIBS:=-l:./libSpinlabHyperfine.a -lhwloc -lgsl -lgslcblas -lm -lz

.PHONY: all clean libs AeF-hyperfine-structure.inl
all: libs aef_hyperfine_structure
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

clean:
	$(RM) aef_hyperfine_structure AeF-hyperfine-structure.inl $(LSPHF_OBJ) SpinlabHyperfineLib/include/pch.h.gch AeF-hyperfine-structure.o libSpinlabHyperfine.*
