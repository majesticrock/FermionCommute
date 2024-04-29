CXX = g++

# -isystem rather than -I as this supresses warnings that occur during
# the compilation of the eigen library (I cant fix them anyways)
INCLUDEFLAGS = -isystem ~/usr/local/include

CXXFLAGS = -g $(WARNINGS) -std=c++17 $(OPT) -fopenmp

LDLIBS = -L/sw/gcc/5.3.0/rtf/lib64 -L/home/joshua/usr/local/include/boost_lib/ -lboost_serialization

WARNINGS = -Wall -Wno-sign-compare

OPT = -march=native -O3

DEF_SRCS=Definitions/Hubbard.cpp Definitions/Continuum.cpp Definitions/StandardOperators.cpp
PART_SRCS=Coefficient.cpp IndexWrapper.cpp Momentum.cpp MomentumList.cpp Operator.cpp Term.cpp WickCleaner.cpp WickOperator.cpp WickOperatorTemplate.cpp WickTerm.cpp 
SRCS=$(PART_SRCS) $(DEF_SRCS) FermionCommute.cpp

OBJS=$(addprefix build/, $(subst .cpp,.o,$(SRCS)))

all: build build/main ../commutators XP std

XP: build build/main ../commutators
	./build/main XP hubbard

std: build build/main ../commutators
	./build/main std hubbard

test: build build/main
	./build/main test

debug: build build/main
	./build/main debug

build/main: $(OBJS) | build
	$(CXX) $(INCLUDEFLAGS) -o build/main $(OBJS) $(CXXFLAGS) $(LDLIBS)

build/%.o: sources/%.cpp
	$(CXX) $(INCLUDEFLAGS) $< -o $@ -c $(CXXFLAGS)

build:
	mkdir -p build
	mkdir -p build/Definitions

../commutators:
	mkdir -p ../commutators

clean:
	rm -rf build
	rm -rf ../commutators/*
