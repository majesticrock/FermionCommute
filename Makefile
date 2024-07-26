BUILD_DIR = build

all: XP std

build/FermionCommute: $(BUILD_DIR)/Makefile
	@$(MAKE) -C $(BUILD_DIR)

$(BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ ..

XP: build/FermionCommute ../commutators
	./build/FermionCommute XP hubbard
	./build/FermionCommute XP hubbard_dispersions
	./build/FermionCommute XP continuum

std: build/FermionCommute ../commutators
	./build/FermionCommute std hubbard
	./build/FermionCommute std hubbard_dispersions
	./build/FermionCommute std continuum

continuum_xp: build/FermionCommute ../commutators
	./build/FermionCommute XP continuum

hubbard_xp: build/FermionCommute ../commutators
	./build/FermionCommute XP hubbard

test: build/FermionCommute
	./build/FermionCommute test

debug: build/FermionCommute
	./build/FermionCommute debug

../commutators:
	mkdir -p ../commutators

clean:
	rm -rf $(BUILD_DIR)
	rm -rf ../commutators/*
