BUILD_DIR = build

all: build/FermionCommute | ../commutators
	./$(BUILD_DIR)/FermionCommute XP continuum
	./$(BUILD_DIR)/FermionCommute std continuum
	./$(BUILD_DIR)/FermionCommute XP hubbard
	./$(BUILD_DIR)/FermionCommute std hubbard
	./$(BUILD_DIR)/FermionCommute std hubbard_dispersions
	./$(BUILD_DIR)/FermionCommute XP lattice_cut

$(BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ ..

build/FermionCommute: $(BUILD_DIR)/Makefile
	@$(MAKE) -C $(BUILD_DIR)

test debug: $(BUILD_DIR)/Makefile
	@$(MAKE) -C $(BUILD_DIR)
	./$(BUILD_DIR)/FermionCommute $@ hubbard_dispersions

../commutators:
	mkdir -p ../commutators

clean:
	rm -rf $(BUILD_DIR)
	rm -rf ../commutators/*

.PHONY: all clean test debug ../commutators build/FermionCommute
