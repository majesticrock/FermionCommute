BUILD_DIR = build

all: $(BUILD_DIR)/Makefile XP std
	@$(MAKE) -C $(BUILD_DIR)

$(BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ ..

XP: $(BUILD_DIR)/Makefile ../commutators
	@$(MAKE) -C $(BUILD_DIR)
	./$(BUILD_DIR)/FermionCommute XP hubbard
	./$(BUILD_DIR)/FermionCommute XP hubbard_dispersions
	./$(BUILD_DIR)/FermionCommute XP continuum

std: $(BUILD_DIR)/Makefile ../commutators
	@$(MAKE) -C $(BUILD_DIR)
	./$(BUILD_DIR)/FermionCommute std hubbard
	./$(BUILD_DIR)/FermionCommute std hubbard_dispersions
	./$(BUILD_DIR)/FermionCommute std continuum

continuum_xp: $(BUILD_DIR)/Makefile ../commutators
	@$(MAKE) -C $(BUILD_DIR)
	./$(BUILD_DIR)/FermionCommute XP continuum

hubbard_xp: $(BUILD_DIR)/Makefile ../commutators
	@$(MAKE) -C $(BUILD_DIR)
	./$(BUILD_DIR)/FermionCommute XP hubbard

test: $(BUILD_DIR)/Makefile
	@$(MAKE) -C $(BUILD_DIR)
	./$(BUILD_DIR)/FermionCommute test

debug: $(BUILD_DIR)/Makefile
	@$(MAKE) -C $(BUILD_DIR)
	./$(BUILD_DIR)/FermionCommute debug

../commutators:
	mkdir -p ../commutators

clean:
	rm -rf $(BUILD_DIR)
	rm -rf ../commutators/*
