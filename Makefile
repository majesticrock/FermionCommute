# Directory for build files
BUILD_DIR = build

# Default target: Build the program and run specified commands
all: $(BUILD_DIR)/Makefile
	@$(MAKE) -C $(BUILD_DIR)
	# Run the program with the required command-line parameters
	./$(BUILD_DIR)/FermionCommute XP continuum
	./$(BUILD_DIR)/FermionCommute std continuum
	./$(BUILD_DIR)/FermionCommute XP hubbard
	./$(BUILD_DIR)/FermionCommute std hubbard
	./$(BUILD_DIR)/FermionCommute std hubbard_dispersions

# Rule to generate the build directory and run CMake
$(BUILD_DIR)/Makefile: CMakeLists.txt
	@mkdir -p $(BUILD_DIR)
	@cd $(BUILD_DIR) && cmake -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ ..

# Test and debug targets
test debug: $(BUILD_DIR)/Makefile
	@$(MAKE) -C $(BUILD_DIR)
	./$(BUILD_DIR)/FermionCommute $@ hubbard_dispersions

# Ensure commutators directory exists
../commutators:
	mkdir -p ../commutators

# Clean up build and commutator directories
clean:
	rm -rf $(BUILD_DIR)
	rm -rf ../commutators/*

# Disable the built-in rule for targets to avoid conflicts
.PHONY: all clean test debug ../commutators
