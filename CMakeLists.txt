# FermionCommute/CMakeLists.txt

cmake_minimum_required(VERSION 3.30)

project(FermionCommute VERSION 1.0)

# Custome package install location
list(APPEND CMAKE_PREFIX_PATH "$ENV{HOME}/usr/local")
set(Boost_DIR /usr/lib64/openmpi/lib/cmake/Boost-1.78.0)

# Custom architecture target
set(FERMION_ARCH
    "native"
    CACHE STRING
    "Architecture passed to compiler as -march=<arch> (empty disables)"
)

# Joined compile options
add_library(FermionCommute_options INTERFACE)
target_compile_features(FermionCommute_options
    INTERFACE
        cxx_std_20
)
target_compile_options(FermionCommute_options INTERFACE
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wall>
    $<$<CXX_COMPILER_ID:GNU,Clang>:-Wextra>
)
if(FERMION_ARCH)
    target_compile_options(FermionCommute_options INTERFACE
        $<$<CXX_COMPILER_ID:GNU,Clang>:-march=${FERMION_ARCH}>
    )
endif()


find_package(Boost REQUIRED COMPONENTS serialization)


add_subdirectory(sources/Definitions)
add_executable(FermionCommute 
    sources/FermionCommute.cpp
)

target_link_libraries(FermionCommute
    PRIVATE
        FermionCommute_options
        Definitions
        Boost::serialization
)