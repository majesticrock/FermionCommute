# cmake/fermion_commute_options.cmake

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

target_include_directories(FermionCommute_options INTERFACE EXTRA_INCLUDE_DIRS)
