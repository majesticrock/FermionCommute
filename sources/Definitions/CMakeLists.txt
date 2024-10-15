set(DEF_SOURCES
    Continuum.cpp
    Hubbard.cpp
    HubbardDispersions.cpp
    DefinitionsBase.cpp
)

add_library(Definitions STATIC ${DEF_SOURCES})

target_include_directories(Definitions PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

find_package(SymbolicOperators REQUIRED)
target_include_directories(Definitions PUBLIC ${SymbolicOperators_INCLUDE_DIRS})
target_link_libraries(Definitions PUBLIC SymbolicOperators::SymbolicOperators)

find_package(Utility REQUIRED)
target_include_directories(Definitions PUBLIC ${Utility_INCLUDE_DIRS})
target_link_libraries(Definitions PUBLIC Utility::Utility)
