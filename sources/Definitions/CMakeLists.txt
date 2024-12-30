set(DEF_SOURCES
    Continuum.cpp
    Hubbard.cpp
    HubbardDispersions.cpp
    DefinitionsBase.cpp
)

add_library(Definitions STATIC ${DEF_SOURCES})

target_include_directories(Definitions PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

find_package(mrock REQUIRED)
target_include_directories(Definitions PUBLIC ${mrock_INCLUDE_DIRS})
target_link_libraries(Definitions PUBLIC ${mrock_LIBRARIES})