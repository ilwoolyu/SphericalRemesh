project(SphericalRemesh)
cmake_minimum_required(VERSION 2.8)

find_package(SlicerExecutionModel REQUIRED)

include(${SlicerExecutionModel_USE_FILE})
include(${GenerateCLP_USE_FILE})
include_directories( ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} )

generateCLP(SRemeshCLP ${CMAKE_SOURCE_DIR}/wrapper/SRemesh.xml)

# MESH library
find_package(MeshLib REQUIRED)
include_directories(${Mesh_INCLUDE_DIRS} ${Mesh_INCLUDE_DIRS}/Util)
link_directories(${Mesh_LIBRARY_DIRS})

include_directories(src)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY .)
set (Remesh_SOURCES ${CMAKE_SOURCE_DIR}/src/SphericalRemeshing.cpp)
add_executable(SphericalRemesh ${Remesh_SOURCES} wrapper/main.cpp ${SRemeshCLP})
target_link_libraries(SphericalRemesh Mesh)
