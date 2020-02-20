project(SphericalRemesh)
cmake_minimum_required(VERSION 3.5)

set (CMAKE_CXX_STANDARD 11)

set(CLI11_VER "1.9.0" CACHE STRING "CLI11 Version: see https://github.com/CLIUtils/CLI11/releases")
mark_as_advanced(CLI11_VER)
set(CLI11_DIRS ${CMAKE_BINARY_DIR}/CMakeFiles/CLI11)
file(DOWNLOAD https://github.com/CLIUtils/CLI11/releases/download/v${CLI11_VER}/CLI11.hpp ${CLI11_DIRS}/CLI11.hpp)

if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	message(STATUS "Setting build type to 'Release' as none was specified.")
	set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
	mark_as_advanced(CMAKE_BUILD_TYPE)
	# Set the possible values of build type for cmake-gui
	set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

include_directories( ${CMAKE_SOURCE_DIR} ${CMAKE_BINARY_DIR} )

include_directories(${CLI11_DIRS})

# MESH library
find_package(MeshLib REQUIRED)
include_directories(${Mesh_INCLUDE_DIRS} ${Mesh_INCLUDE_DIRS}/Util)
link_directories(${Mesh_LIBRARY_DIRS})

include_directories(src)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY .)
set (Remesh_SOURCES ${CMAKE_SOURCE_DIR}/src/SphericalRemeshing.cpp)
add_executable(SphericalRemesh ${Remesh_SOURCES} wrapper/main.cpp ${SRemeshCLP})
target_link_libraries(SphericalRemesh Mesh)
