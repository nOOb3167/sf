cmake_minimum_required(VERSION 3.12)
cmake_policy(SET CMP0074 NEW) # PackageName_ROOT variables
project(sf)

# Eigen3_ROOT SFML_ROOT assist find_package searches

set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin")

find_package(Eigen3 CONFIG COMPONENTS Eigen3::Eigen REQUIRED)
find_package(SFML CONFIG COMPONENTS audio network graphics system window REQUIRED)

add_executable(sf
	src/selfup/ns_filesys.cpp
	src/selfup/ns_helpers.cpp
	src/tinyxml2/tinyxml2.cpp
	src/tinyxml2/tinyxml2.h
	src/main.cpp
)
target_include_directories(sf PUBLIC ${CMAKE_SOURCE_DIR}/src)
target_include_directories(sf PUBLIC ${Eigen3_ROOT}/../../../include/eigen3)
target_link_libraries(sf sfml-graphics sfml-system sfml-window)
if (WIN32)
	target_link_libraries(sf shlwapi)
endif ()
