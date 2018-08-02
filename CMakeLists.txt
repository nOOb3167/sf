cmake_minimum_required(VERSION 3.10)
project(sf)

set(EXECUTABLE_OUTPUT_PATH "${CMAKE_SOURCE_DIR}/bin")

set(SF_MODULE_PATH_SFML -NOTFOUND CACHE STRING "path")
set(SF_MODULE_PATH_EIGEN -NOTFOUND CACHE STRING "path")

list(APPEND CMAKE_MODULE_PATH ${SF_MODULE_PATH_SFML} ${SF_MODULE_PATH_EIGEN})

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
target_include_directories(sf PUBLIC ${SF_MODULE_PATH_EIGEN}/../../../include/eigen3)
target_link_libraries(sf sfml-graphics sfml-system sfml-window)
if (WIN32)
	target_link_libraries(sf shlwapi)
endif ()
