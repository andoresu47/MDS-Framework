cmake_minimum_required(VERSION 3.22)

# Definitions
set(PROJECT_NAME_STR mds_framework)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/bin")

# Project title and define versioning
project(${PROJECT_NAME_STR})

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Load machine-specific paths (e.g. Boost_INCLUDE_DIRS) from ./cmake/local.cmake. 
if(EXISTS "${CMAKE_SOURCE_DIR}/cmake/local.cmake")
  include("${CMAKE_SOURCE_DIR}/cmake/local.cmake")
endif()

# Project includes 
set(PROJECT_INCLUDE_DIR ${PROJECT_INCLUDE_DIR} ${PROJECT_SOURCE_DIR}/include)

# Sources: grab everything under src, then exclude main.cpp
file(GLOB_RECURSE SRC_FILES "${PROJECT_SOURCE_DIR}/src/*/*.cpp")

set(MAIN_FILE "${PROJECT_SOURCE_DIR}/src/main.cpp")
set(SOURCE_FILES ${MAIN_FILE} ${SRC_FILES})

# Executable
add_executable(${PROJECT_NAME_STR} ${SOURCE_FILES})
target_include_directories(${PROJECT_NAME_STR} PRIVATE 
    "${PROJECT_INCLUDE_DIR}"
)
if(DEFINED Boost_INCLUDE_DIRS)
  target_include_directories(${PROJECT_NAME_STR} PRIVATE "${Boost_INCLUDE_DIRS}")
endif()

# Unit tests subdir
add_subdirectory(test)
