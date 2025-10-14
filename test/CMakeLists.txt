cmake_minimum_required(VERSION 3.22)

# Project settings
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ".")
set(PROJECT_UNIT_TESTS_NAME ${PROJECT_NAME_STR}_unit_tests)
#set(Boost_INCLUDE_DIRS /home/andoresu/Workspace/lib/boost_1_89_0)

# Sources
file(GLOB UNIT_TESTS_SRC_FILES ${PROJECT_SOURCE_DIR}/test/*.cpp)

# Test executable (also links project sources if provided from parent via SRC_FILES)
add_executable(${PROJECT_UNIT_TESTS_NAME}
    ${UNIT_TESTS_SRC_FILES}
    ${SRC_FILES}
)

# Includes
target_include_directories(${PROJECT_UNIT_TESTS_NAME} PRIVATE
    ${PROJECT_INCLUDE_DIR}
    ${Boost_INCLUDE_DIRS}
)

# Include paths if you dont want to expose the path, Boost_INCLUDE_DIRS will be provided by the user at configure time;
# e.g., cmake -S . -B build -DBoost_INCLUDE_DIRS=/home/andoresu/Workspace/lib/boost_1_89_0
#target_include_directories(${PROJECT_UNIT_TESTS_NAME} PRIVATE
#    ${PROJECT_INCLUDE_DIR}
#    ${Boost_INCLUDE_DIRS}
#)

# Enable and register CTest integration
enable_testing()
add_test(NAME MaxDisjointDemoTest COMMAND ${PROJECT_UNIT_TESTS_NAME} --log_level=all)
