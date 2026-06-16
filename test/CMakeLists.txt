cmake_minimum_required(VERSION 3.14)

project(fdaPDE-cpp-testing)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_WARNINGS_FLAGS "")
set(CMAKE_CXX_FLAGS "-O2 -march=native -std=c++20 ${CMAKE_WARNINGS_FLAGS}")

if(APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	add_compile_definitions(_Static_assert=static_assert)
endif()

find_package(Eigen3 3.4 REQUIRED NO_MODULE)

include_directories(../src/)
include_directories(../../fdaPDE/)
include_directories(../)
include_directories(../fdaPDE/core/)

# Ipopt
include_directories(/opt/homebrew/opt/ipopt/include/coin-or)
link_directories(/opt/homebrew/opt/ipopt/lib)

include(FetchContent)
FetchContent_Declare(
		googletest
		URL https://github.com/google/googletest/archive/f8d7d77c06936315286eb55f8de22cd23c188571.zip
)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(fdapde_test main.cpp)

target_link_libraries(fdapde_test Eigen3::Eigen)
target_link_libraries(fdapde_test gtest_main)
target_link_libraries(fdapde_test ipopt)

add_test(NAME fdapde_test
		COMMAND $<TARGET_FILE:fdapde_test>
)

include(GoogleTest)
gtest_discover_tests(fdapde_test)
