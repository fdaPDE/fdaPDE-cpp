cmake_minimum_required(VERSION 3.14)

project(fdaPDE-cpp-testing)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_WARNINGS_FLAGS "")
set(CMAKE_CXX_FLAGS "-O2 -march=native -std=c++20 ${CMAKE_WARNINGS_FLAGS}")

option(FDAPDE_ENABLE_COUT "Enable fdapde::cout diagnostic output" OFF)
if(FDAPDE_ENABLE_COUT)
	add_compile_definitions(FDAPDE_ENABLE_COUT)
endif()

if(APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	add_compile_definitions(_Static_assert=static_assert)
endif()

find_package(Eigen3 3.4 REQUIRED NO_MODULE)

include_directories(../src/)
include_directories(../../fdaPDE/)
include_directories(../)
include_directories(../fdaPDE/core/)

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

if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/frgcca_covmax_timing.cpp")
    add_executable(frgcca_covmax_timing frgcca_covmax_timing.cpp)
    target_compile_definitions(frgcca_covmax_timing PRIVATE FDAPDE_ENABLE_COUT)
    target_link_libraries(frgcca_covmax_timing Eigen3::Eigen)
endif()

add_executable(rgcca_bootstrap_reproducibility rgcca_bootstrap_reproducibility.cpp)
target_link_libraries(rgcca_bootstrap_reproducibility Eigen3::Eigen)

add_test(NAME fdapde_test
		COMMAND $<TARGET_FILE:fdapde_test>
)

add_test(NAME rgcca_bootstrap_worker_reproducibility
    COMMAND ${CMAKE_COMMAND}
        -DREPRO_EXECUTABLE=$<TARGET_FILE:rgcca_bootstrap_reproducibility>
        -P ${CMAKE_CURRENT_SOURCE_DIR}/compare_rgcca_bootstrap_reproducibility.cmake.in
)

include(GoogleTest)
gtest_discover_tests(fdapde_test)
