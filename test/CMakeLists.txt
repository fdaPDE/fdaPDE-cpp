cmake_minimum_required(VERSION 3.14)

project(TEST)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_FLAGS "-O0 -march=native -std=c++17")

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(.)
include_directories(../src/core/)

# enable testing functionality using gtest
include(FetchContent)
FetchContent_Declare(
  googletest
  URL https://github.com/google/googletest/archive/609281088cfefc76f9d0ce82e1ff6c30cc3591e5.zip
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(fdaPDE_test MainTest.cpp)
target_link_libraries (fdaPDE_test Eigen3::Eigen)
target_link_libraries (fdaPDE_test gtest_main)

add_test(NAME fdaPDE_test
	 COMMAND $<TARGET_FILE:fdaPDE_test>
  )

include(GoogleTest)
gtest_discover_tests(fdaPDE_test)
