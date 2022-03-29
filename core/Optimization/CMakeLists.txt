cmake_minimum_required(VERSION 3.14)

project(Optimization)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package (Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(.)

add_executable(Main
	       Main.cpp)
target_link_libraries (Main Eigen3::Eigen)
