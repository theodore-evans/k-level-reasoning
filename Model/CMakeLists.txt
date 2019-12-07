cmake_minimum_required(VERSION 3.10)

project(k-level-reasoning VERSION 1.0)

include_directories(${PROJECT_SOURCE_DIR})
add_executable(k-level-reasoning main.cpp player.cpp config.cpp matrix.cpp) 
