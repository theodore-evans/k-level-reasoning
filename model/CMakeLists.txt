cmake_minimum_required(VERSION 3.1.3)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

project(k-level-reasoning VERSION 1.0)

# import library target `json`
find_package(PkgConfig REQUIRED)
pkg_check_modules(JSONCPP jsoncpp)
link_libraries(${JSONCPP_LIBRARIES})

include_directories("${JSONCPP_INCLUDE_DIRS}")
include_directories(${PROJECT_SOURCE_DIR})

add_executable(k-level-reasoning.app main.cpp player.cpp config.cpp matrix.cpp)
target_link_libraries(k-level-reasoning.app ${JSONCPP_LIBRARIES})
