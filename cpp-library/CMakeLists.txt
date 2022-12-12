# version 3.11 or later of CMake or needed later for installing GoogleTest
# so let's require it now.
cmake_minimum_required(VERSION 3.11-3.18)

project(cpp-library)

add_library(cpp-library src/test-cpp-lib-hello.cpp)

# PUBLIC needed to make both test-cpp-lib-hello.h and test-cpp-library library available elsewhere in project
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include)

# Tell compiler to use C++20 features. The code doesn't actually use any of them.
target_compile_features(cpp-library PUBLIC cxx_std_20)