# version 3.11 or later of CMake or needed later for installing GoogleTest
# so let's require it now.
cmake_minimum_required(VERSION 3.11-3.18)

project(main)

#LINK_DIRECTORIES(C:/Aljosa/Development/Unity-Github/vcpkg/installed/x64-windows/lib)

add_executable(main main.cpp)

# We need test-cpp-lib-hello.h and the hello library
target_link_libraries(main PRIVATE cpp-library ${JSONCPP_LIBRARIES})

# Tell C++ compiler to use C++20 features. We don't actually use any of them.
target_compile_features(main PUBLIC cxx_std_20)