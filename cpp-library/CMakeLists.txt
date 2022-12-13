# version 3.11 or later of CMake or needed later for installing GoogleTest
# so let's require it now.
cmake_minimum_required(VERSION 3.11-3.18)

project(cpp-library)

add_library(cpp-library 
    src/test-cpp-lib-hello.cpp 
    src/StringUtils.cpp
    src/HighResolutionTimeManager.cpp
    src/Constants.cpp
)
# add_library(cpp-library src/StrUtils.cpp) - this doesnt work, add all libs to the already existing add_library function call.

# PUBLIC needed to make both test-cpp-lib-hello.h and test-cpp-library library available elsewhere in project
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include)

# Tell compiler to use C++20 features. The code doesn't actually use any of them.
target_compile_features(cpp-library PUBLIC cxx_std_20)