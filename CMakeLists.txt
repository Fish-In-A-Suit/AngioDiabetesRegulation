cmake_minimum_required(VERSION 3.11-3.18)

# set(CMAKE_FIND_DEBUG_MODE 1)

project(
    AngioDiabetesRegulation
    VERSION 0.1
    DESCRIPTION "Small CMake example built with VSCode")

# need to have include directories so vscode can find vpckg libraries
include_directories(C:/Aljosa/Development/Unity-Github/vcpkg/installed/x64-windows/include)

# find_package checks in C:\Program Files\CMake\share\cmake-3.25\Modules and tried to find FindRapidjson.cmake
# (I snatched a custom rapidjson.cmake that sets Rapidjson_INCLUDE_DIRS variable)
# Rapidjson_INCLUDE_DIRS is set to C:/Aljosa/Development/Unity-Github/vcpkg/installed/x64-windows/include
# find_package(Rapidjson CONFIG REQUIRED)
# include_directories(${Rapidjson_INCLUDE_DIRS})  
# message("rapidjson include dirs = ${Rapidjson_INCLUDE_DIRS}")

add_subdirectory(cpp-library)   # look in "test-cpp-library" subdirectory for CMakeLists.txt to process
add_subdirectory(apps)    # look in "test-apps" subdirectory for CMakeLists.txt to process