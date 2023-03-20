cmake_minimum_required(VERSION 3.11-3.18)

project(cpp-library)

# if there are any errors with JsonObject, also pass "includes/JsonObject.h" here
add_library(cpp-library 
    src/test-cpp-lib-hello.cpp 
    src/StringUtils.cpp
    src/HighResolutionTimeManager.cpp
    src/HighResolutionTimeManagerV2.cpp
    src/Constants.cpp
    src/FileUtils.cpp
    src/JsonObject.cpp
    src/Logger.cpp
    src/PermutationUtils.cpp
    src/SequenceComparator.cpp)

# PUBLIC needed to make both /include/*.h and /src/*.cpp available elsewhere in project
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/includes)

# Tell compiler to use C++20 features. The code doesn't actually use any of them.
target_compile_features(cpp-library PUBLIC cxx_std_20)