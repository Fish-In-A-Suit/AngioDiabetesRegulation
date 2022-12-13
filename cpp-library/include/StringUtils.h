// StringUtils.h

#ifndef StringUtils_H
#define StringUtils_H
// alternatively, you can use #pragma once

// class declaration includes are automatically transferred to .cpp implementation
#include <string>
#include <iomanip>
// never use 'using namespace <...>' in .h !!!

// class declaration is made in .h file
class StringUtils{
    int testInt; //just for demonstration

public:
    // StrUtils(int test); // you can define constructor params like so
    StringUtils();
    std::string to_string(double, int);
};

#endif