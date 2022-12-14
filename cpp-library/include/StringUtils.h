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

public:
    // StrUtils(int test); // you can define constructor params like so
    static std::string to_string(double, int); // declare statics here and not in cpp; use StringUtils::to_string to call.

private:
    StringUtils(); // constructor declared private, as StringUtils is meant to provide static string manipulation methods
};

#endif