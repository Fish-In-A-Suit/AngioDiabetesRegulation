// StringUtils.h

#ifndef StringUtils_H
#define StringUtils_H
// alternatively, you can use #pragma once

// class declaration includes are automatically transferred to .cpp implementation
#include <string>
#include <iomanip>
#include <iostream>
#include <vector>
// never use 'using namespace <...>' in .h !!!

// class declaration is made in .h file
class StringUtils{

public:
    // StrUtils(int test); // you can define constructor params like so
    ~StringUtils();
    static std::string to_string(double, int); // declare statics here and not in cpp; use StringUtils::to_string to call.
    static void split(std::string&, std::string, std::vector<std::string>&);
    static void print_vector(std::vector<std::string>&);
    static bool contains(std::string &, std::string &);
    static bool contains(std::string &, const char*);

private:
    StringUtils(); // constructor declared private, as StringUtils is meant to provide static string manipulation methods
};

#endif