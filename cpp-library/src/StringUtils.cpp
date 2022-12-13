// StringUtils.cpp

#include "StringUtils.h"

// class implementation belongs to the .cpp file

// you can implement constructor params like so:
// StringUtils::StringUtils(int test) {
//     testInt = test;
// }

// empy constructor implementation
// note: no ';' after class brackets!
StringUtils::StringUtils() {

}

std::string StringUtils::to_string(double value, int precision) {
    std::ostringstream streamObj;              // create an output string stream
    streamObj << std::fixed;                   // set fixed-point notation
    streamObj << std::setprecision(precision); // set precision
    streamObj << value;                        // add number to stream
    return streamObj.str();                    // get string from output string stream
}

