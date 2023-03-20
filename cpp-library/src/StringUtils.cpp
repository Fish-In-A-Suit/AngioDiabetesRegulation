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

// destructor implementation
StringUtils::~StringUtils() {
    std::cout << "StringUtils destructor called" << std::endl;
}

std::string StringUtils::to_string(double value, int precision) {
    std::ostringstream streamObj;              // create an output string stream
    streamObj << std::fixed;                   // set fixed-point notation
    streamObj << std::setprecision(precision); // set precision
    streamObj << value;                        // add number to stream
    return streamObj.str();                    // get string from output string stream
}

/*
 * Splits 'str' string based on delimiter and populates the remaining split elements inside the supplied
 * std::vector<std::string>
 * 
 * @param str: The string to split
 * @param delimiter: The delimiter string where the splits occur
 * @param line_elements: an empty std::vector<std::string> to populate
 * 
 * @return std::vector<std::string> of remaining elements after split
 */
void StringUtils::split(std::string& str, std::string delimiter, std::vector<std::string>& line_elements) {
    size_t pos = 0;
    std::string token;

    while ((pos = str.find(delimiter)) != std::string::npos) {
        token = str.substr(0, pos);
        line_elements.push_back(token);
        str.erase(0, pos + delimiter.length());
    }
    return;
}

bool StringUtils::contains(std::string& str, std::string& substr) {
    std::cout << "str = " << str << ", substr = " << substr << ", str.find(substr) = " << str.find(substr) << std::endl;
    if(str.find(substr) != std::string::npos) {
        std::cout << "true" << std::endl;
        return true; // found
    } else {
        std::cout << "false" << std::endl;
        return false; // not found
    }
}

bool StringUtils::contains(std::string& str, const char* substr) {
    std::string substr_str = substr; // std::string has a constructor from const char *, so this is legal (https://stackoverflow.com/questions/24127946/converting-a-const-char-to-stdstring)
    contains(str, substr_str);
}

/*
 * Prints the std::vector<std::string> string elements line by line to std::cout. Pass-by-reference is used (&vecString), since the contents of
 * std::vector are only read. 
 * 
 * @param vecString: A vector of strings whose elements to print to std::cout
 */
void StringUtils::print_vector(std::vector<std::string> &vecString) {
    for (std::string str : vecString) {
        std::cout << str << std::endl;
    }
}



