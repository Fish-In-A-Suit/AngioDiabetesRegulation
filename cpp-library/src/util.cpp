#include <string>
#include <iomanip>

class StringUtils {
public:
    std::string to_string(double value, int precision) {
        std::ostringstream streamObj; // create an output string stream
        streamObj << std::fixed; // set fixed-point notation
        streamObj << std::setprecision(precision); // set precision
        streamObj << value; // add number to stream
        return streamObj.str(); // get string from output string stream
    }
};