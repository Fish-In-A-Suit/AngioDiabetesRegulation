// FileUtils.h
#ifndef FileUtils_H
#define FileUtils_H

#include <string>

class FileUtils {
public:
    FileUtils(int);
    static std::pair<std::string, const char*> getAbsoluteFilepath(std::string);
    static std::string getProjectRootPath();
    static void setProjectRootPath(std::string);
private:
    inline static std::string projectRootPath; //check 'inline' here: https://stackoverflow.com/questions/9110487/undefined-reference-to-a-static-member
};

#endif