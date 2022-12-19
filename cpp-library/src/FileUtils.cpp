// FileUtils.cpp

#include "FileUtils.h"


/**
 * FileUtils constructor. It initialises the projectRootPath to the folder which contains the .cpp file from where FileUtils(int) constr is called.
 * Uses directoryClimb to climb into the "true" project directory in case if the calling .cpp file is nested.
*/
FileUtils::FileUtils(int directoryClimb) {
    // find and set project root path
    FileUtils::projectRootPath = std::filesystem::current_path().string();
    std::cout << "project root path before climbing: " << projectRootPath << std::endl;
    for(int i = 0; i < directoryClimb; i++) {
        std::filesystem::path p = FileUtils::projectRootPath; 
        FileUtils::projectRootPath = p.parent_path().string();
    }
    std::cout << "project root path = " << projectRootPath << std::endl;
}

/**
 * Returns the absolute filepath to the relative filepath parameter.
 * @param relativeProjectFilepath: a relative path to a directory or a file on the project 
 * @returns absolute filepath to the specified relative file
 * 
 * Read: https://stackoverflow.com/questions/41872550/returning-two-values-of-two-different-types-in-a-c-function
*/
std::pair<std::string, const char*> FileUtils::getAbsoluteFilepath(std::string relativeProjectFilepath) {
    std::filesystem::path root = FileUtils::projectRootPath;
    std::filesystem::path rel = relativeProjectFilepath;
    std::filesystem::path result = root / rel;
    std::string resultString = result.string(); //TODO: this destroy double slashes --> search how to return str without it destroying double slashes
    std::cout << "Returning " << resultString << " root = " << root << " rel = " << rel << std::endl;
    return std::make_pair(resultString, resultString.c_str());
    //return result.string();
}

/**
 * @returns The root path of the project (if set beofre either via setProjectRootPath or using the constructor)
*/
std::string FileUtils::getProjectRootPath() {
    return FileUtils::projectRootPath;
}

void FileUtils::setProjectRootPath(std::string rootPath) {
    FileUtils::projectRootPath = rootPath;
}