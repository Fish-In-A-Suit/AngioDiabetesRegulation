// #include "test-cpp-lib-hello.h"
#include "StringUtils.h"
#include "HighResolutionTimeManager.h"
#include "HighResolutionTimeManagerV2.h"
#include "FileUtils.h"
#include "JsonObject.h"
#include "Logger.h"
#include "PermutationUtils.h"
#include "SequenceComparator.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <list>
#include <chrono>
#include <json/json.h> //jsoncpp
#include <rapidjson/document.h> //rapidjson
#include "rapidjson/filereadstream.h"
#include <cstdio>
#include <Windows.h>
#include <filesystem>

#include <cuda_runtime.h>

// #include <test-cpp-library/include/test-cpp-lib-hello.h>

using namespace std;


/*
Workflow for miRNA prediction:
  (1) input parameters:
        - productNames: the uniprot ids of products
        - mRNAs: array of mRNA molecules to be analysed (each corresponds to one productName)
        - productScores: the scores of the products
        - miRNALength: the length of the miRNAs that will be computed. Immensely increases required computation power after 14+
        - thresholdToAccept:
        - targetFolder:
        - debugLoopBreak:

  (2) Given miRNALength, generate all possible (unique) mRNA substrings (from all mRNAs) and store them in mRNA_unique_substrings
      Hash the substrings into bits --> mRNA_unique_substrings_hashed. If this is not done, you would have to do complex bit-stepping in part (3)
      (you'd have to step through the mRNA sequence bit by bit and append miRNA at different starting bits, first in a forward direction and then
      backwards too)

  miRNA_mRNA_fitting_dict =
  {
    {
        "miRNA": "xxxxx",
        "mRNAsubstrings": ["xxxxxx", "aaaaaa", "bbbbbb"]
        "mRNAsubstringsScores": [1.0, 0.25, 0.32]
    }
    {
        "miRNA": "xxxxx",
        "mRNAsubstrings": ["xxxxxx", "aaaaaa", "bbbbbb"]
        "mRNAsubstringsScores": [1.0, 0.25, 0.32]
    }
  }

  numberOfPermutations = compute total number of permutations (mathematically)
  int i = 0
  // also compute timing (averageNestedForLoopTime = totalNestedForLoopTime/i) and then find ETA from AverageNestedForLoopTime * (numberOfPermutations - i)
  (3) top level for loop: loop through possible miRNA permutations of ["A", "T", "C", "G"] of given miRNALength
      for each miRNA permutation (in a bit-format):
        matchedSubstrings = []
        matchedScores = [] // the scores of the matchedSubstrings
        for each mRNA substring:
            score = check number of matching nucleotides (if reading from left to right, 2 bits have to be the same consecutively for a single nucleotide match)
            if score > threshold:
                matchedSubstrings.append(mRNA)
                matchedScores.append(score)
        append miRNA, mRNAsubstrings and mRNAsubstringScores into miRNA_mRNA_fitting_dict
        i++
      order & save miRNA_mRNA_fitting_dict
*/


// todo: port into VectorUtils
void printVectorElements(std::vector<int> &vec) {
    for (int i = 0; i < vec.size(); i++) {
        std::cout << vec.at(i) << ", " << std::bitset<6>(vec.at(i)) << std::endl;
    }
} 

int main()
{
    // init functions

    // Warning: if FileUtils is ever needed, make sure to initialise it.
    // this causes INVALID HEAP ARGUMENT: allocated with malloc, freed with operator delete
    // FileUtils fileUtils(0); // directoryClimb = 1, because main.cpp is located in AngioDiabetesRegulation/apps, this climbs up one directory into the project root; directoryClimb no longer required after solving err1
    FileUtils* fileUtils = new FileUtils(0);

    // Logger::setLevel(Constants::LogLevels::DEBUG);

    vector<string> msg{"Hello", "C++", "World", "from", "VS Code", "and the C++ extension!"};
    int i = 0;
    for (const string &word : msg)
    {
        std::cout << word << " ";
        ++i;
    }
    std::cout << endl;

    // get the project source path
    TCHAR szPath[MAX_PATH];
    if (GetModuleFileName(NULL, szPath, MAX_PATH)) {
        std::string path(szPath);
        std::string::size_type pos = path.find_last_of("\\/");
        std::string root = path.substr(0, pos);
        std::cout << "Project root path is: " << root << std::endl;
    }

    // another way of getting the project source path
    // this error is caused by Visual Studio, isn't displayed in VSCode and runs normally in VS also.
    std::cout << "Project root path (using std::filesystem::current_path()): " << std::filesystem::current_path().string() << std::endl;

    // --- RapidJSON parsing ---
    HighResolutionTimeManagerV2 hrtm2;

    // 65565 is enough memory for all jsons in test_run_1
    // JsonObject causes possible heap corruptions -> check
    JsonObject jsonObj("src_data_files/test.json", 65565, false);
    JsonObject mRNAProductsJson("test_run_1/product_mRNA.json", 65565, false);
    JsonObject productScoresJson("test_run_1/product_scores.json", 6556, false);
    JsonObject termsDirectProductsJson("test_run_1/terms_direct_products.json", 65565, false);

    printf("book = %s\n", jsonObj.getValue("book"));
    // --- ---

    // --- Sequence comparisons ---
    SequenceComparator sequenceComparator("src_data_files/miRNAdbs/mirbase_miRNA_hsa-only.txt", "test_run_2/product_mRNAs_cpp.txt");
    // --- ---

    std::vector<int> permutations = PermutationUtils::generatePermutations(3, {0b00, 0b01, 0b10, 0b11});
    cout << "Number of permutations: " << permutations.size() << endl;
    //PermutationUtils::printPermutations(permutations);
    //printVectorElements(permutations);
    // s

    // --- BILLION COUNTING CODE ---
    // hrtm2.setStartTime();
    // takes 3,3 s; 0b(31 x 1)
    // for (uint32_t i; i < 0b1111111111111111111111111111111; i++ ){
    // }

    // takes ~450 years; 0b(63 x 1)
    /*
    for (uint64_t i; i < 0b111111111111111111111111111111111111111111111111111111111111111; i++) {
        if (i % 1000000000 == 0) {
            // print each billion
            std::cout << (i / 1000000000) << " bil" << std::endl;
        }
    }
    */

    /*
    for (uint64_t i = 0; i < 18446744073709551615ull; ++i)
    {
        // Do something with i
    }
    */

    /* enable this for showcasing billion count speeds to Ladi & Umek
    int bilcount = 0;
    for (uint64_t i = 0; i < 274877906943; i++) {
        if (i % 1000000000 == 0)
        {
            bilcount++;
            std::cout << bilcount << " bil." << std::endl;
        }
    }
    */
    // hrtm2.getElapsedTime(Constants::TimeUnits::MILLISECONDS, true);
    // --- END OF BILLION COUNTING CODE ---

    // call destructors
    delete fileUtils;

    return 0;
}

/**
 * A "functional-programming" way of coding to compute elapsed time -> ported to OOP with the HighResolutionTimeManager class.
 *
 * Computes elapsed program runtime by subtracting current time (computed at this function call) from the supplied start time,
 * returning the amount of units (seconds, microseconds, etc.) passed.
 *
 * @param start: the start time computed by std::chrono::high_resolution_clock::now()
 * @param unit: either "nanoseconds", "microseconds", "milliseconds", "seconds", "minutes, "hours" from Constants.TimeUnits
 *
long compute_elapsed_time(std::chrono::high_resolution_clock::time_point start, Constants::TimeUnits timeUnit)
{
    auto stop = std::chrono::high_resolution_clock::now();
    int i = 4;
    auto duration = std::chrono::duration_cast<nanoseconds>(stop - start); // initialise auto var
    // namespaces used for code portability
    // Constants::TimeUnits timeUnit = timeUnit;
    switch (timeUnit)
    {
    case Constants::TimeUnits::NANOSECONDS:
        duration = std::chrono::duration_cast<std::chrono::nanoseconds>(stop - start);
        break;
    case Constants::TimeUnits::MICROSECONDS:
        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        break;
    case Constants::TimeUnits::MILLISECONDS:
        duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
        break;
    case Constants::TimeUnits::SECONDS:
        duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
        break;
    case Constants::TimeUnits::HOURS:
        duration = std::chrono::duration_cast<std::chrono::hours>(stop - start);
        break;
    }
    return duration.count();
}
*/