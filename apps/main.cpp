#include "test-cpp-lib-hello.h"
#include "StringUtils.h"
#include "HighResolutionTimeManager.h"
#include "HighResolutionTimeManagerV2.h"
#include "FileUtils.h"
#include "JsonObject.h"
#include "Logger.h"
#include "PermutationUtils.h"

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

// #include <test-cpp-library/include/test-cpp-lib-hello.h>

using namespace std;

int d_callcount = 0;
int d_perms = 0;

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

// todo: port into OOP style
void predict_miRNA_values(list<string> t)
{
}

long test_sequence_container_speed()
{
    // Use auto keyword to avoid typing long type definitions
    auto start = std::chrono::high_resolution_clock::now();
    int i = 111111111;

    // todo: test vectors, lists, forward lists and deques for speed comparisons
    return -1;
}

// this function is buggy, only works for length = 4 
std::vector<std::string> generatePermutations(std::string sequence, int length) {
    std::vector<std::string> permutations;
    std::string permutation;

    // generate all permutations of the given length
    for (int a = 0; a < sequence.length(); ++a) {
        permutation += sequence[a];
        
        for (int t = 0; t < sequence.length(); ++t) {
            permutation += sequence[t];

            for (int c = 0; c < sequence.length(); ++c) {
                permutation += sequence[c];

                for (int g = 0; g < sequence.length(); ++g) {
                    permutation += sequence[g];

                    // add the permutation to the vector if it is the desired length
                    if(permutation.length() == length) {
                        // std::cout << permutation << std::endl;
                        permutations.push_back(permutation);
                    }
                    permutation.pop_back();
                }
                permutation.pop_back();
            }
            permutation.pop_back();
        }
        permutation.pop_back();
    }
    return permutations;
}

/**
 * Creates permutations of input string S of length n
 * 
 * @param s: Input string, the characters of which will be used in the permutation
 * @param n: The length of the permutations
 * @param count: Pointer to the count variable (should be initialised before this function is called)
 * @param t: Used for recursion, do not set it when calling the function.
 * 
 * Example call: 
 * int count = 0;
 * enumerate("ATCG", 8, count);
*/
void enumerate(const string& s, int n, int &count, string t = "") {
    if (n == 0) {
        ++count;
    } else {
        for (char c : s) {
            ++count;
            enumerate(s, n-1, count, t+c);
        }
    }
}

// Recursive function to generate permutations
void generatePermutationsV1(std::vector<int> permutation, int length)
{
    ++d_callcount;
    // * long long d_size; // enable this to check size of permutation in debug view; setting 'length' to either 6 or 12, or 22 keeps the permutation vector at 24bytes in memory

    // Base case: if the permutation has the desired length, print or store it
    if (permutation.size() == length)
    {
        ++d_perms;
        for (int b : permutation)
        {
            std::cout << std::bitset<2>(b);
        }
        std::cout << std::endl;
        // * d_size = sizeof(permutation);
        return;
    }

    // Recursive case: generate permutations by appending each of the possible values
    for (int b : {0b00, 0b01, 0b10, 0b11})
    {
        std::vector<int> newPermutation = permutation;
        newPermutation.push_back(b);
        generatePermutationsV1(newPermutation, length);
    }
}

/**
 * WARNING: When changing length, make sure to also change std::bitset<VALUE> to length * 2 (C++ doesn't allow it to be dynamic)
*/
void generatePermutationsV2(int permutation, int length) {
    // use for loop
    for (int i = 0; i <= length; ++i) {
        for (int bit : {0b00, 0b01, 0b10, 0b11}) {
            // set the last 2 bits in permutation to 'bit'
            permutation |= bit;

            // TODO: also set the "farther" bits ie after first loop, you need to set bits 3 and 4 from the left 

            // store the bits
            std::cout << std::bitset<12>(permutation) << std::endl;

            // clear the bits - this is necessary, otherwise if 01 is set and then 10, the end will be 11 (should be 10)
            permutation &= 0b00;
        }
        // bit shift for 2 left
        permutation <<= 2;
    }
}

void printVectorElements(std::vector<int> &vec) {
    for (int i = 0; i < vec.size(); i++) {
        std::cout << vec.at(i) << std::endl;
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

    // *** RapidJSON parsing ***
    // // HighResolutionTimeManager hrtm;
    // // cout << "Took " << hrtm.getElapsedTime(Constants::MICROSECONDS) << " to parse file." << endl;
    HighResolutionTimeManagerV2 hrtm2;

    // 65565 is enough memory for all jsons in test_run_1

    // JsonObject causes possible heap corruptions -> check
    JsonObject jsonObj("src_data_files/test.json", 65565, false);
    JsonObject mRNAProductsJson("test_run_1/product_mRNA.json", 65565, false);
    JsonObject productScoresJson("test_run_1/product_scores.json", 6556, false);
    JsonObject termsDirectProductsJson("test_run_1/terms_direct_products.json", 65565, false);
    // // Logger::checkType(&jsonObj);
    // // Logger::debug(termsDirectProductsJson.toString(false));

    printf("book = %s\n", jsonObj.getValue("book"));
    // // printf("characters = %s\n", jsonObj.getValue("characters"));
    // // hrtm2.getElapsedTime(Constants::TimeUnits::MILLISECONDS, true);

    hrtm2.setStartTime();
    // ATTEMPT 1
    // // vector<std::string> permutations = generatePermutations("ATCG", 8);
    // // cout << permutations.size() << endl;

    //ATTEMPT 2: len(12) -> 2,64s   len(14) -> 42,1s 
    // int count = 0;
    // enumerate("ATCG", 4, count);
    // Logger::debug(std::to_string(count));

    // ATTEMPT 3: recursive function; WORKS!
    //generatePermutationsV1({}, 22);
    //std::cout << "The function was called " << d_callcount << " times." << std::endl;
    //std::cout << "Number of permutations: " << d_perms << std::endl;

    // ATTEMPT 4: for-loop function;
    // generatePermutationsV2(0b000000000000, 6);

    // ATTEMPT 5:
    // std::vector<int> permutations = generatePermutationsV3(3);
    std::vector<int> permutations = PermutationUtils::generatePermutations(3, {0b00, 0b01, 0b10, 0b11});
    cout << "Number of permutations: " << permutations.size() << endl;
    printVectorElements(permutations);

    // this causes code to break smh
    // hrtm2.getElapsedTime(Constants::TimeUnits::MILLISECONDS, true);

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