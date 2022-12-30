// PermutationUtils.cpp

#include "PermutationUtils.h"


PermutationUtils::~PermutationUtils() {
    // destructor implementation
    std::cout << "PermutationUtils destructor called." << std::endl;
}

// todo: remove this, doesnt work
std::vector<int> PermutationUtils::generatePermutations(int length, std::vector<int> bitPairValues = {0b00, 0b01, 0b10, 0b11})
{
    std::vector<int> permutations = {};
    permutations = bitShiftPrimary(bitPairValues, 2, {0b00, 0b01, 0b10, 0b11});
    for (int i = 0; i < (length - 2); ++i)
    {
        permutations = bitShiftPrimary(permutations, 2, {0b00,0b01,0b10,0b11});
    }
    return permutations;
}

/**
 * Example call: generatePermutationsRecursively({}, 12);
 */
void PermutationUtils::generatePermutationsRecursively(std::vector<int> permutation, int length) {
    // Base case: if the permutation has the desired length, print or store it
    if (permutation.size() == length)
    {
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
        generatePermutationsRecursively(newPermutation, length);
    }
}

/**
 * The length of permutation can be seen nucleotide-wise or bit-wise (which is 2 times the length of nucleotides).
 * This function excepts the permutationLength to be specified bit-wise (ie. AAAA -> len = 8 bits)
 * 
 * WARNING: Make sure to adjust the value in std::bitset<VALUE> according to the amount of bits you want to represent.
*/
void PermutationUtils::printPermutations(std::vector<int> permutations) {
    for (int i = 0; i < permutations.size(); i++)
    {
        std::cout << permutations.at(i) << ", " << std::bitset<6>(permutations.at(i)) << std::endl;
    }
}

std::vector<int> PermutationUtils::bitShiftPrimary(int permutation, int bitShift = 2, std::vector<int> primaryBitPairValues = {0b00, 0b01, 0b10, 0b11})
{
    std::vector<int> result = {};

    permutation <<= bitShift; // bit shift original permutation by two (to preserve the original bits)
    for (int bitPair : primaryBitPairValues)
    {
        permutation |= bitPair; // set the last bit pair
        result.push_back(permutation);

        // permutation &= 0b00; // clear the last bit pair; this is wrong

        // clear the last 2 bits
        permutation &= ~(1 << 0);
        permutation &= ~(1 << 1);
    }
    return result;
}

std::vector<int> PermutationUtils::bitShiftPrimary(std::vector<int> permutations, int bitShift = 2, std::vector<int> primaryBitPairValues = {0b00, 0b01, 0b10, 0b11}) 
{
    std::vector<int> result = {};
    for (int permutation : permutations)
    {
        std::vector<int> shiftedPermutations = bitShiftPrimary(permutation, bitShift, primaryBitPairValues);
        result.insert(result.end(), shiftedPermutations.begin(), shiftedPermutations.end());
    }
    return result;
}

