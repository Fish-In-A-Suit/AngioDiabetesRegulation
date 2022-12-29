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
    permutations = bitShiftPrimary(bitPairValues);
    for (int i = 0; i < (length - 2); ++i)
    {
        permutations = bitShiftPrimary(permutations);
    }
    return permutations;
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

