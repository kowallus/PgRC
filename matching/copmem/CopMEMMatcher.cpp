/* ================================================================= *
*  CopMEMMatcher is  based on copMEM:                               *
*                                                                   *
*  copMEM is a program for efficient computation of MEMs            *
*  (Maximal Exact Matches) in a pair of genomes.                    *
*  Its main algorithmic idea requires that two internal parameters  *
*  (k1 and k2) are coprime, hence the name.                         *
*                                                                   *
*                                                                   *
*  Copyright (c) 2018, Szymon Grabowski and Wojciech Bieniecki      *
*  All rights reserved                                              *
*                                                                   *
*  This program is free software: you can redistribute it and/or    *
*  modify it under the terms of the GNU General Public License as   *
*  published by the Free Software Foundation, either version 3 of   *
*  the License, or (at your option) any later version.              *
*                                                                   *
*  This program is distributed in the hope that it will be useful,  *
*  but WITHOUT ANY WARRANTY; without even the implied warranty of   *
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    *
*  GNU General Public License for more details.                     *
*                                                                   *
*  You should have received a copy of the GNU General Public        *
*  License along with this program.                                 *
*                                                                   *
*  This file is subject to the terms and conditions defined in the  *
*  file 'license', which is part of this source code package.       *
* ================================================================= */

#include "CopMEMMatcher.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>
#include <cassert>
#include <cmath>

#define _prefetch(x,y) __builtin_prefetch(x,1,(4-y))

#include "Hashes.h"

//////////////////// GLOBAL CONSTS //////////////////////////
const uint64_t NOT_MATCHED_POSITION = UINT64_MAX;

//////////////////// GLOBAL VARS ////////////////////////////

#define INIT_HASH_FUNC(n) hashFuncMatrix[n][1] = maRushPrime1HashSimplified<n>; hashFuncMatrix[n][2] = xxhash32<n>; hashFuncMatrix[n][3] = maRushPrime1HashSparsified<n>;hashFuncMatrix[n][4] = metroHash64<n>; hashFuncMatrix[n][5] = cityHash64<n>;
void CopMEMMatcher::initHashFuncMatrix() {
    INIT_HASH_FUNC(20);
    INIT_HASH_FUNC(24);
    INIT_HASH_FUNC(28);
    INIT_HASH_FUNC(32);
    INIT_HASH_FUNC(36);
    INIT_HASH_FUNC(40);
    INIT_HASH_FUNC(44);
    INIT_HASH_FUNC(56);
}


//////////////////// GLOBALS ////////////////////////////

std::ostream *v1logger;

void CopMEMMatcher::initParams(uint32 minMatchLength) {
    if (L > 110) K = 56;
    else if (L > 62) K = 44;
    else if (L > 53) K = 40;
    else if (L > 46) K = 36;
    else if (L > 42) K = 32;
    else if (L > 32) K = 28;
    else K = (L / 4 - 1) * 4;
    if (minMatchLength < 24) {
        cout << "Error: Minimal matching length too short!" << endl;
        exit(EXIT_FAILURE);
    }
    int KmmL = (minMatchLength / 4 - 1) * 4;
    if (KmmL < K) K = KmmL;
    calcCoprimes();
    hashFunc32 = hashFuncMatrix[K][H];
    v1logger = PgSAHelpers::logout;

    LK2 = (L - K) / 2;
    LK2_MINUS_4 = LK2 - 4;
    K_PLUS_LK24 = K + LK2_MINUS_4;

    uint8_t i = HASH_SIZE_MIN_ORDER;
    do {
        hash_size = ((uint32_t) 1) << (i++);
    } while (i <= HASH_SIZE_MAX_ORDER && hash_size < N / k1);
    hash_size_minus_one = hash_size - 1;
}

void CopMEMMatcher::displayParams() {
	std::cout << "copMEM PARAMETERS: ";
	std::cout << "l = " << L << "; ";
	std::cout << "K = " << K << "; ";
	std::cout << "HASH_SIZE = " << hash_size << "; ";
	std::cout << "k1 = " << k1 << "; ";
	std::cout << "k2 = " << k2 << std::endl;
	*PgSAHelpers::logout << "Hash function: maRushPrime1HashSparsified" << std::endl;
    *PgSAHelpers::logout << "Hash collisions per position limit: " << HASH_COLLISIONS_PER_POSITION_LIMIT << std::endl;
    *PgSAHelpers::logout << "Average hash collisions per position limit (in approx mode): " << AVERAGE_HASH_COLLISIONS_PER_POSITION_LIMIT << std::endl;
    *PgSAHelpers::logout << "Unlimited number of hash collisions per position (in approx mode): " << UNLIMITED_NUMBER_OF_HASH_COLLISIONS_PER_POSITION << std::endl;
}

void CopMEMMatcher::calcCoprimes()
{
	/* setting k1 and k2 */
	int tempVar = L - K + 1;
	if (tempVar <= 0) {
		std::cerr << "\nL and K mismatch.\n";
		exit(EXIT_FAILURE);
	}
	if (tempVar >= 20) {
        k1 = (int) (pow(tempVar, 0.5)) + 1;
        k2 = k1 - 1;
        if (k1 * k2 > tempVar) {
            --k2;
            --k1;
        }
    } else if (tempVar >= 15) {
	    k1 = 5; k2 = 3;
	} else if (tempVar >= 12) {
        k1 = 4; k2 = 3;
    } else if (tempVar >= 10) {
        k1 = 5; k2 = 2;
    } else if (tempVar >= 6) {
        k1 = 3; k2 = 2;
    } else {
        k1 = tempVar; k2 = 1;
    }
}

template<typename MyUINT1, typename MyUINT2>
void CopMEMMatcher::genCumm(size_t N, const char* gen, MyUINT2* cumm, vector<MyUINT1> &skippedList) {
	const size_t MULTI1 = 128;
	const size_t k1MULTI1 = k1 * MULTI1;

	uint32_t hashPositions[MULTI1];
	size_t i;

	for (i = 0; i + K + k1MULTI1 < N ; i += k1MULTI1) {
		const char* tempPointer = gen + i;
		for (size_t temp = 0; temp < MULTI1; ++temp) {
			hashPositions[temp] = hashFunc(tempPointer) + 2;
			tempPointer += k1;
			_prefetch((char*)(cumm + hashPositions[temp]), 1);
		}

		for (size_t temp = 0; temp < MULTI1; ++temp) {
			if (cumm[hashPositions[temp]] <= HASH_COLLISIONS_PER_POSITION_LIMIT)
			    ++cumm[hashPositions[temp]];
			else
			    skippedList.push_back(i + k1 * temp);
		}
	}

	//////////////////// processing the end part of R  //////////////////////
	for (; i < N - K + 1; i += k1) {
		uint32_t h = hashFunc(gen + i) + 2;
		if (cumm[h] <= HASH_COLLISIONS_PER_POSITION_LIMIT)
		    ++cumm[h];
        else
            skippedList.push_back(i);
	}
	//////////////////// processing the end part of R //////////////////////
	std::partial_sum(cumm, cumm + hash_size + 2, cumm);
	skippedList.push_back(N);
}

template<typename MyUINT1, typename MyUINT2>
HashBuffer<MyUINT1, MyUINT2> CopMEMMatcher::processRef() {

	const unsigned int MULTI2 = 128;
	const unsigned int k1MULTI2 = k1 * MULTI2;

	MyUINT2* cumm = new MyUINT2[hash_size + 2]();
	vector<MyUINT1> skippedList;
	genCumm(N, start1, cumm, skippedList);
    const size_t hashCount = cumm[hash_size + 1];
    MyUINT1* sampledPositions = new MyUINT1[hashCount + 2];
    *v1logger << "Hash count = " << hashCount << std::endl;

	uint32_t hashPositions[MULTI2];
	MyUINT1 i1;

	size_t s = 0;
	for (i1 = 0; i1 + K + k1MULTI2 < N ; i1 += k1MULTI2) {
		const char* tempPointer = start1 + i1;
		for (unsigned int temp = 0; temp < MULTI2; ++temp) {
		    hashPositions[temp] = hashFunc(tempPointer) + 1;
			tempPointer += k1;
			_prefetch((char*)(cumm + hashPositions[temp]), 1);
		}

		for (unsigned int temp = 0; temp < MULTI2; ++temp) {
			_prefetch((char*)(sampledPositions + *(cumm + hashPositions[temp])), 1);
		}

		MyUINT1 i2 = i1;
		for (size_t temp = 0; temp < MULTI2; ++temp, i2 += k1) {
            if (skippedList[s] == i2) {
                s++;
                continue;
            }
			sampledPositions[cumm[hashPositions[temp]]] = i2;
			++cumm[hashPositions[temp]];
		}
	}

	//////////////////// processing the end part of R
	for (; i1 < N - K + 1; i1 += k1) {
        if (skippedList[s] == i1) {
            s++;
            continue;
        }
	    uint32_t h = hashFunc(start1 + i1) + 1;
		sampledPositions[cumm[h]] = i1;
		++cumm[h];
	}

	return { sampledPositions, cumm };
}

template <class MyUINT1, class MyUINT2>
void CopMEMMatcher::deleteHashBuffer(HashBuffer<MyUINT1, MyUINT2> & buf) {
    delete[] buf.first;
    delete[] buf.second;
}

template<typename MyUINT1, typename MyUINT2>
void CopMEMMatcher::processExactMatchQueryTight(HashBuffer<MyUINT1, MyUINT2> buffer, vector<TextMatch> &resMatches,
                                                const string &destText, bool destIsSrc, bool revComplMatching,
                                                uint32_t minMatchLength){
	const unsigned int MULTI = 256;
	const unsigned int k2MULTI = k2 * MULTI;

	MyUINT1* sampledPositions = buffer.first;
	MyUINT2* cumm = buffer.second;

	uint32_t hArray[MULTI];
	MyUINT2 posArray[MULTI * 2];

	std::uint32_t l1 = 0, l2 = 0, r1 = 0, r2 = 0;

	size_t charExtensions = 0ULL;

    size_t N2 = destText.length();
    const char* start2 = destText.data();

    const int skip = K / k1 - k2;
    const int skipK2 = skip * k2;
    const bool MULTI_MODE = true;
    *v1logger << "Minimal matching length = " << minMatchLength << "; ";
    *v1logger << "Skip factor = " << skip << "; ";
    *v1logger << "Multi-mode = " << (MULTI_MODE?"true":"false") << std::endl;

    size_t i1 = 0;
    const char* end1 = start1 + N;
    const char* end2 = start2 + N2;
    if (MULTI_MODE) {
        for (i1 = 0; i1 + K + k2MULTI < N2 + 1; i1 += k2MULTI) {
            const char *curr2 = start2 + i1;
            size_t tempCount = 0;
            for (size_t i2 = 0; i2 < MULTI; ++i2) {
                hArray[tempCount++] = hashFunc(curr2);
                curr2 += k2;
            }
            for (size_t i2 = 0; i2 < tempCount; ++i2) {
                memcpy(posArray + i2 * 2, cumm + hArray[i2], sizeof(MyUINT2) * 2);
            }

            curr2 = start2 + i1;
            for (size_t i2 = 0; i2 < tempCount; ++i2) {
                if (posArray[i2 * 2] == posArray[i2 * 2 + 1]) {
                    curr2 += k2;
                    continue;
                }

                if (curr2 - LK2 >= start2) memcpy(&l2, curr2 - LK2, sizeof(std::uint32_t));
                if (curr2 + K_PLUS_LK24 + sizeof(std::uint32_t) <= end2) memcpy(&r2, curr2 + K_PLUS_LK24, sizeof(std::uint32_t));

                for (MyUINT1 j = posArray[i2 * 2]; j < posArray[i2 * 2 + 1]; ++j) {
                    ++charExtensions;
                    const char *curr1 = start1 + sampledPositions[j];

                    uint64_t tmpMatchSrcPos = sampledPositions[j];
                    uint64_t tmpMatchDestPos = curr2 - start2;
                    if (destIsSrc && (revComplMatching ? destText.length() - tmpMatchSrcPos < tmpMatchDestPos
                                                       : curr2 - start2 >= tmpMatchSrcPos))
                        continue;
                    if (resMatches.size() > 0 &&
                        tmpMatchDestPos - tmpMatchSrcPos == resMatches.back().posDestText - resMatches.back().posSrcText
                        && tmpMatchDestPos + K < resMatches.back().posDestText + resMatches.back().length) {
                        curr2 += skipK2;
                        i2 += skip;
                        break;
                    }

                    memcpy(&l1, curr1 - LK2, sizeof(std::uint32_t));
                    memcpy(&r1, curr1 + K_PLUS_LK24, sizeof(std::uint32_t));

                    if (r1 == r2 || l1 == l2) {
                        const char *p1 = curr1 + K - 1;
                        const char *p2 = curr2 + K - 1;
                        while (++p1 != end1 && ++p2 != end2 && *p1 == *p2);
                        const char *right = p1;
                        p1 = curr1;
                        p2 = curr2;
                        while (p1 != start1 && p2 != start2 && *--p1 == *--p2);

                        if (right - p1 > minMatchLength && memcmp(curr1, curr2, K) == 0) {
                            resMatches.push_back(TextMatch(p1 + 1 - start1, right - p1 - 1, (p2 + 1 - start2)));

                            curr2 += skipK2;
                            i2 += skip;
                            break;
                        }
                    }
                }
                curr2 += k2;
            }
        }
    }
    //////////////////// processing the end part of Q  //////////////////////
    const char* curr2 = start2 + i1;
    for (; i1 + K < N2 + 1; i1 += k2) {
        memcpy(posArray, cumm + hashFunc(curr2), sizeof(MyUINT2) * 2);

        if (posArray[0] == posArray[1]) {
            curr2 += k2;
            continue;
        }

        if (curr2 - LK2 >= start2) memcpy(&l2, curr2 - LK2, sizeof(std::uint32_t));
        if (curr2 + K_PLUS_LK24 + sizeof(std::uint32_t) <= end2) memcpy(&r2, curr2 + K_PLUS_LK24, sizeof(std::uint32_t));

        for (MyUINT1 j = posArray[0]; j < posArray[1]; ++j) {
            ++charExtensions;
            const char* curr1 = start1 + sampledPositions[j];

            uint64_t tmpMatchSrcPos = sampledPositions[j];
            uint64_t tmpMatchDestPos = curr2 - start2;
            if (destIsSrc && (revComplMatching ? destText.length() - tmpMatchSrcPos < tmpMatchDestPos
                                               : curr2 - start2 >= tmpMatchSrcPos))
                continue;
            if (resMatches.size() > 0 &&
                tmpMatchDestPos - tmpMatchSrcPos == resMatches.back().posDestText - resMatches.back().posSrcText
                && tmpMatchDestPos + K < resMatches.back().posDestText + resMatches.back().length) {
                curr2 += skipK2;
                i1 += skipK2;
                break;
            }

            memcpy(&l1, curr1 - LK2, sizeof(std::uint32_t));
            memcpy(&r1, curr1 + K_PLUS_LK24, sizeof(std::uint32_t));

            if (r1 == r2 || l1 == l2) {
                const char* p1 = curr1 + K - 1;
                const char* p2 = curr2 + K - 1;
                while (++p1 != end1 && ++p2 != end2 && *p1 == *p2);
                const char *right = p1;
                p1 = curr1;
                p2 = curr2;
                while (p1 != start1 && p2 != start2 && *--p1 == *--p2);

                if (right - p1 > minMatchLength && memcmp(curr1, curr2, K) == 0) {
                    resMatches.push_back(TextMatch(p1 + 1 - start1, right - p1 - 1, (p2 + 1 - start2)));
                    curr2 += skipK2;
                    i1 += skipK2;
                    break;
                }
            }
        }
        curr2 += k2;
    }
    //////////////////// processing the end part of Q  //////////////////////

	*v1logger << "Character extensions = " << charExtensions <<  "\n";
}

template<typename MyUINT1, typename MyUINT2>
uint64_t CopMEMMatcher::processApproxMatchQueryTight(HashBuffer<MyUINT1, MyUINT2> buffer, const char *start2,
                                                     const uint_read_len_max N2, uint8_t maxMismatches,
                                                     uint8_t minMismatches, uint8_t &mismatchesCount,
                                                     uint64_t& betterMatchCount, uint64_t& falseMatchCount) {
    if (mismatchesCount < maxMismatches)
        maxMismatches = mismatchesCount - 1;
    MyUINT1* sampledPositions = buffer.first;
    MyUINT2* cumm = buffer.second;

    MyUINT2 posArray[2];

    const uint_read_len_max N2trim8 = (N2 / sizeof(uint64_t)) * sizeof(uint64_t);
    uint_read_len_max k2positionsCount = (N2 + 1 - K) / k2;
    uint64_t falseMatchCountLimit = falseMatchCount + k2positionsCount * AVERAGE_HASH_COLLISIONS_PER_POSITION_LIMIT;

    uint64_t matchPosition = NOT_MATCHED_POSITION;
    size_t i1 = 0;
    const char* curr2 = start2 + i1;
    for (; i1 + K < N2 + 1; i1 += k2) {
        memcpy(posArray, cumm + hashFunc(curr2), sizeof(MyUINT2) * 2);

        if (posArray[0] == posArray[1]) {
            curr2 += k2;
            continue;
        }
        if (falseMatchCountLimit < falseMatchCount) {
            MyUINT2 tmp = posArray[0] + UNLIMITED_NUMBER_OF_HASH_COLLISIONS_PER_POSITION;
            if (posArray[1] > tmp)
                posArray[1] = tmp;
        }
        for (MyUINT1 j = posArray[0]; j < posArray[1]; ++j) {
            const uint_read_len_max positionShift = curr2 - start2;
            if (positionShift > sampledPositions[j])
                continue;
            if (sampledPositions[j] - positionShift + N2 > N)
                continue;
            const char* curr1 = start1 + sampledPositions[j];

            uint64_t temp1, temp2;
            uint8_t res = 0;
            const char* patternPtr = start2;
            const char* patGuardPtr = start2 + N2trim8;
            const char* textPtr = curr1 - positionShift;
            while(res <= maxMismatches && patGuardPtr != patternPtr) {
                memcpy(&temp1, patternPtr, sizeof(uint64_t));
                memcpy(&temp2, textPtr, sizeof(uint64_t));
                uint64_t temp3 = ((temp1 ^ temp2) + 0x7F7F7F7F7F7F7F7F) & 0x8080808080808080;
                res += __builtin_popcountll(temp3);
                patternPtr += sizeof(uint64_t);
                textPtr += sizeof(uint64_t);
            }
            if (res > maxMismatches) {
                falseMatchCount++;
                continue;
            }
            while (patternPtr != start2 + N2) {
                if (*patternPtr++ != *textPtr++) {
                    if (res++ >= maxMismatches) {
                        falseMatchCount++;
                        break;
                    }
                }
            }
            if (res > maxMismatches){
                falseMatchCount++;
                continue;
            }
            if (mismatchesCount != UINT8_MAX)
                betterMatchCount++;
            mismatchesCount = res;
            matchPosition = curr1 - start1 - positionShift;
            if (res <= minMismatches)
                return matchPosition;
            maxMismatches = res - 1;
        }
        curr2 += k2;
    }

    return matchPosition;
}


using namespace std;

CopMEMMatcher::CopMEMMatcher(const char *srcText, const size_t srcLength, const uint32_t targetMatchLength, uint32_t minMatchLength)
    : start1(srcText), N(srcLength), L(targetMatchLength) {
    initHashFuncMatrix();
    if (minMatchLength > targetMatchLength)
        minMatchLength = targetMatchLength;
    initParams(minMatchLength);
    displayParams();

    if ((N) / k1 >= (1ULL << 32)) {
        bigRef = 2;  // huge Reference
        *v1logger  << "WARNING - LARGE reference file (SIZE / k1 > 4GB), 64-bit arrays\n";
        buffer2 = processRef<std::uint64_t, std::uint64_t>();
    } else if (N >= (1ULL << 32)) {
        bigRef = 1;  // large Reference
        *v1logger << "WARNING - BIG reference file (>4GB), 64-bit arrays\n";
        buffer1 = processRef<std::uint64_t, std::uint32_t>();
    } else {
        bigRef = 0;  // small Reference
        buffer0 = processRef<uint32_t, uint32_t>();
    }
}

CopMEMMatcher::~CopMEMMatcher() {
    if (bigRef == 2)
        deleteHashBuffer(buffer2);
    else if (bigRef == 1)
        deleteHashBuffer(buffer1);
    else if (bigRef == 0)
        deleteHashBuffer(buffer0);
}


void
CopMEMMatcher::matchTexts(vector <TextMatch> &resMatches, const string &destText, bool destIsSrc, bool revComplMatching,
                          uint32_t minMatchLength) {
    if (minMatchLength < K) {
        fprintf(stderr, "Minimal matching length cannot be smaller than K (%d < %d)\n\n", minMatchLength, K);
        exit(EXIT_FAILURE);
    }
    resMatches.clear();
    if (bigRef == 2) {
        processExactMatchQueryTight<std::uint64_t, std::uint64_t>(buffer2, resMatches, destText, destIsSrc,
                                                                  revComplMatching, minMatchLength);
    } else if (bigRef == 1) {
        processExactMatchQueryTight<std::uint64_t, std::uint32_t>(buffer1, resMatches, destText, destIsSrc,
                                                                  revComplMatching, minMatchLength);
    }
    else {
        processExactMatchQueryTight<std::uint32_t, std::uint32_t>(buffer0, resMatches, destText, destIsSrc,
                                                                  revComplMatching, minMatchLength);
    }
}

uint64_t CopMEMMatcher::approxMatchPattern(const char *pattern, const uint_read_len_max length, uint8_t maxMismatches, uint8_t minMismatches,
        uint8_t &mismatchesCount, uint64_t& multiMatchCount, uint64_t& falseMatchCount) {
    if (bigRef == 2) {
        return processApproxMatchQueryTight<std::uint64_t, std::uint64_t>(buffer2, pattern, length, maxMismatches,
                                                                          minMismatches, mismatchesCount,
                                                                          multiMatchCount, falseMatchCount);
    } else if (bigRef == 1) {
        return processApproxMatchQueryTight<std::uint64_t, std::uint32_t>(buffer1, pattern, length, maxMismatches,
                                                                          minMismatches, mismatchesCount,
                                                                          multiMatchCount, falseMatchCount);
    }
    else {
        return processApproxMatchQueryTight<std::uint32_t, std::uint32_t>(buffer0, pattern, length, maxMismatches,
                                                                          minMismatches, mismatchesCount,
                                                                          multiMatchCount, falseMatchCount);
    }
}



