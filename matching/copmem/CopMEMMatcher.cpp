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
#include <list>

#define _prefetch(x,y) __builtin_prefetch(x,1,(4-y))

#include "Hashes.h"

/////////////////// OWN TYPES ///////////////////////////

typedef std::tuple<std::string, char*, size_t> SequenceItem2;
typedef std::vector<SequenceItem2> SequenceVector2;

typedef std::tuple <std::size_t, std::size_t, std::vector<size_t>> BlockItem;
typedef std::vector<BlockItem> BlockVector;

//////////////////// GLOBAL CONSTS //////////////////////////
const std::uint32_t HASH_SIZE = 1U << 29;

//////////////////// GLOBAL VARS ////////////////////////////

char* blockBuffer; //buffer used in buffered reading of some Query sequences
size_t blockBufferSize;
BlockVector blockVector;

void CopMEMMatcher::initGlobals() {
	H = 1;
    K = 36;
	L = 47;
	k1 = 4;
	k2 = 3;
	hashFunc = hashFuncMatrix[K][H];
	isVerbose = v1;
}


void CopMEMMatcher::initHashFuncMatrix() {
	hashFuncMatrix[36][1] = maRushPrime1HashSimplified<36>;
	hashFuncMatrix[36][2] = xxhash32<36>;
	hashFuncMatrix[36][3] = xxhash64<36>;
	hashFuncMatrix[36][4] = metroHash64<36>;
	hashFuncMatrix[36][5] = cityHash64<36>;
	hashFuncMatrix[44][1] = maRushPrime1HashSimplified<44>;
	hashFuncMatrix[44][2] = xxhash32<44>;
	hashFuncMatrix[44][3] = xxhash64<44>;
	hashFuncMatrix[44][4] = metroHash64<44>;
	hashFuncMatrix[44][5] = cityHash64<44>;
	hashFuncMatrix[56][1] = maRushPrime1HashSimplified<56>;
	hashFuncMatrix[56][2] = xxhash32<56>;
	hashFuncMatrix[56][3] = xxhash64<56>;
	hashFuncMatrix[56][4] = metroHash64<56>;
	hashFuncMatrix[56][5] = cityHash64<56>;
}


//////////////////// GLOBALS ////////////////////////////

class NullBuffer : public std::streambuf
{
public:
	int overflow(int c) { return c; }
};


NullBuffer null_buffer;
std::ostream null_stream(&null_buffer);

std::ostream *v1logger;

void CopMEMMatcher::displayParams() {
	std::cout << "PARAMETERS" << std::endl;
	std::cout << "l = " << L << "; ";
	std::cout << "K = " << K << "; ";
	std::cout << "HASH_SIZE = " << HASH_SIZE << "; ";
	std::cout << "k1 = " << k1 << "; ";
	std::cout << "k2 = " << k2 << std::endl;
	std::cout << "Hash function: ";
	if (hashFunc == hashFuncMatrix[36][1] || hashFunc == hashFuncMatrix[44][1] || hashFunc == hashFuncMatrix[56][1]) std::cout << "maRushPrime1HashSimplified\n";
	if (hashFunc == hashFuncMatrix[36][2] || hashFunc == hashFuncMatrix[44][2] || hashFunc == hashFuncMatrix[56][2]) std::cout << "xxhash32\n";
	if (hashFunc == hashFuncMatrix[36][3] || hashFunc == hashFuncMatrix[44][3] || hashFunc == hashFuncMatrix[56][3]) std::cout << "xxhash64\n";
	if (hashFunc == hashFuncMatrix[36][4] || hashFunc == hashFuncMatrix[44][4] || hashFunc == hashFuncMatrix[56][4]) std::cout << "metroHash64\n";
	if (hashFunc == hashFuncMatrix[36][5] || hashFunc == hashFuncMatrix[44][5] || hashFunc == hashFuncMatrix[56][5]) std::cout << "cityHash64\n";
}


void CopMEMMatcher::initParams() {
    assert(K % 4 == 0);
    v1logger = &std::cout;
    //v1logger = &null_stream;
    //hashFunc = hashFuncMatrix[K][H];
}

void CopMEMMatcher::calcCoprimes()
{
	/* setting k1 and k2 */
	int tempVar = L - K + 1;
	if (tempVar <= 0) {
		std::cerr << "\nL and K mismatch.\n";
		exit(1);
	}
	k1 = (int)(pow(tempVar, 0.5)) + 1;
	k2 = k1 - 1;
	if (k1 * k2 > tempVar) {
        --k1;
        --k2;
    }
}

template <class MyUINT>
void CopMEMMatcher::genCumm(size_t N, const char* gen, MyUINT* cumm) {
	const size_t MULTI1 = 128;
	const size_t k1MULTI1 = k1 * MULTI1;

	memset(cumm, sizeof(MyUINT)*(HASH_SIZE + 2), (MyUINT)0);
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
			++cumm[hashPositions[temp]];
		}
	}

	//////////////////// processing the end part of R  //////////////////////
	for (; i < N - K + 1; i += k1) {
		uint32_t h = hashFunc(gen + i) + 2;
		++cumm[h];
	}
	//////////////////// processing the end part of R //////////////////////
	std::partial_sum(cumm, cumm + HASH_SIZE + 2, cumm); 
}

template<typename MyUINT1, typename MyUINT2>
HashBuffer<MyUINT1, MyUINT2> CopMEMMatcher::processRef() {

    size_t N = srcText.length();
    const char* start = srcText.data();

	const size_t hashCount = (N - K + 1) / k1;
	const unsigned int MULTI2 = 128;
	const unsigned int k1MULTI2 = k1 * MULTI2;

	*v1logger << "Hash count = " << hashCount << std::endl;

	MyUINT1* sampledPositions = new MyUINT1[hashCount + 2];
	MyUINT2* cumm = new MyUINT2[HASH_SIZE + 2];
	genCumm(N, start, cumm);

	uint32_t hashPositions[MULTI2];
	MyUINT1 i1;

	for (i1 = 0; i1 + K + k1MULTI2 < N ; i1 += k1MULTI2) {
		const char* tempPointer = start + i1;
		for (unsigned int temp = 0; temp < MULTI2; ++temp) {
			hashPositions[temp] = hashFunc(tempPointer) + 1;
			tempPointer += k1;
			_prefetch((char*)(cumm + hashPositions[temp]), 1);
		}

		for (unsigned int temp = 0; temp < MULTI2; ++temp) {
			_prefetch((char*)(sampledPositions + *(cumm + hashPositions[temp])), 1);
		}

		MyUINT1 i2 = i1;
		for (unsigned int temp = 0; temp < MULTI2; ++temp) {
			sampledPositions[cumm[hashPositions[temp]]] = i2;
			++cumm[hashPositions[temp]];
			i2 += k1;
		}
	}

	//////////////////// processing the end part of R
	for (; i1 < N - K + 1; i1 += k1) {
		uint32_t h = hashFunc(start + i1) + 1;
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
void CopMEMMatcher::processQueryTight(HashBuffer<MyUINT1, MyUINT2> buffer, vector <TextMatch> &resMatches,
        const string &destText, bool destIsSrc, bool revComplMatching, uint32_t minMatchLength){
    const char* start1 = srcText.data();

    const int L_PLUS_ONE = L + 1;
	const int LK2 = (L - K) / 2;
	const int LK2_MINUS_4 = LK2 - 4;
	const int K_PLUS_LK24 = K + LK2_MINUS_4;

	const unsigned int MULTI = 256;
	const unsigned int k2MULTI = k2 * MULTI;

	MyUINT1* sampledPositions = buffer.first;
	MyUINT2* cumm = buffer.second;

	uint32_t hArray[MULTI];
	MyUINT2 posArray[MULTI * 2];

	std::uint32_t l1, l2, r1, r2;

	size_t charExtensions = 0ULL;

    size_t N2 = destText.length();
    const char* start2 = destText.data();

    const int skip = K / k1 - k2;
    const int skipK2 = skip * k2;
    const bool MULTI_MODE = true;

    size_t i1 = 0;
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

                memcpy(&l2, curr2 - LK2, sizeof(std::uint32_t));
                memcpy(&r2, curr2 + K_PLUS_LK24, sizeof(std::uint32_t));

                for (MyUINT1 j = posArray[i2 * 2]; j < posArray[i2 * 2 + 1]; ++j) {
                    ++charExtensions;
                    const char *curr1 = start1 + sampledPositions[j];

                    uint64_t tmpMatchSrcPos = sampledPositions[j];
                    uint64_t tmpMatchDestPos = curr2 - start2;
                    if (destIsSrc && (revComplMatching ? destText.length() - tmpMatchSrcPos < tmpMatchDestPos
                                                       : curr2 - start2 >= tmpMatchSrcPos))
                        continue;

                    memcpy(&l1, curr1 - LK2, sizeof(std::uint32_t));
                    memcpy(&r1, curr1 + K_PLUS_LK24, sizeof(std::uint32_t));

                    if (r1 == r2 || l1 == l2) {
                        const char *p1 = curr1 + K - 1;
                        const char *p2 = curr2 + K - 1;
                        while (*++p1 == *++p2);
                        const char *right = p1;
                        p1 = curr1;
                        p2 = curr2;
                        while (*--p1 == *--p2);

                        if (right - p1 > L && memcmp(curr1, curr2, K) == 0) {
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

        memcpy(&l2, curr2 - LK2, sizeof(std::uint32_t));
        memcpy(&r2, curr2 + K_PLUS_LK24, sizeof(std::uint32_t));

        for (MyUINT1 j = posArray[0]; j < posArray[1]; ++j) {
            ++charExtensions;
            const char* curr1 = start1 + sampledPositions[j];

            uint64_t tmpMatchSrcPos = sampledPositions[j];
            uint64_t tmpMatchDestPos = curr2 - start2;
            if (destIsSrc && (revComplMatching ? destText.length() - tmpMatchSrcPos < tmpMatchDestPos
                                               : curr2 - start2 >= tmpMatchSrcPos))
                continue;

            memcpy(&l1, curr1 - LK2, sizeof(std::uint32_t));
            memcpy(&r1, curr1 + K_PLUS_LK24, sizeof(std::uint32_t));

            if (r1 == r2 || l1 == l2) {
                const char* p1 = curr1 + K - 1;
                const char* p2 = curr2 + K - 1;
                while (*++p1 == *++p2);
                const char* right = p1;
                p1 = curr1;
                p2 = curr2;
                while (*--p1 == *--p2) ;

                if (right - p1 > L && memcmp(curr1, curr2, K) == 0) {
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

using namespace std;

CopMEMMatcher::CopMEMMatcher(const string &srcText, const uint32_t targetMatchLength) : srcText(srcText),
    targetMatchLength(targetMatchLength) {
    initHashFuncMatrix();
    initGlobals();
    initParams();
    calcCoprimes();
    displayParams();

    if ((srcText.length()) / k1 >= (1ULL << 32)) {
        bigRef = 2;  // huge Reference
        *v1logger  << "WARNING - LARGE reference file (SIZE / k1 > 4GB), 64-bit arrays\n";
        buffer2 = processRef<std::uint64_t, std::uint64_t>();
    } else if (srcText.length() >= (1ULL << 32)) {
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
    resMatches.clear();
    if (bigRef == 2) {
        processQueryTight<std::uint64_t, std::uint64_t>(buffer2, resMatches, destText, destIsSrc, revComplMatching, minMatchLength);
    } else if (bigRef == 1) {
        processQueryTight<std::uint64_t, std::uint32_t>(buffer1, resMatches, destText, destIsSrc, revComplMatching, minMatchLength);
    }
    else {
        processQueryTight<std::uint32_t, std::uint32_t>(buffer0, resMatches, destText, destIsSrc, revComplMatching, minMatchLength);
    }
}



