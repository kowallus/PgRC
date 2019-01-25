/* ================================================================= *
*  CopMEM.cpp : Main file                                           *
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

#include <algorithm>
#include <fstream>
#include <iostream>
#include <numeric>
#include <string>
#include <tuple>
#include <vector>
#include <cassert>
#include <cmath>


#if (defined(linux) || defined(__linux) || defined(__linux__))
#define _prefetch(x,y) __builtin_prefetch(x,1,(4-y))
#else
#include <xmmintrin.h>
#define _prefetch(x,y) _mm_prefetch(x,y)
#endif

#include "StopWatch.h"
#include "Hashes.h"


/////////////////// OWN TYPES ///////////////////////////
typedef std::pair<std::string, size_t> SequenceItem;
typedef std::vector<SequenceItem> SequenceVector;

typedef std::tuple<std::string, char*, size_t> SequenceItem2;
typedef std::vector<SequenceItem2> SequenceVector2;

typedef std::tuple <std::size_t, std::size_t, std::vector<size_t>> BlockItem;
typedef std::vector<BlockItem> BlockVector;

typedef std::tuple<size_t, char*, SequenceVector> GenomeData; //<size of buffer, memory pointer, starting pointer, sequence list>
template<class MyUINT1, class MyUINT2>
 using HashBuffer =  std::pair<MyUINT1*, MyUINT2* >;

enum verbosity { v0, v1, v2 };
enum reverse { no, yes, both };


/////////////////// FUNCTIONS ///////////////////////////
void initHashFuncMatrix();
void displayHelp(const char* progName);
void displayParams();
void initGlobals();
void processCmd(int argc, char* argv[]);
GenomeData readMultiFasta(std::string fn, const char paddingChar, bool removeNs, const char* seqType);
SequenceVector2 readBlock(std::ifstream& f, BlockItem& bi, const char paddingChar, bool removeNs);

//////////////////// GLOBAL CONSTS //////////////////////////
const std::uint32_t HASH_SIZE = 1U << 29;
const size_t MATCHES_BUFFER_SIZE = 1ULL << 24;


//////////////////// GLOBAL VARS ////////////////////////////
int K;
int H;
int L;
int k1;
int k2;

std::uint32_t(*hashFunc)(const char*);
std::uint32_t(*hashFuncMatrix[64][6])(const char*);

std::string matchesFN;
std::string matchesBuffer;

std::ofstream f_matches;
std::string R_FN;
std::string Q_FN;

verbosity isVerbose;
reverse isRC;
bool isFast;

char* blockBuffer; //buffer used in buffered reading of some Query sequences
size_t blockBufferSize;
BlockVector blockVector;


void initGlobals() {
	K = 44;
	H = 1;
	L = 100;
	k1 = 8;
	k2 = 7;
	hashFunc = hashFuncMatrix[K][H];
	isVerbose = v1;
	isRC = no;
	isFast = false;
	matchesBuffer.reserve(MATCHES_BUFFER_SIZE);
}


void initHashFuncMatrix() {
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

std::ostream *v1logger, *v2logger;


void displayHelp(const char* progName) {
	std::cout << "copMEM 0.2, by Szymon Grabowski and Wojciech Bieniecki, July 2018." << std::endl;
	std::cout << "Usage: " << progName << " [-l n] [-H l] [-K 36|44|56] [-e] [-v|-q] [-f]|[-b]|[-r] <-o MEMs_file> <Ref_genome> <Query_genome>\n";
	std::cout << "Attention: -o is a required parameter. l is optional (default: 100).\n";
	std::cout << "-o MEMs_file - the output file with matches.\n";
	std::cout << "-v - verbose mode. Display more details.\n";
	std::cout << "-q - quiet mode. No screen output.\n";
	std::cout << "-f - fast mode. Requires additional memory. May speed up I/O.\n";
	std::cout << "-b - compute forward and reverse-complement matches. Not available with -f.\n";
	std::cout << "-r - compute only reverse-complement matches. Not available with -f.\n";
	std::cout << "-l n - minimal length of matches. Default value is 100.\n";
	std::cout << "-H n - Hash Function. 1: maRushPrime1HashSimplified (default), 2: xxhash32, 3: xxhash64, 4: metroHash64, 5: cityHash64.\n";
	std::cout << "-K 36|44|56 Default K is 44.\n";
	std::cout << "-e - simulates E-MEM (forces k_2 = 1). If you don't understand it, don't use!\n";
}


void displayParams() {
	std::cout << "PARAMETERS" << std::endl;
	std::cout << "Reference filename: " << R_FN << std::endl;
	std::cout << "Query filename: " << Q_FN << std::endl;
	std::cout << "l = " << L << std::endl;
	std::cout << "K = " << K << std::endl;
	std::cout << "HASH_SIZE = " << HASH_SIZE << std::endl;
	std::cout << "k1 = " << k1 << std::endl;
	std::cout << "k2 = " << k2 << std::endl;
	std::cout << "Fast mode = " << (isFast? "Yes" : "No") << std::endl;
	std::cout << "Reverse Complement = " << ((isRC==both)? "Both" : (isRC == yes) ? "Yes" : "No") << std::endl;
	std::cout << "Hash function: ";
	if (hashFunc == hashFuncMatrix[36][1] || hashFunc == hashFuncMatrix[44][1] || hashFunc == hashFuncMatrix[56][1]) std::cout << "maRushPrime1HashSimplified\n";
	if (hashFunc == hashFuncMatrix[36][2] || hashFunc == hashFuncMatrix[44][2] || hashFunc == hashFuncMatrix[56][2]) std::cout << "xxhash32\n";
	if (hashFunc == hashFuncMatrix[36][3] || hashFunc == hashFuncMatrix[44][3] || hashFunc == hashFuncMatrix[56][3]) std::cout << "xxhash64\n";
	if (hashFunc == hashFuncMatrix[36][4] || hashFunc == hashFuncMatrix[44][4] || hashFunc == hashFuncMatrix[56][4]) std::cout << "metroHash64\n";
	if (hashFunc == hashFuncMatrix[36][5] || hashFunc == hashFuncMatrix[44][5] || hashFunc == hashFuncMatrix[56][5]) std::cout << "cityHash64\n";
}


void processCmd(int argc, char* argv[]) {
	assert(K % 4 == 0);
	v1logger = &std::cout;
	v2logger = &null_stream;
	bool isOset = false;
	bool isEmemLike = false;
	const char* incompleteCmd = " option requires one integer argument.\n";
	for (int i = 1; i < argc - 2; ++i) {
		std::string arg = argv[i];

		if (arg == "-o") {
			if (i + 1 < argc) {
				matchesFN = argv[++i];
				isOset = true;
			}
			else {
				std::cerr << "-o requires file name.\n";
				exit(1);
			}
		}

		if (arg == "-v") {
			if (isVerbose == v0) {
				std::cerr << "-v and -q parameters are mutually exclusive.\n";
				exit(1);
			}
			isVerbose = v2;
			v2logger = &std::cout;
		}

		if (arg == "-e") {
			isEmemLike = true;
		}
		if (arg == "-q") {
			if (isVerbose == v2) {
				std::cerr << "-v and -q parameters are mutually exclusive.\n";
				exit(1);
			}
			isVerbose = v0;
			v1logger = &null_stream;
		}
		if (arg == "-f") {
			isFast = true;
		}
		if (arg == "-b") {
			if (isRC == yes) {
				std::cerr << "-b and -r parameters are mutually exclusive.\n";
				exit(1);
			}
			isRC = both;
		}
		if (arg == "-r") {
			if (isRC == both) {
				std::cerr << "-b and -r parameters are mutually exclusive.\n";
				exit(1);
			}
			isRC = yes;
		}

		if (arg == "-h") {
			displayHelp("");
			exit(0);
		}
		if (arg == "-l") {
			if (i + 1 < argc) {
				L = atoi(argv[++i]);
				if (L < 50) {
					std::cerr << "Incorrect L value (must be >= 50).\n";
					exit(1);
				}
			}
			else {
				std::cerr << "L" << incompleteCmd;
				exit(1);
			}
		}
		if (arg == "-K") {
			if (i + 1 < argc) {
				K = atoi(argv[++i]);
				if ((K != 44) && (K != 36) && (K != 56)) {
					std::cerr << "Incorrect K value (must be 36 or 44, or 56).\n";
					exit(1);
				}
			}
			else {
				std::cerr << "K" << incompleteCmd;
				exit(1);
			}
		}
		if (arg == "-H") {
			if (i + 1 < argc) {
				H = atoi(argv[++i]);
				if (H < 1 || H > 5) {
					std::cerr << "Incorrect H value (must be from 1 to 5).\n";
					exit(1);
				}
			}
			else {
				std::cerr << "H" << incompleteCmd;
				exit(1);
			}
		}
	}

	if (isOset == false) {
		std::cerr << "-o not given or specified correctly.\n";
		exit(1);
	}
	if (isFast && (isRC != no)) {
		std::cerr << "-r or -b parameters don't go with -f.\n";
		exit(1);
	}
	//setting hash func
	hashFunc = hashFuncMatrix[K][H];

	//touching files:
	R_FN = argv[argc - 2];
	Q_FN = argv[argc - 1];
	std::ifstream f;

	f.open(R_FN);
	if (f.fail()) {
		std::cerr << "\nReference file '" << R_FN << "' does not exist. Quit.";
		exit(1);
	}
	f.close();
	f.open(Q_FN);
	if (f.fail()) {
		std::cerr << "\nQuery file '" << Q_FN << "' does not exist. Quit.";
		exit(1);
	}
	f.close();
	f_matches.open(matchesFN);
	if (f_matches.fail()) {
		std::cerr << "\nFAILED. The -o parameter specifies a file that cannot be created.\n";
		exit(1);
	}
	f_matches.close();

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

	if (isEmemLike == true) {
		k1 = L - K + 1;
		k2 = 1;
	}
}


bool arrayComp(size_t* p1, size_t* p2) {
	if (p1[1] < p2[1]) return true;
	if (p1[1] > p2[1]) return false;
	return p1[0] < p2[0];
}


template <class MyUINT>
void genCumm(GenomeData& genome, MyUINT* cumm) {
	const size_t MULTI1 = 128;
	const size_t k1MULTI1 = k1 * MULTI1;

	size_t N = std::get<0>(genome);
	char* gen = std::get<1>(genome);

	//std::fill(cumm, cumm + HASH_SIZE + 2, 0);
	memset(cumm, sizeof(MyUINT)*(HASH_SIZE + 2), (MyUINT)0);
	uint32_t hashPositions[MULTI1];
	size_t i;

	for (i = 0; i + K + k1MULTI1 < N ; i += k1MULTI1) {
		char* tempPointer = gen + i;
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


void dumpMEM(SequenceItem& item1, SequenceItem& item2, size_t* match) {
	size_t baseindex1 = match[0] - item1.second;
	size_t baseindex2 = match[1] -item2.second;
	matchesBuffer.append(" ");
	matchesBuffer.append(item1.first);
	matchesBuffer.append("\t");
	matchesBuffer.append(std::to_string(baseindex1));
	matchesBuffer.append("\t");
	matchesBuffer.append(std::to_string(baseindex2));
	matchesBuffer.append("\t");
	matchesBuffer.append(std::to_string(match[2]));
	matchesBuffer.append("\n");
}


void dumpMEMTight(SequenceItem& item1, size_t* match, size_t counter) {
	size_t baseindex1 = match[0] - item1.second;
	size_t baseindex2 = match[1] - counter;
	matchesBuffer.append(" ");
	matchesBuffer.append(item1.first);
	matchesBuffer.append("\t");
	matchesBuffer.append(std::to_string(baseindex1));
	matchesBuffer.append("\t");
	matchesBuffer.append(std::to_string(baseindex2));
	matchesBuffer.append("\t");
	matchesBuffer.append(std::to_string(match[2]));
	matchesBuffer.append("\n");
}


bool lowerBoundComp(SequenceItem lhs, SequenceItem rhs) {
	return lhs.second < rhs.second;
}


SequenceItem findSeqDesc(size_t index, SequenceVector& seq) {
	SequenceItem dummySequenceItem = { "", index };
	SequenceItem item = seq[0];
	auto lower = std::lower_bound(seq.begin(), seq.end(), dummySequenceItem, lowerBoundComp);
	size_t diff = lower - seq.begin();
	return seq[diff - 1];
}


void displayMatchInfo(std::string& name, size_t count) {
	switch (count) {
		case 0:  *v2logger << name << ": no matches.\n"; break;
		case 1:  *v2logger << name << ": match.\n"; break;
		default: *v2logger << name << ": matches.\n"; break;
	}
}


void postProcess(std::vector<size_t*> &matches, GenomeData& t1, GenomeData& t2, size_t index2) {
	SequenceVector& seq2 = std::get<2>(t2);
	SequenceItem seq2item = seq2[index2];
	matchesBuffer.append("> ");
	matchesBuffer.append(seq2item.first);
	matchesBuffer.append("\n");

	if (matches.size() == 0) {
		displayMatchInfo(seq2item.first, 0);
		return;
	}
	char* gen1 = std::get<1>(t1);
	SequenceVector& seq1 = std::get<2>(t1);
	if (matches.size() > 1)
		std::sort(matches.begin(), matches.end(), arrayComp);

	SequenceItem seq1item;

	size_t* prev_match = matches[0];
	seq1item = findSeqDesc(prev_match[0], seq1);

	auto foundPos = std::find(gen1 + prev_match[0], gen1 + prev_match[0] + prev_match[2], 'N');
	if (foundPos == gen1 + prev_match[0] + prev_match[2])
		dumpMEM(seq1item, seq2item, prev_match);

	std::uint64_t count = 1ULL;

	for (auto match: matches) {
		if (memcmp(prev_match, match, 3 * sizeof(size_t)) != 0) {
			auto foundPos = std::find(gen1 + match[0], gen1 + match[0] + match[2], 'N');
			if (foundPos != gen1 + match[0] + match[2])
				continue;
			if (prev_match[0] != match[0])
				seq1item = findSeqDesc(match[0], seq1);
			dumpMEM(seq1item, seq2item, match);
			++count;
		}
		prev_match = match;
	}
	displayMatchInfo(seq2item.first, count);
	for (auto match: matches)
		delete[] match;
	matches.clear();
}


void postProcessTight(std::vector<size_t*> &matches, GenomeData& t1, std::string seqName, size_t counter) {
	matchesBuffer.append("> ");
	matchesBuffer.append(seqName);
	matchesBuffer.append("\n");

	if (matches.size() == 0) {
		displayMatchInfo(seqName, 0);
		return;
	}
	char* gen1 = std::get<1>(t1);
	SequenceVector& seq1 = std::get<2>(t1);
	if (matches.size() > 1)
		std::sort(matches.begin(), matches.end(), arrayComp);

	SequenceItem seq1item;

	size_t* prev_match = matches[0];
	seq1item = findSeqDesc(prev_match[0], seq1);

	auto foundPos = std::find(gen1 + prev_match[0], gen1 + prev_match[0] + prev_match[2], 'N');
	if (foundPos == gen1 + prev_match[0] + prev_match[2])
		dumpMEMTight(seq1item,  prev_match, counter);

	std::uint64_t count = 1ULL;

	for (auto match: matches) {
		if (memcmp(prev_match, match, 3 * sizeof(size_t)) != 0) {
			auto foundPos = std::find(gen1 + match[0], gen1 + match[0] + match[2], 'N');
			if (foundPos != gen1 + match[0] + match[2])
				continue;
			if (prev_match[0] != match[0])
				seq1item = findSeqDesc(match[0], seq1);
			dumpMEMTight(seq1item,  match, counter);
			++count;
		}
		prev_match = match;
	}
	displayMatchInfo(seqName, count);
	for (auto match: matches)
		delete[] match;
	matches.clear();
}


template<class MyUINT1, class MyUINT2>
HashBuffer<MyUINT1, MyUINT2> processRef(GenomeData& rGenome) {
	CStopWatch stopwatch;
	stopwatch.start();
	size_t N = std::get<0>(rGenome);
	char* start = std::get<1>(rGenome);

	const size_t hashCount = (N - K + 1) / k1;
	const unsigned int MULTI2 = 128;
	const unsigned int k1MULTI2 = k1 * MULTI2;

	*v1logger << "Hash count = " << hashCount << std::endl;

	MyUINT1* sampledPositions = new MyUINT1[hashCount + 2];
	MyUINT2* cumm = new MyUINT2[HASH_SIZE + 2];
	*v2logger << "\tprocessRef: init = " << stopwatch.stop() << std::endl;
	stopwatch.resume();
	genCumm(rGenome, cumm);

	*v2logger << "\tprocessRef: cumm(1) = " << stopwatch.stop() << std::endl;
	stopwatch.resume();

	uint32_t hashPositions[MULTI2];
	MyUINT1 i1;

	for (i1 = 0; i1 + K + k1MULTI2 < N ; i1 += k1MULTI2) {
		char* tempPointer = start + i1;
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
		char* tempPointer = start + i1;
		uint32_t h = hashFunc(start + i1) + 1;
		sampledPositions[cumm[h]] = i1;
		++cumm[h];
	}
	*v2logger << "\tprocessRef: cumm(2) = " << stopwatch.stop() << std::endl;

	return { sampledPositions, cumm };
}

template<class MyUINT1, class MyUINT2>
void processQuery(HashBuffer<MyUINT1, MyUINT2> buffer, GenomeData& rGenome, GenomeData& qGenome) {
	const int L_PLUS_ONE = L + 1;
	const int L2 = L / 2;
	const int LK2 = (L - K) / 2;
	const int LK2_MINUS_2 = LK2 - 2;
	const int LK2_MINUS_4 = LK2 - 4;
	const int K_PLUS_LK22 = K + LK2_MINUS_2;
	const int K_PLUS_LK24 = K + LK2_MINUS_4;

	size_t N1 = std::get<0>(rGenome);
	size_t N2 = std::get<0>(qGenome);
	char* start1 = std::get<1>(rGenome);
	char* start2 = std::get<1>(qGenome);
	SequenceVector seq1 = std::get<2>(rGenome);
	SequenceVector seq2 = std::get<2>(qGenome);

	const size_t hashCount = (N1 - K + 1) / k1;
	const std::uint32_t MULTI = 256;
	const std::uint32_t k2MULTI = k2 * MULTI;

	char* curr1 = start1;
	char* curr2 = start2;

	MyUINT1* sampledPositions = buffer.first;
	MyUINT2* cumm = buffer.second;

	std::vector<size_t*> matches;

	std::uint32_t hArray[MULTI];
	MyUINT2 posArray[MULTI * 2];

	size_t charExtensions = 0ULL;

	size_t currSeqIndex = 0;
	size_t nextSeqBeginPos = (seq2.size() > 1 ? seq2[currSeqIndex + 1].second : SIZE_MAX);

	f_matches.open(matchesFN);

	std::uint32_t l1, l2, r1, r2;

	size_t i1;
	for (i1 = 0; i1 + K + k2MULTI < N2 + 1 ; i1 += k2MULTI) {
		size_t i1temp = i1;
		char* prev = start2 + i1;
		curr2 = start2 + i1;
		size_t tempCount = 0;
		for (size_t i2 = 0; i2 < MULTI; ++i2) {
			if (*curr2 != 'N') {
				hArray[tempCount++] = hashFunc(curr2);
			}
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
				curr1 = start1 + sampledPositions[j];

				memcpy(&l1, curr1 - LK2, sizeof(std::uint32_t));
				memcpy(&r1, curr1 + K_PLUS_LK24, sizeof(std::uint32_t));

				if (r1 == r2 || l1 == l2) {
					char* p1 = curr1 + K;
					char* p2 = curr2 + K;

					while (*p1 == *p2) {
						++p1;
						++p2;
					}
					char* right = p1;

					p1 = curr1 - 1;
					p2 = curr2 - 1;
					while (*p1 == *p2) {
						--p1;
						--p2;
					}

					if (right - p1 >= L_PLUS_ONE && memcmp(curr1, curr2, K) == 0) {
						size_t* tempMatch = new size_t[3];
						tempMatch[0] = p1 + 1 - start1;
						tempMatch[1] = p2 + 1 - start2;
						tempMatch[2] = right - p1 - 1;
						while ((size_t)(p2 + 1 - start2) >= nextSeqBeginPos) {
							if (currSeqIndex + 1 < seq2.size()) {
								postProcess(matches, rGenome, qGenome, currSeqIndex);
								if (matchesBuffer.size() > (size_t)(MATCHES_BUFFER_SIZE * 0.95)) { f_matches << matchesBuffer; matchesBuffer.clear(); }
								++currSeqIndex;
								nextSeqBeginPos = seq2[currSeqIndex + 1].second;
							}
							else
								break;
						}
						matches.push_back(tempMatch);
					}
				}
			}
			curr2 += k2;
		}
	}

	//////////////////// processing the end part of Q  //////////////////////
	for (; i1 + K < N2 + 1; i1 += k2) {
		char* prev = start2 + i1;
		curr2 = start2 + i1;
		size_t tempCount = 0;

		if (*curr2 != 'N') {
			memcpy(posArray + tempCount * 2, cumm + hashFunc(curr2), sizeof(MyUINT2) * 2);
			++tempCount;
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
				curr1 = start1 + sampledPositions[j];

				memcpy(&l1, curr1 - LK2, sizeof(std::uint32_t));
				memcpy(&r1, curr1 + K_PLUS_LK24, sizeof(std::uint32_t));

				if (r1 == r2 || l1 == l2) {
					char* p1 = curr1 + K;
					char* p2 = curr2 + K;

					while (*p1 == *p2) {
						++p1;
						++p2;
					}
					char* right = p1;

					p1 = curr1 - 1;
					p2 = curr2 - 1;
					while (*p1 == *p2) {
						--p1;
						--p2;
					}

					if (right - p1 >= L_PLUS_ONE && memcmp(curr1, curr2, K) == 0) {
						size_t* tempMatch = new size_t[3];
						tempMatch[0] = p1 + 1 - start1;
						tempMatch[1] = p2 + 1 - start2;
						tempMatch[2] = right - p1 - 1;

						while ((size_t)(p2 + 1 - start2) > nextSeqBeginPos) {
							if (currSeqIndex + 1 < seq2.size()) {
								postProcess(matches, rGenome, qGenome, currSeqIndex);
								if (matchesBuffer.size() > (size_t)(MATCHES_BUFFER_SIZE * 0.95)) { f_matches << matchesBuffer; matchesBuffer.clear(); }
								++currSeqIndex;
								nextSeqBeginPos = seq2[currSeqIndex + 1].second;
							}
							else
								break;
						}
						matches.push_back(tempMatch);
					}
				}
			}
			curr2 += k2;
		}
	}
	//////////////////// processing the end part of Q  //////////////////////
	while (currSeqIndex + 2 <= seq2.size()) {
		postProcess(matches, rGenome, qGenome, currSeqIndex);
		++currSeqIndex;
		nextSeqBeginPos = seq2[currSeqIndex + 1].second;
	}

	postProcess(matches, rGenome, qGenome, currSeqIndex);
	f_matches << matchesBuffer;
	matchesBuffer.clear();
	f_matches.close();
	*v1logger << "Character extensions = " << charExtensions << "\n";
}


void reverseComplement(char* start, const char* lut, const std::size_t N) {
	char* left = start + 1; // sequence starts from paddingChar
	char* right = start + N - 1;
	while (right > left) {
		char tmp = lut[*left];
		*left = lut[*right];
		*right = tmp;
		++left;
		--right;
	}
	if (left == right)
		*left = lut[*left];
}


template<class MyUINT1, class MyUINT2>
void processQueryTight(HashBuffer<MyUINT1, MyUINT2> buffer, GenomeData& rGenome){
	const int L_PLUS_ONE = L + 1;
	const int L2 = L / 2;
	const int LK2 = (L - K) / 2;
	const int LK2_MINUS_2 = LK2 - 2;
	const int LK2_MINUS_4 = LK2 - 4;
	const int K_PLUS_LK22 = K + LK2_MINUS_2;
	const int K_PLUS_LK24 = K + LK2_MINUS_4;
	const int paddingSize = L_PLUS_ONE;
	const char paddingChar = 125;
	const bool removeNs = true;

	size_t N1 = std::get<0>(rGenome);

	char* start1 = std::get<1>(rGenome);
	SequenceVector seq1 = std::get<2>(rGenome);
	const size_t hashCount = (N1 - K + 1) / k1;
	const unsigned int MULTI = 256;
	const unsigned int k2MULTI = k2 * MULTI;

	char* curr1 = start1;
	char complement[256];
	memset(complement, paddingChar, 256);
	complement['A'] = 'T';
	complement['a'] = 'T';
	complement['C'] = 'G';
	complement['c'] = 'G';
	complement['G'] = 'C';
	complement['g'] = 'C';
	complement['T'] = 'A';
	complement['t'] = 'A';
	complement['N'] = 'N';
	complement['n'] = 'N';


	MyUINT1* sampledPositions = buffer.first;
	MyUINT2* cumm = buffer.second;

	std::vector<size_t*> matches;

	uint32_t hArray[MULTI];
	MyUINT2 posArray[MULTI * 2];

	f_matches.open(matchesFN);
	std::ifstream qf(Q_FN, std::ios::binary);

	std::uint32_t l1, l2, r1, r2;

	size_t i1;

	size_t counter = 0ULL;
	size_t charExtensions = 0ULL;

	for (auto bv: blockVector) {
		SequenceVector2 sv2 = readBlock(qf, bv, paddingChar, true);
		for (auto si3 : sv2) {
			size_t N2 = std::get<2>(si3);
			char* start2 = std::get<1>(si3);
			if (isRC != yes) {
				char* curr2 = start2;
				i1 = 0;
				for (i1 = 0; i1 + K + k2MULTI < N2 + 1; i1 += k2MULTI) {
					size_t i1temp = i1;
					char* prev = start2 + i1;
					curr2 = start2 + i1;
					size_t tempCount = 0;
					for (size_t i2 = 0; i2 < MULTI; ++i2) {
						hArray[tempCount++] = hashFunc(curr2);
						curr2 += k2;
					}
					for (size_t i2 = 0; i2 < tempCount; ++i2) {
						memcpy(posArray + i2 * 2, cumm + hArray[i2], sizeof(MyUINT2) * 2);
					}

					curr2 = start2 + i1;	// !!!
					for (size_t i2 = 0; i2 < tempCount; ++i2) {
						if (posArray[i2 * 2] == posArray[i2 * 2 + 1]) {
							curr2 += k2;
							continue;
						}
						
						memcpy(&l2, curr2 - LK2, sizeof(std::uint32_t));
						memcpy(&r2, curr2 + K_PLUS_LK24, sizeof(std::uint32_t));
						
						for (MyUINT1 j = posArray[i2 * 2]; j < posArray[i2 * 2 + 1]; ++j) {
							++charExtensions;
							curr1 = start1 + sampledPositions[j];
							
							memcpy(&l1, curr1 - LK2, sizeof(std::uint32_t));
							memcpy(&r1, curr1 + K_PLUS_LK24, sizeof(std::uint32_t));
							
							if (r1 == r2 || l1 == l2) {
								char* p1 = curr1 + K;
								char* p2 = curr2 + K;
								
								while (*p1 == *p2) {
									++p1;
									++p2;
								}
								char* right = p1;
								
								p1 = curr1 - 1;
								p2 = curr2 - 1;
								while (*p1 == *p2) {
									--p1;
									--p2;
								}
								
								if (right - p1 >= L_PLUS_ONE && memcmp(curr1, curr2, K) == 0) {
									size_t* tempMatch = new size_t[3];
									tempMatch[0] = p1 + 1 - start1;
									tempMatch[1] = (p2 + 1 - start2) + counter;
									tempMatch[2] = right - p1 - 1;
									matches.push_back(tempMatch);
								}
							}
						}
						curr2 += k2;
					}
				}
				//////////////////// processing the end part of Q  //////////////////////
				for (; i1 + K < N2 + 1; i1 += k2) {
					char* prev = start2 + i1;
					curr2 = start2 + i1;
					size_t tempCount = 0;
					memcpy(posArray + tempCount * 2, cumm + hashFunc(curr2), sizeof(MyUINT2) * 2);
					++tempCount;
					
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
							curr1 = start1 + sampledPositions[j];
							
							memcpy(&l1, curr1 - LK2, sizeof(std::uint32_t));
							memcpy(&r1, curr1 + K_PLUS_LK24, sizeof(std::uint32_t));
							
							if (r1 == r2 || l1 == l2) {
								char* p1 = curr1 + K;
								char* p2 = curr2 + K;
								
								while (*p1 == *p2) {
									++p1;
									++p2;
								}
								char* right = p1;
								
								p1 = curr1 - 1;
								p2 = curr2 - 1;
								while (*p1 == *p2) {
									--p1;
									--p2;
								}
								
								if (right - p1 >= L_PLUS_ONE && memcmp(curr1, curr2, K) == 0) {
									size_t* tempMatch = new size_t[3];
									tempMatch[0] = p1 + 1 - start1;
									tempMatch[1] = (p2 + 1 - start2) + counter;
									tempMatch[2] = right - p1 - 1;
									matches.push_back(tempMatch);
								}
							}
						}
						curr2 += k2;
					}
				}
				//////////////////// processing the end part of Q  //////////////////////
				postProcessTight(matches, rGenome, std::get<0>(si3), counter);
				if (matchesBuffer.size() > (size_t)(MATCHES_BUFFER_SIZE * 0.95)) { f_matches << matchesBuffer; matchesBuffer.clear(); }
			}
			if (isRC != no) {
				reverseComplement(start2, complement, N2);
				char* curr2 = start2;
				i1 = 0;
				for (i1 = 0; i1 + K + k2MULTI < N2 + 1; i1 += k2MULTI) {
					size_t i1temp = i1;
					char* prev = start2 + i1;
					curr2 = start2 + i1;
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
							curr1 = start1 + sampledPositions[j];

							memcpy(&l1, curr1 - LK2, sizeof(std::uint32_t));
							memcpy(&r1, curr1 + K_PLUS_LK24, sizeof(std::uint32_t));

							if (r1 == r2 || l1 == l2) {
								char* p1 = curr1 + K;
								char* p2 = curr2 + K;

								while (*p1 == *p2) {
									++p1;
									++p2;
								}
								char* right = p1;

								p1 = curr1 - 1;
								p2 = curr2 - 1;
								while (*p1 == *p2) {
									--p1;
									--p2;
								}

								if (right - p1 >= L_PLUS_ONE && memcmp(curr1, curr2, K) == 0) {
									size_t* tempMatch = new size_t[3];
									tempMatch[0] = p1 + 1 - start1;
									tempMatch[1] = (p2 + 1 - start2) + counter;
									tempMatch[2] = right - p1 - 1;
									matches.push_back(tempMatch);
								}
							}
						}
						curr2 += k2;
					}
				}
				//////////////////// processing the end part of Q  //////////////////////
				for (; i1 + K < N2 + 1; i1 += k2) {
					char* prev = start2 + i1;
					curr2 = start2 + i1;
					size_t tempCount = 0;
					memcpy(posArray + tempCount * 2, cumm + hashFunc(curr2), sizeof(MyUINT2) * 2);
					++tempCount;

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
							curr1 = start1 + sampledPositions[j];

							memcpy(&l1, curr1 - LK2, sizeof(std::uint32_t));
							memcpy(&r1, curr1 + K_PLUS_LK24, sizeof(std::uint32_t));

							if (r1 == r2 || l1 == l2) {
								char* p1 = curr1 + K;
								char* p2 = curr2 + K;

								while (*p1 == *p2) {
									++p1;
									++p2;
								}
								char* right = p1;

								p1 = curr1 - 1;
								p2 = curr2 - 1;
								while (*p1 == *p2) {
									--p1;
									--p2;
								}

								if (right - p1 >= L_PLUS_ONE && memcmp(curr1, curr2, K) == 0) {
									size_t* tempMatch = new size_t[3];
									tempMatch[0] = p1 + 1 - start1;
									tempMatch[1] = (p2 + 1 - start2) + counter;
									tempMatch[2] = right - p1 - 1;
									matches.push_back(tempMatch);
								}
							}
						}
						curr2 += k2;
					}
				}
				//////////////////// processing the end part of Q  //////////////////////
				postProcessTight(matches, rGenome, std::get<0>(si3) + " Reverse", counter);
				if (matchesBuffer.size() > (size_t)(MATCHES_BUFFER_SIZE * 0.95)) { f_matches << matchesBuffer; matchesBuffer.clear(); }
			}
			counter += N2;
		}
	}
	f_matches << matchesBuffer;
	matchesBuffer.clear();
	f_matches.close();
	*v1logger << "Character extensions = " << charExtensions <<  "\n";
}


void replaceBadSymbol(char* gen, char* dst, char symbol, char paddingChar) {
	char* movingPtr = gen;
	while (1) {
		char* tempPtr = std::find(movingPtr, dst, symbol);
		while (*tempPtr == symbol) {
			*tempPtr = paddingChar;
			++tempPtr;
		}
		if (tempPtr == dst)
			break;
		movingPtr = tempPtr + 1;
	}
}


GenomeData readMultiFasta(std::string fn, const char paddingChar, bool removeNs, const char* seqType) {
	CStopWatch stopWatch;
	stopWatch.start();
	const char beginChar = '>';
	const char terminatorChar = 0;
	const char eolChar1 = 10;
	const char eolChar2 = 13;
	const char spaceChar = ' ';
	const char filterChar = (char)0xDF;
	const int paddingSize = L + 1;

	*v1logger << "Reading " << seqType << " genome ...\n";  // seqType is "Reference" or "Query"

	//create a buffer for the whole file + padding at left and right
	std::ifstream f(fn, std::ios::ate | std::ios::binary);
	if (f.fail()) {
		std::cerr << "\nFile '" << fn << "' does not exist. Quit.";
		exit(1);
	}
	size_t N = f.tellg();
	f.seekg(0, std::ios::beg);
	char* buf1 = new char[N + 2 * paddingSize];
	if (buf1 == nullptr) {
		std::cerr << "\nFile '" << fn << "' is too large. Quit.";
		exit(1);
	}
	memset(buf1, paddingChar, paddingSize);
	f.read(buf1 + paddingSize, N);
	*v2logger << "\treadMultiFasta: Reading file from disk " << stopWatch.stop() << "\n";
	stopWatch.resume();
	
	buf1[paddingSize + N] = terminatorChar; // null-terminate the string
	memset(buf1 + paddingSize + N + 1, paddingChar, paddingSize - 1);
	memset(buf1 + paddingSize + N, terminatorChar, 10);  // >= sizeof(uint_64) is enough
	f.close();

	char* gen = buf1 + paddingSize;
	SequenceVector seq;

	char* dst = gen;
	char* src = gen;

	char tempLine[512];  // FASTA lines shouldn't be that long (even 128 should be ok)

	while (1) {
		if (*src == beginChar) {
			size_t idx = 0;
			while (*src != eolChar1 && *src != eolChar2 && *src != ' ' && *src != terminatorChar) {
				tempLine[idx] = *src++;
				++idx;
			}
			tempLine[idx] = 0;
			seq.push_back({ tempLine + 1, (dst - gen) });  // + 1, as we omit the starting '>'
														   //search for EOL
			while (*src != eolChar1 && *src != eolChar2)
				src++;

			*dst++ = paddingChar;
		}
		else {
			while (*src == eolChar1 || *src == eolChar2) {
				++src;
			}
			if (*src == beginChar)
				continue;
			if (*src == terminatorChar)
				break;

			uint64_t temp2;
			memcpy(&temp2, src, 8);
			while (((~temp2) & 0x4040404040404040) == 0) {
				temp2 = temp2 & 0xDFDFDFDFDFDFDFDF;
				memcpy(dst, &temp2, 8);
				dst += 8;
				src += 8;
				memcpy(&temp2, src, 8);
			}
			while (((~(*src)) & (char)0x40) == 0) {
				*dst++ = (*src++) & (char)0xDF;
			}
		}
	}

	memset(dst, paddingChar, N + paddingSize - (dst - gen));

	if (removeNs == true) {
		replaceBadSymbol(gen, dst, 'N', paddingChar);
		replaceBadSymbol(gen, dst, 'n', paddingChar);
	}
	*v2logger << "\treadMultiFasta: Analysing file " << stopWatch.stop() << "\n";
	*v2logger << "\treadMultiFasta: Genome data size " << (dst - gen) << std::endl;
	*v2logger << "\treadMultiFasta: Sequences " << seq.size() << std::endl;
	
	return { (dst - gen), gen, seq };
}


SequenceVector2 readBlock(std::ifstream& f, BlockItem& bi, const char paddingChar, bool removeNs) {
	const char terminatorChar = 0;
	const char eolChar1 = 10;
	const char eolChar2 = 13;
	const char beginChar = '>';
	const int paddingSize = L + 1;
	size_t N = std::get<1>(bi);
	std::vector<size_t>& sv2 = std::get<2>(bi);

	char* gen = nullptr;
	char* genTemp = blockBuffer + paddingSize;
	//reading buffer to memory
	f.seekg(std::get<0>(bi), std::ios::beg);
	f.read(genTemp, std::get<1>(bi));
	memset(genTemp + std::get<1>(bi), 0, paddingSize);
	
	gen = genTemp;
	char* dst = gen;
	char* src = gen;
	SequenceVector2 sv3;

	char tempLine[512];  // FASTA lines shouldn't be that long (even 128 should be ok)

	//processing sequences: removing header, new lines
	for (size_t i = 0; i < sv2.size(); i++) {
		size_t startHeader = sv2[i];
		size_t nextSeqence = (sv2.size()>(i+1) ? sv2[i + 1]: std::get<1>(bi));
		
		std::string seqName;
		size_t seqStart;
		size_t seqSize;
		//check, if the sequence starts from '>'
		src = gen + startHeader;
		char* lastChar = gen + (nextSeqence);
		if (*src != beginChar) {
			std::cerr << "Invalid fasta File\n";
			exit(0);
		}
		size_t idx = 0;
		while (*src != eolChar1 && *src != eolChar2 && *src != ' ' && src != lastChar) { //reading header line
			tempLine[idx] = *src++;
			++idx;
		}
		tempLine[idx] = 0;
		seqName = tempLine + 1;
		while (*src != eolChar1 && *src != eolChar2) //search for EOL after header
			src++;
		seqStart = dst - gen;
		char *startPtr = dst;
		*dst++ = paddingChar;  //padding char at the begin to make indexing matches from 1

		while (1) { //scanning the sequence
			while (*src == eolChar1 || *src == eolChar2) {
				if (src == lastChar)
					break;
				++src;
			}
			if (src == lastChar)
				break;

			uint64_t temp2; //scan 8 bytes at once
			memcpy(&temp2, src, 8);
			while (((~temp2) & 0x4040404040404040) == 0) {
				temp2 = temp2 & 0xDFDFDFDFDFDFDFDF;
				memcpy(dst, &temp2, 8);
				dst += 8;
				src += 8;
				memcpy(&temp2, src, 8);
			}
			while (((~(*src)) & (char)0x40) == 0) {
				*dst++ = (*src++) & (char)0xDF;
			}
		}
		seqSize = dst - startPtr;
		sv3.push_back({ seqName, startPtr, seqSize });
	}
	*dst = paddingChar;
	if (removeNs == true) {
		replaceBadSymbol(gen, dst, 'N', paddingChar);
		replaceBadSymbol(gen, dst, 'n', paddingChar);
	}

	return sv3;
}


void createBlockBuffer(std::vector<size_t>& seqStarts) {
	//1. Obtain maximum sequence size
	size_t seqSize = seqStarts.size();
	size_t maxSeqLen = 1ULL << 26; //Block length not less than 64 MB

	for (size_t i = 0; i < seqSize - 1; i++) {
		size_t le = seqStarts[i + 1] - seqStarts[i];
		if (maxSeqLen < le)
			maxSeqLen = le;
	}
	// obtain blocks
	// each block consists of:
	// starting ptr in file
	// length of the block
	// list of sequences
	// pointers in sequences must be shifted to beginning of the block


	BlockItem bi;
	for (size_t i = 0; i < seqSize - 1; i++) {
		size_t le = seqStarts[i + 1] - seqStarts[i];

		if (std::get<1>(bi) + le > maxSeqLen) {  //Will not fit - push block to the list and start new block
			blockVector.push_back(bi);
			std::get<0>(bi) = seqStarts[i];
			std::get<1>(bi) = le;
			std::get<2>(bi).clear();
		}
		else
			std::get<1>(bi) += le;
		std::get<2>(bi).push_back(seqStarts[i] - std::get<0>(bi));
	}
	blockVector.push_back(bi);
	seqStarts.clear();

	// reserve memory for block buffer
	blockBufferSize = maxSeqLen + 2 * (L + 1);
	blockBuffer = new char[blockBufferSize];
	memset(blockBuffer, 0, blockBufferSize);
	*v2logger << "\tcreateBlockBuffer: Sequences " << seqSize << std::endl;
	*v2logger << "\tcreateBlockBuffer: Sequence blocks " << blockVector.size() << std::endl;
	*v2logger << "\tcreateBlockBuffer: Max block size " << maxSeqLen << std::endl;
}


void deleteBlockBuffer() {
	blockVector.clear();
	delete[] blockBuffer;
}


void scanMultiFasta(std::string fn){
	*v1logger << "Scanning Query genome ...\n";
	CStopWatch stopWatch;
	stopWatch.start();
	const char beginChar = '>';
	const size_t scanBufSize = 1ULL << 26;
	char* buf1 = new char[scanBufSize];
	if (buf1 == nullptr) {
		std::cerr << "\nFile '" << fn << "' is too large. Quit.";
		exit(1);
	}

	std::ifstream f(fn, std::ios::binary);
	if (f.fail()) {
		std::cerr << "\nFile '" << fn << "' does not exist. Quit.";
		exit(1);
	}
	size_t bytesRead = 0;
	size_t currPos = 0;
	std::vector<size_t> seqStarts;
	while (f) {
		f.read(buf1, scanBufSize);
		if (f)
			bytesRead = scanBufSize;
		else
			bytesRead = f.gcount();
		char* foundPos = buf1;
		while (1) {
			foundPos = std::find(foundPos, buf1 + bytesRead, beginChar);
			if (foundPos != (buf1 + bytesRead)) {
				seqStarts.push_back(currPos + (foundPos - buf1));
				foundPos++;
			}
			else
				break;
		}
		currPos += bytesRead;
	}
	seqStarts.push_back(currPos); //pointer to the end of file
	f.close();
	delete[] buf1;
	createBlockBuffer(seqStarts);
	seqStarts.clear();
	
	*v2logger << "\tscanMultiFasta: Analysing file " << stopWatch.stop() << "\n";
	*v2logger << "\tscanMultiFasta: Sequence blocks " << blockVector.size() << std::endl;	
}


template <class MyUINT1, class MyUINT2>
void deleteHashBuffer(HashBuffer<MyUINT1, MyUINT2> & buf) {
	delete[] buf.first;
	delete[] buf.second;
}


void deleteReading(GenomeData & r) {
	delete[] (std::get<1>(r) - (L + 1));
}


int main(int argc, char* argv[]) {
	if (argc < 3) {
		displayHelp(argv[0]);
		return 1;
	}
	initHashFuncMatrix();
	initGlobals();
	processCmd(argc, argv);

	if (isVerbose > v0) {
		displayParams();
		std::cout.setf(std::ios_base::unitbuf);
	}
	CStopWatch stopwatch;
	stopwatch.start();
	std::size_t maxSeqSize = 0;
	std::size_t N1 = 0;
	std::size_t N2 = 0;
	GenomeData  rGenome, qGenome;
	if (isFast)
		qGenome = readMultiFasta(Q_FN, 125, true, "Query");
	else
		scanMultiFasta(Q_FN);
	rGenome = readMultiFasta(R_FN, 123, true, "Reference");
	
	*v1logger << "Time of I/O = " << stopwatch.stop() << std::endl;
	stopwatch.resume();
	
	
	int bigRef = 0;  // small Reference
	if ((std::get<0>(rGenome)) >= (1ULL << 32)) {  
		bigRef = 1;  // large Reference
		if ((std::get<0>(rGenome)) / k1 >= (1ULL << 32))
			bigRef = 2;  // huge Reference
	}

	if (bigRef == 2) {
		*v1logger  << "WARNING - LARGE reference file (SIZE / k1 > 4GB), 64-bit arrays\n";
		std::pair<std::uint64_t*, std::uint64_t*> buffer = processRef<std::uint64_t, std::uint64_t>(rGenome);
		
		*v1logger  << "Time of processRef = " << stopwatch.stop() << std::endl;
		if(isFast)
			processQuery<std::uint64_t, std::uint64_t>(buffer, rGenome, qGenome);
		else
			processQueryTight<std::uint64_t, std::uint64_t>(buffer, rGenome);
			
		*v1logger << "Time of processQuery = " << stopwatch.stop() << std::endl;
		stopwatch.resume();
		if (isFast == false)
			deleteBlockBuffer();
		deleteHashBuffer(buffer);
	} else
	if (bigRef == 1) {
		*v1logger << "WARNING - BIG reference file (>4GB), 64-bit arrays\n";
		std::pair<std::uint64_t*, std::uint32_t*> buffer = processRef<std::uint64_t, std::uint32_t>(rGenome);
		*v1logger << "Time of processRef = " << stopwatch.stop() << std::endl;
		stopwatch.resume();
		if(isFast)
			processQuery<std::uint64_t, std::uint32_t>(buffer, rGenome, qGenome);
		else
			processQueryTight<std::uint64_t, std::uint32_t>(buffer, rGenome);
			
		*v1logger << "Time of processQuery = " << stopwatch.stop() << std::endl;
		stopwatch.resume();
		if (isFast == false)
			deleteBlockBuffer();
		deleteHashBuffer(buffer);
	}
	else {
		std::pair<std::uint32_t*, std::uint32_t*> buffer = processRef<uint32_t, uint32_t>(rGenome);
		*v1logger << "Time of processRef = " << stopwatch.stop() << std::endl;
		stopwatch.resume();
		if(isFast)
			processQuery<std::uint32_t, std::uint32_t>(buffer, rGenome, qGenome);
		else
			processQueryTight<std::uint32_t, std::uint32_t>(buffer, rGenome);
			
		*v1logger << "Time of processQuery = " << stopwatch.stop() << std::endl;
		stopwatch.resume();
		if (isFast == false)
			deleteBlockBuffer();
		deleteHashBuffer(buffer);
	}
	deleteReading(rGenome);
	if(isFast)
		deleteReading(qGenome);

	*v1logger << "Time of deleting = " << stopwatch.stop() << "\nTotal time = " << stopwatch.totalTime() << "\nFINISHED\n";
	return 0;
}
