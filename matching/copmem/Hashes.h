#pragma once
#ifndef _HASHES_H
#define _HASHES_H
#include <cstdint>
#define XXH_INLINE_ALL
#define XXH_PRIVATE_API
#include "xxhash.h"
#include "metrohash64.h"
#include "city.h"

template <size_t K>
inline std::uint32_t xxhash32(const char* key) {
	XXH32_hash_t result = XXH32((const void*)key, K, (uint32_t) 9876543210UL);
	return (std::uint32_t)(result);
}

template <size_t K>
inline std::uint32_t xxhash64(const char* key) {
	XXH64_hash_t result = XXH64((const void*)key, K, 9876543210UL);
	return (std::uint32_t)(result);
}

unsigned long long result = 0ULL;
template <size_t K>
inline std::uint32_t metroHash64(const char* key) {
	MetroHash64::Hash((const uint8_t*)key, K, (uint8_t*)(&result), 9876543210ULL);
	return (std::uint32_t)(result);
}

template <size_t K>
inline std::uint32_t cityHash64(const char* key) {
	return (std::uint32_t)(CityHash64(key, K));
}

// based on http://www.amsoftware.narod.ru/algo2.html
template <size_t K>
inline std::uint32_t maRushPrime1HashSimplified(const char *str) {
	std::uint64_t hash = K;
	for (std::uint32_t j = 0; j < K/4; ) {
		std::uint32_t k;
		memcpy(&k, str, 4);
		k += j++;
		hash ^= k;
		hash *= 171717;
		str += 4;
	}
	return (std::uint32_t)(hash);
}

const uint32_t SPARSIFY_MASK_A = 0x00FFFFFF;
const int SPARSIFY_MASK_A_COUNT = 3;
const uint32_t SPARSIFY_MASK_B = 0x0000FFFF;
// based on http://www.amsoftware.narod.ru/algo2.html
template <size_t K>
inline std::uint32_t maRushPrime1HashSparsified(const char *str) {
	std::uint64_t hash = K;
	std::uint32_t j = 0;
	std::uint32_t k;
	while (j < SPARSIFY_MASK_A_COUNT) {
		memcpy(&k, str, 4);
		k &= SPARSIFY_MASK_A;
		k += j++;
		hash ^= k;
		hash *= 171717;
		str += 4;
	}
	while (j < K/4) {
		memcpy(&k, str, 4);
		k &= SPARSIFY_MASK_B;
		k += j++;
		hash ^= k;
		hash *= 171717;
		str += 4;
	}
	return (std::uint32_t)(hash);
}

#endif //!_HASHES_H
