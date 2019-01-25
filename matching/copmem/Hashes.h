#pragma once
#ifndef _HASHES_H
#define _HASHES_H
#include <cstdint>
#define XXH_INLINE_ALL
#define XXH_PRIVATE_API
#include "xxhash.h"
#include "metrohash64.h"
#include "city.h"
extern const std::uint32_t HASH_SIZE;
const std::uint32_t HASH_SIZE_MINUS_ONE = HASH_SIZE - 1;
template <size_t K>
inline std::uint32_t xxhash32(const char* key) {
	XXH32_hash_t result = XXH32((const void*)key, K, 9876543210UL);
	return (std::uint32_t)(result & HASH_SIZE_MINUS_ONE);
}

template <size_t K>
inline std::uint32_t xxhash64(const char* key) {
	XXH64_hash_t result = XXH64((const void*)key, K, 9876543210UL);
	return (std::uint32_t)(result & HASH_SIZE_MINUS_ONE);
}

unsigned long long result = 0ULL;
template <size_t K>
inline std::uint32_t metroHash64(const char* key) {
	MetroHash64::Hash((const uint8_t*)key, K, (uint8_t*)(&result), 9876543210ULL);
	return (std::uint32_t)(result & HASH_SIZE_MINUS_ONE);
}

template <size_t K>
inline std::uint32_t cityHash64(const char* key) {
	return (std::uint32_t)(CityHash64(key, K) & HASH_SIZE_MINUS_ONE);
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
	return (std::uint32_t)(hash & HASH_SIZE_MINUS_ONE);
}

#endif //!_HASHES_H
