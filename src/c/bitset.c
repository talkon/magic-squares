#include "bitset.h"

bitset_t bitset_create(size_t size) {
  bitset_t bitset;
  bitset.size = (size + sizeof(uint64_t) * 8 - 1) / (sizeof(uint64_t) * 8);
  bitset.array = calloc(bitset.size, sizeof(uint64_t));
  if (bitset.array == NULL) {
    return bitset;
  }
  return bitset;
}

void bitset_free(bitset_t bitset) {
  free(bitset.array);
}

size_t bitset_maximum(const bitset_t bitset) {
  for (size_t k = bitset.size; k > 0; k--) {
    uint64_t w = bitset.array[k - 1];
    if (w != 0) {
      return 63 - __builtin_clzll(w) + (k - 1) * 64;
    }
  }
  return 0;
}

void bitset_inplace_xor(bitset_t b1, const bitset_t b2) {
  for (size_t k = 0; k < b1.size; ++k) {
    b1.array[k] ^= b2.array[k];
  }
}
