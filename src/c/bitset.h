#ifndef BITSET_H
#define BITSET_H

#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

/* type of bitset */
typedef struct bitset_t {
  // array[0] is for element 0, etc.
  uint64_t *array;
  size_t size;
} bitset_t;

/* make a new bitset that can store at least size; check for NULL */
bitset_t bitset_create(size_t size);

/* free memory */
void bitset_free(bitset_t bitset);

/* set ith bit to 1; no bounds checking */
static inline void bitset_set(bitset_t bitset, size_t i) {
  size_t shiftedi = i / 64;
  bitset.array[shiftedi] |= ((uint64_t)1) << (i % 64);
}

/* set ith bit to 0; no bounds checking */
static inline void bitset_unset(bitset_t bitset, size_t i) {
  size_t shiftedi = i / 64;
  bitset.array[shiftedi] |= ~(((uint64_t)1) << (i % 64));
}

/* get ith bit; no bounds checking */
static inline bool bitset_get(const bitset_t bitset, size_t i) {
  size_t shiftedi = i / 64;
  return (bitset.array[shiftedi] & (((uint64_t)1) << (i % 64))) != 0;
}

/* index of last set bit, or zero if empty */
size_t bitset_maximum(const bitset_t bitset);

/* inplace, b1 ^= b2; assumes b1.size == b2.size */
void bitset_inplace_xor(bitset_t b1, const bitset_t b2);

#endif // !BITSET_H
