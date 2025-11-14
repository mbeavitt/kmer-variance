#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

// Map a character to 2-bit code
#define NUC2BIT(c) \
    ((c) == 'A' ? 0 : \
     (c) == 'C' ? 1 : \
     (c) == 'G' ? 2 : \
     (c) == 'T' ? 3 : 0)

// Map 4-char string to 0..255
#define KMER4_TO_BYTE(s) \
    ((NUC2BIT((s)[0]) << 6) | \
     (NUC2BIT((s)[1]) << 4) | \
     (NUC2BIT((s)[2]) << 2) | \
     (NUC2BIT((s)[3])))

// Check if AVX2 is available
#ifdef __AVX2__
#include <x86intrin.h>

// Define a 256-bit unsigned vector (4 × 64-bit = 32 bytes)
typedef unsigned long long v4du __attribute__((vector_size(32), aligned(32)));

typedef union {
    v4du vec;
    uint64_t v[4];
} __attribute__((aligned(32))) bit256_t;

// Compute Hamming distance between two 256-bit values (AVX2 version)
static inline int hamming256(bit256_t a, bit256_t b) {
    v4du x = a.vec ^ b.vec; // bitwise XOR
    return __builtin_popcountll(x[0]) +
           __builtin_popcountll(x[1]) +
           __builtin_popcountll(x[2]) +
           __builtin_popcountll(x[3]);
}

#else
// Portable version without AVX2

typedef struct {
    uint64_t v[4];
} bit256_t;

// Compute Hamming distance between two 256-bit values (portable version)
static inline int hamming256(bit256_t a, bit256_t b) {
    int dist = 0;
    for (int i = 0; i < 4; ++i)
        dist += __builtin_popcountll(a.v[i] ^ b.v[i]);
    return dist;
}

#endif

void bit256_set(bit256_t *x, int i) {
    int element = i / 64;
    int bit = i % 64;
    x->v[element] |= (1ULL << bit);
}

void find_4mers(const char *s, bit256_t *a) {
    unsigned char kmer_mapping;
    char kmer[5];

    for (int i = 0; i <= 178 - 4; i++) {
        for (int j = 0; j < 4; j++) {
            kmer[j] = s[i + j];
        }
        kmer[4] = '\0';
        kmer_mapping = KMER4_TO_BYTE(kmer);
        bit256_set(a, kmer_mapping);
    }
}

// Function to be called from Python
// sequences: flattened array of sequences (num_sequences * 178 bytes)
// num_sequences: number of sequences
// window_size: size of sliding window
// iterations: number of times to run the analysis
// results: output array for diversity values (will be allocated and returned)
double* run_sliding_window(const char *sequences, int num_sequences, int window_size, int iterations) {
    // Allocate repeat array with proper alignment for AVX2
    size_t alloc_size = num_sequences * sizeof(bit256_t);
    // Ensure size is a multiple of 32 for aligned_alloc
    if (alloc_size % 32 != 0) {
        alloc_size = ((alloc_size + 31) / 32) * 32;
    }
    bit256_t *repeat_array = (bit256_t*)aligned_alloc(32, alloc_size);
    if (!repeat_array) {
        return NULL;
    }
    memset(repeat_array, 0, alloc_size);

    // Convert sequences to bit256_t format
    for (int i = 0; i < num_sequences; i++) {
        find_4mers(sequences + i * 178, &repeat_array[i]);
    }

    int num_windows = num_sequences - window_size + 1;
    if (num_windows <= 0) {
        free(repeat_array);
        return NULL;
    }

    // Allocate results array
    double *results = (double*)malloc(num_windows * sizeof(double));
    if (!results) {
        free(repeat_array);
        return NULL;
    }

    // Run the analysis 'iterations' times, keeping only the last result
    for (int iteration = 0; iteration < iterations; iteration++) {
        // Variables to track the running sum
        int total_distance = 0;
        int total_pairs = (window_size * (window_size - 1)) / 2;

        // Compute the first window completely
        for (int j = 0; j < window_size; j++) {
            for (int k = j + 1; k < window_size; k++) {
                total_distance += hamming256(repeat_array[j], repeat_array[k]);
            }
        }

        // Store first window result
        if (iteration == iterations - 1) {
            results[0] = (double)total_distance / total_pairs / 256.0;
        }

        // Slide the window incrementally
        for (int i = 1; i < num_windows; i++) {
            int leaving_idx = i - 1;
            int entering_idx = i + window_size - 1;

            // Remove distances from leaving sequence to all sequences in the old window
            for (int j = 1; j < window_size; j++) {
                total_distance -= hamming256(repeat_array[leaving_idx], repeat_array[leaving_idx + j]);
            }

            // Add distances from entering sequence to all sequences in the new window
            for (int j = 0; j < window_size - 1; j++) {
                total_distance += hamming256(repeat_array[entering_idx], repeat_array[i + j]);
            }

            // Store result on last iteration
            if (iteration == iterations - 1) {
                results[i] = (double)total_distance / total_pairs / 256.0;
            }
        }
    }

    free(repeat_array);
    return results;
}

// Helper function to free results array from Python
void free_results(double *results) {
    free(results);
}
