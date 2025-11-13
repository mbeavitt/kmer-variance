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

void bit256_clear(bit256_t *x, int i) {
    int element = i / 64;
    int bit = i % 64;
    x->v[element] &= ~(1ULL << bit);
}

int bit256_get(const bit256_t *x, int i) {
    int element = i / 64;
    int bit = i % 64;
    return (x->v[element] >> bit) & 1;
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
    for (int i = 0; i < 256; i++) {
    }
}

void skim_4mers(const char *s, bit256_t *a) {
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

// Calculate average pairwise Hamming distance for a window
// Uses incremental sliding: compute first window fully, then slide by adding/removing edges
double sliding_window_diversity_allpairs(const bit256_t *repeat_array, int start_idx, int window_size) {
    int total_pairs = 0;
    int total_distance = 0;

    // Calculate all pairwise distances for the window
    for (int j = 0; j < window_size; j++) {
        for (int k = j + 1; k < window_size; k++) {
            int dist = hamming256(repeat_array[start_idx + j], repeat_array[start_idx + k]);
            total_distance += dist;
            total_pairs++;
        }
    }

    // Return average distance normalized by total number of possible k-mers (256)
    if (total_pairs == 0) return 0.0;
    return (double)total_distance / total_pairs / 256.0;
}

// Calculate average consecutive Hamming distance for a window
// Compares each sequence to the next one in the window
double sliding_window_diversity_consecutive(const bit256_t *repeat_array, int start_idx, int window_size) {
    if (window_size < 2) return 0.0;

    int total_distance = 0;
    int num_comparisons = window_size - 1;

    // Calculate consecutive distances
    for (int j = 0; j < num_comparisons; j++) {
        int dist = hamming256(repeat_array[start_idx + j], repeat_array[start_idx + j + 1]);
        total_distance += dist;
    }

    // Return average distance normalized by total number of possible k-mers (256)
    return (double)total_distance / num_comparisons / 256.0;
}

int main(int argc, char **argv) {
    if (argc != 3) {
        fprintf(stderr, "Usage: %s <num_sequences> <file.fa>\n", argv[0]);
        return 1;
    }

    int num_sequences = atoi(argv[1]);
    const char *filename = argv[2];
    char temp[178];
    int idx;
    bit256_t repeat_array[num_sequences];

    memset(repeat_array, 0, sizeof(repeat_array));

    FILE *fp = fopen(filename, "r");
    if (!fp) {
        perror("Error opening file");
        return 1;
    }

    idx = 0;
    while (fread(temp, 1, 178, fp) == 178) {
        if (idx >= num_sequences) {
            fprintf(stderr, "Warning: more sequences in file than expected\n");
            break;
        }
        find_4mers(temp, &repeat_array[idx]);
        idx++;
    }

    fclose(fp);

    // Sliding window analysis with incremental updates - run 1000 times for benchmarking
    int window_size = 100;
    int num_windows = idx - window_size + 1;
    long long checksum = 0;

    if (num_windows > 0) {
        for (int iteration = 0; iteration < 1000; iteration++) {
            // Variables to track the running sum
            int total_distance = 0;
            int total_pairs = (window_size * (window_size - 1)) / 2;

            // Compute the first window completely
            for (int j = 0; j < window_size; j++) {
                for (int k = j + 1; k < window_size; k++) {
                    total_distance += hamming256(repeat_array[j], repeat_array[k]);
                }
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

                checksum += total_distance;
            }
        }
        printf("Checksum: %lld\n", checksum);
    } else {
        fprintf(stderr, "Not enough sequences (%d) for window size %d\n", idx, window_size);
    }

    return 0;
}
