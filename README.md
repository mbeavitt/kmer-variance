# kmer-variance

Fast sliding window kmer diversity analysis using incremental Hamming distances.

## Installation

```bash
make
pip install .
```

## Usage

```python
import numpy as np
from kmer_variance import calculate_diversity

# Generate random DNA sequences
bases = np.array([65, 67, 71, 84], dtype=np.uint8)  # A, C, G, T
sequences = bases[np.random.randint(0, 4, (1000, 178))]

# Calculate diversity across sliding windows
diversity = calculate_diversity(sequences, window_size=100)

print(f"Mean diversity: {diversity.mean():.4f}")
```

## API

### `calculate_diversity(sequences, window_size=100)`

Calculate kmer diversity across a sliding window of sequences.

**Parameters:**
- `sequences` (numpy.ndarray): Array of shape `(num\_sequences, 178)` with dtype uint8
- `window_size` (int, optional): Size of the sliding window. Default: 100

**Returns:**
- `numpy.ndarray`: Array of diversity values, one per window position

**Raises:**
- `ValueError`: If sequences has wrong shape or `window_size` is invalid
- `TypeError`: If sequences is not a numpy array
- `RuntimeError`: If the C library fails

## Algorithm

The function computes the average pairwise distance of repeats in an array across a sliding window.

For each window position:
1. First window: Calculate all O(n²) pairwise distances
2. Subsequent windows: Remove distances for leaving sequence, add distances for entering sequence
3. This reduces complexity from O(W × n²) to O(W × n) where W is number of windows and n is window size

Each sequence is converted to a 256-bit vector representing the presence/absence of all 256 possible 4-mers (4⁴ = 256) to take advantage of popcount assembly instructions `__builtin_popcount(x ^ y)`.
On platforms with AVX2 support, the 256 bit width of the bitstrings can allow SIMD operations, specifically the expression:

`v4du x = a.vec ^ b.vec; // bitwise XOR`

which processes each of the four 64-bit chunks in the bitstring simultaneously. The following expressions:

```
return __builtin_popcountll(x[0]) +
       __builtin_popcountll(x[1]) +
       __builtin_popcountll(x[2]) +
       __builtin_popcountll(x[3]);
```

Do not benefit from SIMD in theory, but in LLVM/Clang >= 17.0.0, Muła's
algorithm is automatically applied to roughly half the number of CPU cycles
required for this operation, which improves the performance of the program by
about 1.67x. (see https://arxiv.org/abs/1611.07612)

## Requirements

- Python >= 3.7
- NumPy >= 1.19.0
- C compiler (clang or gcc) with support for builtins (clang >= 17.0.0 HIGHLY recommended for better optimisation of AVX2 related expressions)

## License

MIT
