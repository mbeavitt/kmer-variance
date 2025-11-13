#!/usr/bin/env python3
"""
Example usage of the kmer_variance package
"""

import numpy as np
from kmer_variance import calculate_diversity

# Example 1: Random test data
print("Example 1: Random test sequences")
print("-" * 50)
sequences = np.random.randint(65, 85, (200, 178), dtype=np.uint8)  # Random A-T
diversity = calculate_diversity(sequences, window_size=50)
print(f"Input: {len(sequences)} sequences")
print(f"Window size: 50")
print(f"Output: {len(diversity)} diversity values")
print(f"Mean diversity: {diversity.mean():.6f}")
print(f"Min diversity: {diversity.min():.6f}")
print(f"Max diversity: {diversity.max():.6f}")
print()

# Example 2: Reading from a file
print("Example 2: Real data from file")
print("-" * 50)
with open('../test_sampled.fa', 'rb') as f:
    data = f.read()

num_sequences = len(data) // 178
sequences = np.frombuffer(data[:num_sequences * 178], dtype=np.uint8).reshape(num_sequences, 178)

print(f"Loaded {num_sequences} sequences from file")

import time
start = time.time()
diversity = calculate_diversity(sequences, window_size=100)
elapsed = time.time() - start

print(f"Window size: 100")
print(f"Computation time: {elapsed:.4f} seconds")
print(f"Output: {len(diversity)} diversity values")
print(f"Mean diversity: {diversity.mean():.6f}")
print(f"First 10 values: {diversity[:10]}")
