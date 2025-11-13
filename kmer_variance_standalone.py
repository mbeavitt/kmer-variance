#!/usr/bin/env python3
"""
Python wrapper for kmer_variance C library
"""

import ctypes
import numpy as np
from pathlib import Path

# Load the shared library
lib_path = Path(__file__).parent / "libkmer_variance.so"
lib = ctypes.CDLL(str(lib_path))

# Define function signatures
lib.run_sliding_window.argtypes = [
    ctypes.POINTER(ctypes.c_char),  # sequences
    ctypes.c_int,                    # num_sequences
    ctypes.c_int,                    # window_size
    ctypes.c_int                     # iterations
]
lib.run_sliding_window.restype = ctypes.POINTER(ctypes.c_double)

lib.free_results.argtypes = [ctypes.POINTER(ctypes.c_double)]
lib.free_results.restype = None


def read_sequences(filename, num_sequences=None):
    """
    Read sequences from a binary file where each sequence is 178 bytes.

    Args:
        filename: Path to the input file
        num_sequences: Number of sequences to read (None = read all)

    Returns:
        numpy array of shape (num_sequences, 178) with dtype uint8
    """
    with open(filename, 'rb') as f:
        data = f.read()

    # Calculate number of sequences
    total_sequences = len(data) // 178
    if num_sequences is None:
        num_sequences = total_sequences
    else:
        num_sequences = min(num_sequences, total_sequences)

    # Reshape into array of sequences
    sequences = np.frombuffer(data[:num_sequences * 178], dtype=np.uint8)
    sequences = sequences.reshape(num_sequences, 178)

    return sequences


def run_sliding_window(sequences, window_size=100, iterations=1):
    """
    Run sliding window diversity analysis on sequences.

    Args:
        sequences: numpy array of shape (num_sequences, 178) with dtype uint8
        window_size: Size of the sliding window
        iterations: Number of times to repeat the calculation (for benchmarking)

    Returns:
        numpy array of diversity values for each window position
    """
    num_sequences = sequences.shape[0]
    num_windows = num_sequences - window_size + 1

    if num_windows <= 0:
        raise ValueError(f"Not enough sequences ({num_sequences}) for window size {window_size}")

    # Ensure sequences is contiguous
    sequences = np.ascontiguousarray(sequences, dtype=np.uint8)

    # Call C function
    results_ptr = lib.run_sliding_window(
        sequences.ctypes.data_as(ctypes.POINTER(ctypes.c_char)),
        num_sequences,
        window_size,
        iterations
    )

    if not results_ptr:
        raise RuntimeError("Failed to run sliding window analysis")

    # Copy results to numpy array
    results = np.array([results_ptr[i] for i in range(num_windows)])

    # Free C memory
    lib.free_results(results_ptr)

    return results


def main():
    import sys
    import time

    if len(sys.argv) != 3:
        print(f"Usage: {sys.argv[0]} <num_sequences> <file.fa>", file=sys.stderr)
        sys.exit(1)

    num_sequences = int(sys.argv[1])
    filename = sys.argv[2]

    # Read sequences
    print(f"Reading {num_sequences} sequences from {filename}...")
    sequences = read_sequences(filename, num_sequences)
    print(f"Loaded {len(sequences)} sequences")

    # Run analysis with 1000 iterations for benchmarking
    window_size = 100
    iterations = 1000

    print(f"Running sliding window analysis (window_size={window_size}, iterations={iterations})...")
    start_time = time.time()
    results = run_sliding_window(sequences, window_size=window_size, iterations=iterations)
    elapsed = time.time() - start_time

    print(f"Completed in {elapsed:.3f} seconds")
    print(f"Results shape: {results.shape}")
    print(f"First 5 diversity values: {results[:5]}")
    print(f"Last 5 diversity values: {results[-5:]}")


if __name__ == "__main__":
    main()
