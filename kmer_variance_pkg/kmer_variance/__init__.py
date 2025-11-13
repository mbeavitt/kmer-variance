"""
kmer_variance - Fast sliding window kmer diversity analysis

Example usage:
    import numpy as np
    from kmer_variance import calculate_diversity

    # sequences: numpy array of shape (num_sequences, 178) with dtype uint8
    sequences = np.random.randint(65, 85, (1000, 178), dtype=np.uint8)

    # Run analysis with window size 100
    diversity = calculate_diversity(sequences, window_size=100)
"""

import ctypes
import numpy as np
from pathlib import Path

__version__ = "0.1.0"

# Load the shared library
_lib_path = Path(__file__).parent / "libkmer_variance.so"
_lib = ctypes.CDLL(str(_lib_path))

# Define function signatures
_lib.run_sliding_window.argtypes = [
    ctypes.POINTER(ctypes.c_char),  # sequences
    ctypes.c_int,                    # num_sequences
    ctypes.c_int,                    # window_size
    ctypes.c_int                     # iterations
]
_lib.run_sliding_window.restype = ctypes.POINTER(ctypes.c_double)

_lib.free_results.argtypes = [ctypes.POINTER(ctypes.c_double)]
_lib.free_results.restype = None


def calculate_diversity(sequences, window_size=100):
    """
    Calculate kmer diversity across a sliding window of sequences.

    This function computes the average pairwise Hamming distance of 4-mer
    presence/absence patterns across a sliding window. The algorithm uses
    incremental updates for efficiency: it computes the first window fully,
    then slides by adding/removing edge distances.

    Args:
        sequences: numpy array of shape (num_sequences, 178) with dtype uint8
                   Each row is a 178-byte sequence
        window_size: Size of the sliding window (default: 100)
                    Must be > 1 and <= num_sequences

    Returns:
        numpy array of diversity values, one per window position
        Shape: (num_sequences - window_size + 1,)
        Values are normalized to [0, 1] range

    Raises:
        ValueError: If sequences has wrong shape or window_size is invalid
        TypeError: If sequences is not a numpy array
        RuntimeError: If the C library fails

    Example:
        >>> import numpy as np
        >>> from kmer_variance import calculate_diversity
        >>>
        >>> # Create random sequences (A, C, G, T as ASCII)
        >>> sequences = np.random.randint(65, 85, (1000, 178), dtype=np.uint8)
        >>>
        >>> # Calculate diversity
        >>> diversity = calculate_diversity(sequences, window_size=100)
        >>> print(f"Diversity shape: {diversity.shape}")
        >>> print(f"Mean diversity: {diversity.mean():.4f}")
    """
    # Validate input
    if not isinstance(sequences, np.ndarray):
        raise TypeError("sequences must be a numpy array")

    if sequences.ndim != 2:
        raise ValueError(f"sequences must be 2D, got shape {sequences.shape}")

    if sequences.shape[1] != 178:
        raise ValueError(f"sequences must have 178 columns, got {sequences.shape[1]}")

    num_sequences = sequences.shape[0]

    if window_size < 2:
        raise ValueError(f"window_size must be >= 2, got {window_size}")

    if window_size > num_sequences:
        raise ValueError(
            f"window_size ({window_size}) cannot be larger than "
            f"num_sequences ({num_sequences})"
        )

    num_windows = num_sequences - window_size + 1

    # Ensure sequences is contiguous uint8
    sequences = np.ascontiguousarray(sequences, dtype=np.uint8)

    # Call C function (iterations=1 for single calculation)
    results_ptr = _lib.run_sliding_window(
        sequences.ctypes.data_as(ctypes.POINTER(ctypes.c_char)),
        num_sequences,
        window_size,
        1  # single iteration
    )

    if not results_ptr:
        raise RuntimeError("C library failed to compute sliding window diversity")

    # Copy results to numpy array
    results = np.array([results_ptr[i] for i in range(num_windows)], dtype=np.float64)

    # Free C memory
    _lib.free_results(results_ptr)

    return results


__all__ = ['calculate_diversity']
