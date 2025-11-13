from setuptools import setup
from setuptools.command.build_py import build_py
import subprocess
import os
import sys

class CustomBuild(build_py):
    """Custom build to compile shared library with specific flags"""

    def run(self):
        # Run the standard build first
        super().run()

        # Get the output directory
        package_dir = os.path.join(self.build_lib, 'kmer_variance')
        os.makedirs(package_dir, exist_ok=True)

        # Source and output paths
        source_file = os.path.join('kmer_variance', 'kmer_variance_lib.c')
        output_file = os.path.join(package_dir, 'libkmer_variance.so')

        # Compiler and flags
        compiler = os.environ.get('CC', 'clang')
        flags = ['-O3', '-march=native', '-fPIC', '-shared']

        # Compile command
        cmd = [compiler] + flags + ['-o', output_file, source_file]

        print(f"Compiling: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode != 0:
            print(result.stderr, file=sys.stderr)
            raise RuntimeError(f"Failed to compile {source_file}")

        print(f"Successfully compiled {output_file}")


# Read README if it exists
readme = ""
if os.path.exists("README.md"):
    with open("README.md", "r", encoding="utf-8") as f:
        readme = f.read()

setup(
    name="kmer-variance",
    version="0.1.0",
    author="Your Name",
    author_email="your.email@example.com",
    description="Fast sliding window kmer diversity analysis using incremental Hamming distances",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/kmer-variance",
    packages=["kmer_variance"],
    package_data={
        "kmer_variance": ["*.c", "*.so"],
    },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.7",
    install_requires=[
        "numpy>=1.19.0",
    ],
    cmdclass={
        'build_py': CustomBuild,
    },
    zip_safe=False,
)
