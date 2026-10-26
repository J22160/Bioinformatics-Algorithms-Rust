# Huffman Coding

## What Is This?

DNA sequences are stored as text files (FASTQ, FASTA) that can be
enormous — a human WGS dataset can be 200+ GB. **Huffman coding**
is a lossless compression algorithm that assigns shorter binary
codes to more frequent characters and longer codes to rare characters.

## How It Works

1. Count the frequency of each character in the data
2. Create a leaf node for each character
3. Repeatedly merge the two lowest-frequency nodes into a new parent
4. The path from root to each leaf defines the binary code
   (left = 0, right = 1)

## Example: Compressing a DNA Quality String

Quality string: `IIIIIFFFFFF!!!`

| Char | Freq | Huffman Code |
|------|------|-------------|
| F    | 6    | 0           |
| I    | 5    | 10          |
| !    | 3    | 11          |

Original: 14 chars × 8 bits = 112 bits
Huffman: 6×1 + 5×2 + 3×2 = 22 bits (80% compression)

## Complexity

- **Time**: O(n log n) using a min-heap
- **Space**: O(n) for the tree

## Bioinformatics Connection

Tools like `gzip` and `bgzf` (used in BAM files) employ Huffman coding
as part of the DEFLATE algorithm. Every time you work with compressed
BAM files, Huffman coding is at work behind the scenes.
