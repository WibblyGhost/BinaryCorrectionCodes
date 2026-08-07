#  Forward Error Correction (FEC)

A pure Python implementation of several Forward Error Correction (FEC) schemes used throughout APCO P25 and related digital radio protocols.

## Features

This library currently provides implementations of:

- Binary Hamming Codes
    - Hamming (15,11,3)
    - Shortened Hamming (10,6,3)

- Binary Golay Codes
    - Extended Golay (24,12,8)
    - Perfect Golay (23,12,7)
    - Shortened Golay (18,6,8)

- Binary Cyclic / BCH Codes
    - Shortened BCH (16,8,5)

- Reed-Solomon Codes over GF(2⁶)
    - RS(36,20,17)
    - RS(24,16,9)
    - RS(24,12,13)

## Overview

All implementations are written in pure Python and provide:

- Encoding
- Decoding
- Error detection
- Error correction
- Syndrome-based decoding where applicable
- Support for stream and symbol level operations

The implementations are designed around the requirements of APCO P25 Fixed Station Interface (FSI) messaging and related protocol structures.

---

# Binary Hamming Codes

Implemented Schemes:

| Code | N | K | D | Correction |
|--------|----|----|----|------------|
| Hamming(15,11,3) | 15 | 11 | 3 | 1 bit |
| Hamming(10,6,3) | 10 | 6 | 3 | 1 bit |

## Properties

### Hamming (15,11,3)

- Perfect single-error-correcting code
- Corrects up to 1 bit error
- Detects up to 2 bit errors
- Uses 4 parity bits

### Shortened Hamming (10,6,3)

- Derived from Hamming (15,11,3)
- Corrects up to 1 bit error
- Detects up to 2 bit errors

## Implementation

The Hamming implementation uses:

- Systematic generator matrices
- Syndrome decoding
- Parity-check matrices
- GF(2) arithmetic

Parity bits occupy positions:

```text
1, 2, 4, 8
```

and the syndrome directly identifies the location of a single-bit error.

---

# Binary Golay Codes

Implemented Schemes:

| Code | N | K | D | Correction |
|--------|----|----|----|------------|
| Golay(24,12,8) | 24 | 12 | 8 | 3 bits |
| Golay(23,12,7) | 23 | 12 | 7 | 3 bits |
| Golay(18,6,8) | 18 | 6 | 8 | 3 bits |

## Properties

### Extended Golay (24,12,8)

- Extended binary Golay code
- Minimum distance 8
- Detects up to 7 errors
- Corrects up to 3 errors
- All codewords have even parity

### Perfect Golay (23,12,7)

- Obtained by puncturing the extended code
- Perfect error-correcting code
- Corrects up to 3 errors

### Shortened Golay (18,6,8)

- Derived from Golay (24,12,8)
- Corrects up to 3 errors
- Reuses the parent Golay decoder

## Implementation

Golay decoding is performed using syndrome decoding.

The shortened and punctured variants reuse the G24 decoder through shortening and puncturing transformations, reducing implementation complexity while preserving correctness.

---

# Binary BCH / Cyclic Codes

Implemented Schemes:

| Code | N | K | D | Correction |
|--------|----|----|----|------------|
| BCH(16,8,5) | 16 | 8 | 5 | 2 bits |

## Properties

- Shortened BCH code
- Derived from BCH(31,23,5)
- Minimum Hamming distance of 5
- Corrects up to 2 bit errors
- Detects up to 4 bit errors

## Implementation

The encoder performs:

```text
Message Polynomial
        ↓
Polynomial Division
        ↓
Parity Bits
        ↓
Codeword
```

Decoding uses:

```text
Received Word
        ↓
Syndrome Calculation
        ↓
Lookup Table
        ↓
Error Pattern
        ↓
Correction
```

A syndrome lookup table is generated for all correctable error patterns of weight ≤ 2.

---

# Reed-Solomon Codes

Implemented Schemes:

| Code | Data Symbols | Parity Symbols | Total Symbols |
|--------|--------|--------|--------|
| RS(36,20,17) | 20 | 16 | 36 |
| RS(24,16,9) | 16 | 8 | 24 |
| RS(24,12,13) | 12 | 12 | 24 |

Symbol size:

```text
6 bits
```

Field:

```text
GF(2^6)
```

## Features

- Systematic encoding
- Symbol-level error correction
- Finite field arithmetic using lookup tables
- Efficient decoding algorithms

## Decoder Pipeline

### 1. Syndrome Calculation

Compute syndromes from the received codeword.

### 2. Berlekamp-Massey

Determine the error locator polynomial.

### 3. Chien Search

Locate symbol error positions.

### 4. Forney Algorithm

Calculate error magnitudes.

### 5. Correction

Apply the calculated corrections to recover the original message.

---

# Mathematical Foundations

The implementations use arithmetic over finite fields.

## Binary Codes

Hamming, Golay and BCH codes operate over:

```text
GF(2)
```

where:

```text
0 + 0 = 0
0 + 1 = 1
1 + 1 = 0
```

(addition modulo 2)

## Reed-Solomon Codes

The Reed-Solomon implementation operates over:

```text
GF(2^6)
```

allowing symbol-based correction of 6-bit values.

---

# Common Design

All schemes expose a common interface:

```python
encoded = Scheme.encode(message)
decoded = Scheme.decode(codeword)
```

and support:

- Binary vectors
- Bit streams
- Systematic encoding
- Syndrome-based error correction

---

# Applications

These codes are commonly used throughout digital radio systems and P25 infrastructure for:

- Link Control Data
- Voice Header protection
- Encryption Sync Words
- Low-Speed Data
- Fixed Station Interface messaging
- General digital communications

---

# References

## Hamming Codes

- https://en.wikipedia.org/wiki/Hamming_code
- https://en.wikipedia.org/wiki/Syndrome_decoding

## Golay Codes

- https://en.wikipedia.org/wiki/Binary_Golay_code
- https://personalpages.manchester.ac.uk/staff/yuri.bazlov/code/notes/ch10.pdf

## BCH Codes

- https://en.wikipedia.org/wiki/BCH_code
- https://en.wikipedia.org/wiki/Cyclic_code
- https://web.ntpu.edu.tw/~yshan/BCH_code.pdf

## Reed-Solomon Codes

- https://en.wikiversity.org/wiki/Reed%E2%80%93Solomon_codes_for_coders
- https://en.wikipedia.org/wiki/Berlekamp%E2%80%93Massey_algorithm
- https://en.wikipedia.org/wiki/Chien_search
- https://en.wikipedia.org/wiki/Forney_algorithm


# Development

### Linting and Code Formatters

For all projects uploaded to GitLab we highly recommend running code linters and code formatters to keep all code compliant with the PEP-8 standards and to unify codebases.

For any **Python 3.x** project we require that you use Ruff for formatting and linting.

When uploading any code changes to this repo please run the code linters and formatters before commiting any code.

**Please install `pre-commit` via the command `pre-commit install && pre-commit autoupdate`**

### Versioning Standards & PyPi Publishing

When you create git branches, every commit will run unit tests and build jobs to assure that no ground breaking errors have been uploaded. We require a merge request to merge anything into the master branch.

To trigger a PyPi publishing job you firstly need to create a tag through GitLab, this will trigger an automated build procedure.

**Note:** Please look at the following versioning guidelines to help with tagging versions:

[Semantic Versioning Guidelines](https://semver.org/)

### General Versioning Guidelines

Given a version number `MAJOR.MINOR.PATCH` (e.g. 1.0.0), increment the:

- `MAJOR` version when you make incompatible API changes
- `MINOR` version when you add functionality in a backwards compatible manner
- `PATCH` version when you make backwards compatible bug fixes
