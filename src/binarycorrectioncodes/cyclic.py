"""
Binary cyclic (shortened BCH) encoder/decoder implementations.

This module implements:

- BCH(16,8,5)

Parameters
----------
   (n, k, d) = (16, 8, 5)

    Properties:
        - N = 16 : Codeword length
        - K = 8  : Number of information bits
        - D = 5  : Minimum Hamming distance
        - Corrects up to 2 bit errors
        - Detects up to 4 bit errors
        - Shortened from BCH(31,23,5)

Implementation Details
- The encoder uses polynomial division over GF(2).
- A received codeword is decoded using syndrome
  decoding. A lookup table maps syndromes onto
  all correctable error patterns of weight <= 2.

Terminology
    n      Codeword length
    k      Information bits
    d      Minimum distance
    t      Error correction capability
    GF(2)  Binary field

Notes
-----
    - https://en.wikipedia.org/wiki/Cyclic_code
    - https://en.wikipedia.org/wiki/BCH_code
    - https://web.ntpu.edu.tw/~yshan/BCH_code.pdf
"""

from binarycorrectioncodes.common import EncoderProtocol, Matrix, Scheme, Syndromes

type CyclicSchemes = BCH_16_8_5


class CyclicScheme(Scheme): ...


class BCH_16_8_5(CyclicScheme):
    """
    Shortened BCH code.

    Parameters
    ----------
        - N = 16: Codeword length.
        - K = 8: Number of information bits.
        - D = 5: Minimum Hamming distance.

    Properties:
        - Encodes 8 information bits
        - Produces a 16-bit codeword
        - Corrects 2 bit error
        - Detects 3 bit errors
        - Uses 4 parity bits
    """

    N = 16
    K = 8
    D = 5
    T = (D - 2) // 2

    I = Matrix.identity(K)
    P = Matrix(
        [
            [1, 1, 1, 0, 1, 0, 0, 0],
            [0, 1, 1, 1, 0, 1, 0, 0],
            [0, 0, 1, 1, 1, 0, 1, 0],
            [0, 0, 0, 1, 1, 1, 0, 1],
            [1, 0, 0, 0, 1, 1, 1, 0],
            [0, 1, 0, 0, 0, 1, 1, 1],
            [1, 0, 1, 0, 0, 0, 1, 1],
            [1, 1, 0, 1, 0, 0, 0, 1],
        ]
    )
    G = I | P
    H = P.transpose | Matrix.identity(8)
    SYNDROME_TABLE = Syndromes.build_table(H, N, T)


class Cyclic(EncoderProtocol[CyclicSchemes]):
    """
    High-level interface for Hamming code encoding and decoding.

    Provides a common API consistent with the
    Golay and Reed-Solomon implementations.

    NOTE: Cyclic can only correct a single/dual bit errors
    """
