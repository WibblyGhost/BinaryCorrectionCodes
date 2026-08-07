"""
Binary Hamming encoder/decoder implementations.

This module implements the following binary Hamming codes:

- Hamming H15:
  Parameters:
    (n, k, d) = (15, 11, 3)

  Properties:
    - N = 15 : Codeword length
    - K = 11 : Number of information bits
    - D = 3  : Minimum Hamming distance
    - Corrects up to 1 bit error
    - Detects up to 2 bit errors
    - Perfect single-error-correcting code

- Shortened Hamming H10:
  Parameters:
    (n, k, d) = (10, 6, 3)

  Properties:
    - N = 10 : Codeword length
    - K = 6  : Number of information bits
    - D = 3  : Minimum Hamming distance
    - Corrects up to 1 bit error
    - Detects up to 2 bit errors
    - Obtained by shortening H15

Key implementation details:
- All vectors are represented as Python lists of binary values:
    [0, 1, 1, 0, ...]
- Codewords are encoded in systematic form where the
  information bits appear unchanged in the codeword.
- Parity bits occupy positions that are powers of two:
      1, 2, 4, 8
- Decoding is performed using syndrome decoding.
- The syndrome directly identifies the location of a
  single-bit error and allows correction by flipping
  the corresponding bit.

Terminology
    n      Codeword length
    k      Number of information bits
    d      Minimum Hamming distance
    t      Error correcting capability
    wt(v)  Hamming weight of vector v
    GF(2)  Binary finite field

Notes:
    - https://en.wikipedia.org/wiki/Hamming_code
    - https://en.wikipedia.org/wiki/Syndrome_decoding
"""

from binarycorrectioncodes.common import EncoderProtocol, Matrix, Scheme, Syndromes

# (15,11,3) Hamming code words
# (10,6,3) shortened Hamming code

type HammingSchemes = Hamming_15_11_3 | Hamming_10_6_3


class HammingScheme(Scheme):
    """
    Protocol implemented by all Hamming code schemes.

    All codes are linear block codes over GF(2) and use
    systematic codewords consisting of information bits
    and parity bits.
    """


def _create_error_table(H: Matrix, N: int) -> dict[tuple[int, ...], int]:
    # NOTE: This is here, as dictionary comprehension doesn't
    #   work correcctly in class definitions
    return {tuple(H.column(i)): i for i in range(N)}


class Hamming_15_11_3(HammingScheme):
    """
    Standard binary Hamming code.

    Parameters
    ----------
        - N = 15: Codeword length.
        - K = 11: Number of information bits.
        - D = 3: Minimum Hamming distance.

    Properties:
        - Encodes 11 information bits
        - Produces a 15-bit codeword
        - Corrects 1 bit error
        - Detects 2 bit errors
        - Uses 4 parity bits

    Parity bit positions:
        1, 2, 4, 8

    Information bit positions:
        3, 5, 6, 7,
        9, 10, 11, 12,
        13, 14, 15
    """

    N = 15
    K = 11
    D = 3
    T = (D - 1) // 2

    H = Matrix(
        [
            [1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 0, 1, 0, 0, 0],
            [1, 0, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0],
            [0, 1, 1, 1, 0, 0, 1, 0, 0, 1, 1, 0, 0, 1, 0],
            [1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1],
        ]
    )
    P = Matrix(
        [
            [1, 1, 0, 1],
            [1, 0, 1, 1],
            [0, 1, 1, 1],
            [1, 1, 1, 0],
            [1, 0, 0, 0],
            [0, 1, 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1],
            [1, 1, 0, 0],
            [1, 0, 1, 0],
            [0, 1, 1, 0],
        ]
    )
    G = Matrix.identity(11) | P
    H = P.transpose | Matrix.identity(4)
    SYNDROME_TABLE = Syndromes.build_table(H, N, T)


class Hamming_10_6_3(Hamming_15_11_3):
    """
    Shortened binary Hamming code.

    Parameters
    ----------
        - N = 10: Codeword length.
        - K = 6: Number of information bits.
        - D = 3: Minimum Hamming distance.

    Properties:
        - Encodes 6 information bits
        - Produces a 10-bit codeword
        - Corrects 1 bit error
        - Detects 2 bit errors
        - Derived from Hamming (15,11,3)

    Notes
    -----
    This implementation reuses the parent
    Hamming (15,11,3) encoder and decoder.

    Encoding is performed by padding the
    message with fixed zero bits before
    applying the parent code.

    Decoding performs the inverse operation.
    """

    N = 10
    K = 6
    D = 3
    T = (D - 1) // 2

    P = Matrix([[1, 1, 0, 1], [1, 0, 1, 1], [0, 1, 1, 1], [1, 1, 1, 0], [1, 1, 0, 0], [1, 0, 1, 0]])
    I = Matrix.identity(6)
    G = I | P
    H = P.transpose | Matrix.identity(4)
    SYNDROME_TABLE = Syndromes.build_table(H, N, 1)


class Hamming(EncoderProtocol[HammingSchemes]):
    """
    High-level interface for Hamming code encoding and decoding.

    Provides a common API consistent with the
    Golay and Reed-Solomon implementations.

    # NOTE: Hamming can only correct a single bit error
    """
