"""
Binary Golay encoder/decoder implementations.

This module implements three closely related binary Golay codes:
- Extended Golay G24:
    Parameters:
        (n, k, d) = (24, 12, 8)
    Properties:
        - N = 24: Codeword length.
        - K = 12: Number of information bits.
        - D = 8: Minimum Hamming distance.
        - Detects up to 3 bit errors
        - All codewords have even parity

- Perfect Golay G23:
    Parameters:
        (n, k, d) = (23, 12, 7)
    Properties:
        - N = 23: Codeword length.
        - K = 12: Number of information bits.
        - D = 7: Minimum Hamming distance.
        - Detects up to 3 bit errors
        - Obtained by shortening G24

- Shortened Golay G18:
    Parameters:
        (n, k, d) = (18, 6, 8)
    Properties:
        - N = 18: Codeword length.
        - K = 6: Number of information bits.
        - D = 8: Minimum Hamming distance.
        - Corrects up to 3 bit errors
        - Obtained by shortening G24

Key implementation details:
- All vectors are represented as Python lists of binary values:
    [0, 1, 1, 0, ...]
- Matrices are represented as nested lists and all matrix operations
  are performed over GF(2).

The G23 and G18 decoders reuse the G24 decoder through puncturing
and shortening transformations, avoiding duplication of the core
syndrome-decoding logic.

Terminology
    n      Codeword length.
    k      Number of information bits.
    d      Minimum Hamming distance.
    t      Error correcting capability.
    wt(v)  Hamming weight of vector v.
    GF(2)  Binary finite field with arithmetic modulo two.

Notes:
    - https://personalpages.manchester.ac.uk/staff/yuri.bazlov/code/notes/ch10.pdf
    - https://en.wikipedia.org/wiki/Binary_Golay_code
"""

from typing import override

from binarycorrectioncodes.common import BinaryVector, EncoderProtocol, Matrix, Scheme, Syndromes

type GolaySchemes = Golay_24_12_8 | Golay_23_12_7 | Golay_18_6_8


class GolayScheme(Scheme):
    @classmethod
    def _decode_syndrome(cls, vector: BinaryVector, stage: int = 1) -> BinaryVector:
        """
        Perform syndrome decoding for the extended Golay code.

        Mathematical Formulation:
            Given syndrome
                s = eHᵀ
            the goal is to recover an error vector e having
                wt(e) <= 3
            where wt() denotes the Hamming weight.

            Stage 1:
                Assume errors occurred primarily in the parity half
                of the codeword.

            Stage 2:
                Assume errors occurred primarily in the information
                half of the codeword.

            Let p_j denote the j-th row or column of the parity matrix P.

            The decoder searches for a low-weight error pattern by
            testing:
                wt(s) <= 3 or wt(s + p_j) <= 2
            which corresponds to a valid error configuration of total
            weight at most three bits.
        """
        hamming_weight = 2
        I_rows = cls.I.matrix
        # Base case: small syndrome
        if vector.weight() <= cls.T:
            if stage == 1:
                return [0] * cls.K + vector
            return vector + [0] * cls.K

        # Stage 1: try P columns
        if stage == 1:
            for j in range(cls.K):
                test = vector ^ cls.P.column(j)
                if test.weight() <= hamming_weight:
                    return I_rows[j] + test
            # recurse to stage 2
            s2 = vector * cls.P
            return cls._decode_syndrome(s2, stage=2)

        # Stage 2: try P rows
        for j in range(cls.K):
            test = vector ^ cls.P.matrix[j]
            if test.weight() <= hamming_weight:
                return test + I_rows[j]
        # failure case
        return BinaryVector([0] * cls.N)


class Golay_24_12_8(GolayScheme):  # Extended Decode
    """
    Extended Binary Golay Code G24.

    Parameters
    ----------
        - N = 24: Codeword length.
        - K = 12: Number of information bits.
        - D = 8: Minimum Hamming distance.
        - T = 3: Error correcting capability.

    Notes
    -----
    The extended binary Golay code is a linear block code with
    parameters
        (n, k, d) = (24, 12, 8)
    and therefore corrects
        t = floor((d - 1) / 2) = 3
    bit errors.

    The code is systematic and is generated using
        G = [ I12 | P ]
    where
        I12 : 12x12 identity matrix
        P   : 12x12 parity matrix
    Every codeword is computed as
        c = mG
    over GF(2).
    The parity-check matrix is
        H = [ Pᵀ | I12 ]
    and satisfies
        GHᵀ = 0.
    The syndrome of a received vector r is
        s = rHᵀ
    which depends only on the error pattern.

    Because the minimum distance is 8, the decoder can:
    - Correct up to 3 bit errors.
    - Detect up to 7 bit errors.
    """

    N = 24
    K = 12
    D = 8
    T = (D - 1) // 2

    I = Matrix.identity(K)
    # fmt: off
    P = Matrix([
        [1,1,0,1,1,1,0,0,0,1,0,1],
        [1,0,1,1,1,0,0,0,1,0,1,1],
        [0,1,1,1,0,0,0,1,0,1,1,1],
        [1,1,1,0,0,0,1,0,1,1,0,1],
        [1,1,0,0,0,1,0,1,1,0,1,1],
        [1,0,0,0,1,0,1,1,0,1,1,1],
        [0,0,0,1,0,1,1,0,1,1,1,1],
        [0,0,1,0,1,1,0,1,1,1,0,1],
        [0,1,0,1,1,0,1,1,1,0,0,1],
        [1,0,1,1,0,1,1,1,0,0,0,1],
        [0,1,1,0,1,1,1,0,0,0,1,1],
        [1,1,1,1,1,1,1,1,1,1,1,0]
    ])
    # fmt: on
    G = I | P
    H = P.transpose | I
    SYNDROME_TABLE = Syndromes.build_table(H, N, T)

    @override
    @classmethod
    def encode(cls, vector: BinaryVector) -> BinaryVector:
        """
        Encode a 12-bit message using the extended Golay code.

        Mathematical Formulation:
            Let m ∈ GF(2)^12 be the information vector.
            The encoded codeword is computed as
                c = mG
            where
                G = [ I12 | P ]
            is the systematic generator matrix.
            The resulting codeword
                c ∈ GF(2)^24
            consists of:
                c = [ message | parity ]
        """
        return super().encode(vector)[: cls.N]

    @override
    @classmethod
    def decode(cls, vector: BinaryVector) -> BinaryVector:
        """
        Decode a received 24-bit Golay codeword.

        Mathematical Formulation:
            The received vector is
                r = c + e
            where
                c = transmitted codeword
                e = error vector
            The syndrome is computed as
                s = rHᵀ
            Since every valid codeword satisfies
                cHᵀ = 0
            we obtain
                s = eHᵀ
            Therefore the syndrome depends only on the error pattern.

            The decoder applies the standard Golay syndrome decoding
            algorithm to determine the most likely error vector e and
            computes c = r + e over GF(2).
        """
        return super().decode(vector)[: cls.K]


class Golay_23_12_7(Golay_24_12_8):  # Standard Decode
    """
    Binary Golay Code G23.

    Parameters
    ----------
        - N = 23: Codeword length.
        - K = 12: Number of information bits.
        - D = 7: Minimum Hamming distance.
        - T = 3: Error correcting capability.

    Notes
    -----
    The standard binary Golay code is a perfect linear block code with
    parameters
        (n, k, d) = (23, 12, 7)
    and can therefore correct
        t = floor((d - 1) / 2) = 3
    bit errors.

    The code is obtained by shortning the extended Golay code G24,
    where the 24th parity coordinate is removed.

    Since G23 is a perfect code, every received vector is within
    Hamming distance <= 3 of exactly one valid codeword.
    """

    N = 23
    K = 12
    D = 7
    T = (D - 1) // 2

    I = Matrix.identity(K)
    P = Golay_24_12_8.P
    # Remove last column -> 23 columns
    G = (I | P).shorten_col(-1)
    H = (P.transpose | I).shorten_all(-1)
    SYNDROME_TABLE = Syndromes.build_table(H, N, T)

    @override
    @classmethod
    def encode(cls, vector: BinaryVector) -> BinaryVector:
        """
        Encode a 12-bit message using the binary Golay (23,12,7) code.

        Mathematical Formulation:
            Let m ∈ GF(2)^12 denote the information vector.
            Using the systematic generator matrix
                G = [ I12 | P ]
            the extended codeword is
                c24 = mG
            where
                c24 ∈ GF(2)^24.
        """
        if len(vector) != cls.K:
            raise ValueError("Message length mismatch")
        # Encode as extended Golay
        cw24 = Golay_24_12_8.encode(vector)
        # Remove overall parity bit
        return cw24[:-1]

    @override
    @classmethod
    def decode(cls, vector: BinaryVector) -> BinaryVector:
        """
        Decode a received Golay (23,12,7) codeword.

        Instead of implementing a dedicated G23 syndrome decoder,
        we reconstruct the missing parity coordinate, convert the
        received word back into a candidate G24 codeword, and reuse
        the extended Golay decoder.

        Mathematical Formulation:
            The received vector is
                r = c + e
            where
                c ∈ G23
                e = error vector
            and
                wt(e) <= 3
            for guaranteed correction.
            The parity bit is selected such that
                wt(c24) mod 2 = 0
            because every codeword in G24 has even parity.
        """
        if len(vector) != cls.N:
            raise ValueError("Codeword length mismatch")
        # reconstruct missing parity bit
        parity = sum(vector) & 1
        cw24 = BinaryVector([*vector, parity])
        decoded = Golay_24_12_8.decode(cw24)
        return decoded[: cls.K]


class Golay_18_6_8(Golay_24_12_8):  # Shortened Decode
    """
    Shortened Extended Golay Code.

    Parameters
    ----------
        - N = 18: Codeword length.
        - K = 6: Number of information bits.
        - D = 8: Minimum Hamming distance.
        - T = 3: Error correcting capability.

    The code is obtained by shortening the extended Golay code
    G24.

    Notes
    -----
    Starting from a (24,12,8) code:
        1. Fix six information bits to zero.
        2. Encode using G24.
        3. Remove those six known coordinates.

    This produces a shortened code with parameters (18,6,8)
    while preserving the minimum distance.

    Since d = 8, the decoder retains the ability to correct
        t = (8 - 1) // 2 = 3
    bit errors.
    """

    N = 18
    K = 6
    D = 8
    T = (D - 1) // 2

    # NOTE: We aren't shrinking the tables, instead work in the
    #   parent's table space
    SYNDROME_TABLE = Syndromes.build_table(Golay_24_12_8.H, Golay_24_12_8.N, T)

    @override
    @classmethod
    def encode(cls, vector: BinaryVector) -> BinaryVector:
        """
        Encode a 6-bit message using the shortened Golay code.

        Mathematical Formulation:
            Let m ∈ GF(2)^6 be the information vector.
            Construct an equivalent 12-bit message
                m' = [0,0,0,0,0,0,m]
            and encode it using the extended Golay encoder:
                c24 = m'G24
            Since the first six coordinates are known to be zero,
            they are removed to form the transmitted codeword:
                c18 = shorten(c24)
        """
        if len(vector) != cls.K:
            raise ValueError("Message length mismatch")
        # Shortening:
        # [000000 | m0..m5] -> 12-bit message
        full_msg = BinaryVector([0] * 6) + vector
        # Encode using (24,12,8)
        full_codeword = Golay_24_12_8.encode(full_msg)
        # Remove the first 6 positions
        return full_codeword[6:]

    @override
    @classmethod
    def decode(cls, vector: BinaryVector) -> BinaryVector:
        """
        Decode an 18-bit shortened Golay codeword.

        The shortened code is decoded by embedding the received
        word back into the parent (24,12,8) code space.

        Mathematical Formulation:
            Given r18 ∈ GF(2)^18 reconstruct
                r24 = [0,0,0,0,0,0,r18]
            and apply the extended Golay decoder.
            The corrected 12-bit message is
                m' = decode_G24(r24)
            and the original information bits are recovered as
                m = m'[6:]
        """
        if len(vector) != cls.N:
            raise ValueError("Codeword length mismatch")
        # Reinsert the shortened coordinates
        full_codeword = BinaryVector([0] * 6) + vector
        # Decode as a (24,12,8) codeword
        decoded_msg = Golay_24_12_8.decode(full_codeword)
        # Extract original 6 information bits
        return decoded_msg[6:]


class Golay(EncoderProtocol[GolaySchemes]):
    """
    High-level interface for Golay-Codes encoding/decoding.

    Handles:
        - Symbol-level encoding/decoding
        - Stream-to-symbol conversion via HexBits
    """
