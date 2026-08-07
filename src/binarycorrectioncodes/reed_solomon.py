"""
Reed-Solomon encoder/decoder over GF(2^6).

This module implements the following binary Reed-Solomon codes:
    - RS(36, 20, 17)
    - RS(24, 16, 9)
    - RS(24, 12, 13)

This module implements systematic Reed-Solomon (RS) coding using:
- Finite field GF(2^6) with log/antilog tables
- Generator polynomial construction
- Encoding via polynomial division
- Decoding via:
    * Syndrome calculation
    * Berlekamp-Massey algorithm (error locator)
    * Chien search (error positions)
    * Forney algorithm (error magnitudes)

Key implementation details:
- All polynomials are stored in ascending order:
  [a0, a1, a2, ...] = a0 + a1·x + a2·x² + ...
- Field arithmetic uses lookup tables for efficiency
- Symbols are 6-bit values (0-63)

This implementation supports:
- RS(36,20), RS(24,16), RS(24,12)
- Correction of up to nsym//2 symbol errors

Notes:
    Reed-Solomon:
        - https://en.wikiversity.org/wiki/Reed%E2%80%93Solomon_codes_for_coders#Syndrome_calculation
    Horner:
        - https://www.academia.edu/110249244/A_Note_on_Horners_Method
        - https://en.wikipedia.org/wiki/Horner%27s_method
    Berlekamp-Massey
        - https://arxiv.org/pdf/2211.11721
        - https://en.wikipedia.org/wiki/Berlekamp%E2%80%93Massey_algorithm
    Chien
        - https://arxiv.org/pdf/cs/0606035
        - https://en.wikipedia.org/wiki/Chien_search
    Forney
        - https://every-algorithm.github.io/2024/05/01/forney_algorithm.html
        - https://arxiv.org/pdf/2108.02758
        - https://en.wikipedia.org/wiki/Forney_algorithm
"""

from collections.abc import Iterable, Sequence
from typing import ClassVar, Protocol

from bitstructures import BitStream

from binarycorrectioncodes.common import HexBits

type RsSchemes = RS_36_20_17 | RS_24_16_9 | RS_24_12_13

# IMPORTANT: All polynomials are stored in ascending order
# [a0, a1, a2, ...] = a0 + a1 x + a2 x² + ...


class ReedSolomonScheme(Protocol):
    """
    Base class for Reed-Solomon parameter sets.

    Attributes:
        k     : number of information symbols
        nsym  : number of parity symbols

    Codeword length:
        n = k + nsym
    """

    K: int
    nsym: int

    @classmethod
    def encode(cls, symbols: HexBits) -> HexBits:
        return _RsMatrix.encode_poly(symbols, cls.nsym)

    @classmethod
    def decode(cls, symbols: HexBits) -> HexBits:
        return _RsMatrix.decode_poly(symbols, cls.nsym, cls.K)


# ---------------------------
# GF(2^6) arithmetic
# ---------------------------


class GF6:
    """
    Finite field arithmetic for GF(2^6).

    Elements are represented as integers in [0, 63].

    Arithmetic is performed using logarithm (LOG_TABLE) and exponent (EXP_TABLE)
    lookup tables for fast multiplication and division.

    Field properties:
        - Addition: bitwise XOR
        - Multiplication: log/exp modular arithmetic
        - Inverses: a⁻¹ = α^(−log(a))

    The field is generated using a primitive polynomial (implicit in tables).

    MAX_EXP = 63 corresponds to the multiplicative order of the field:
        α^63 = 1
    """

    LOG_0 = 0xFF
    MAX_EXP = 63
    # fmt: off
    EXP_TABLE: ClassVar = [
        0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x03, 0x06,
        0x0c, 0x18, 0x30, 0x23, 0x05, 0x0a, 0x14, 0x28,
        0x13, 0x26, 0x0f, 0x1e, 0x3c, 0x3b, 0x35, 0x29,
        0x11, 0x22, 0x07, 0x0e, 0x1c, 0x38, 0x33, 0x25,
        0x09, 0x12, 0x24, 0x0b, 0x16, 0x2c, 0x1b, 0x36,
        0x2f, 0x1d, 0x3a, 0x37, 0x2d, 0x19, 0x32, 0x27,
        0x0d ,0x1a, 0x34, 0x2b, 0x15, 0x2a, 0x17, 0x2e,
        0x1f, 0x3e, 0x3f, 0x3d, 0x39, 0x31, 0x21, 0x01
    ]
    LOG_TABLE: ClassVar = [
        LOG_0, 0,  1,  6,   2, 12,  7, 26,
        3,     32, 13, 35,  8, 48, 27, 18,
        4,     24, 33, 16, 14, 52, 36, 54,
        9,     45, 49, 38, 28, 41, 19, 56,
        5,     62, 25, 11, 34, 31, 17, 47,
        15,    23, 53, 51, 37, 44, 55, 40,
        10,    61, 46, 30, 50, 22, 39, 43,
        29,    60, 42, 21, 20, 59, 57, 58
    ]
    # fmt: on

    # Mathematical Field Fundamentals

    @classmethod
    def gf_add(cls, a: int, b: int) -> int:
        """
        Add two field elements.

        In GF(2^m), addition is XOR.

        This corresponds to polynomial addition modulo 2.
        """
        return a ^ b

    @classmethod
    def gf_mul(cls, a: int, b: int) -> int:
        """
        Multiply two field elements using log/exp tables.

        Formula:
            a * b = α^(log(a) + log(b) mod 63)

        Special cases:
            - if a == 0 or b == 0 → 0
        """
        if a == 0 or b == 0:
            return 0
        return cls.EXP_TABLE[(cls.LOG_TABLE[a] + cls.LOG_TABLE[b]) % cls.MAX_EXP]

    @classmethod
    def gf_div(cls, a: int, b: int) -> int:
        """
        Divide two field elements.

        Formula:
            a / b = α^(log(a) - log(b) mod 63)
        """
        if b == 0:
            raise ZeroDivisionError()
        if a == 0:
            return 0
        return cls.EXP_TABLE[(cls.LOG_TABLE[a] - cls.LOG_TABLE[b]) % cls.MAX_EXP]

    @classmethod
    def gf_pow(cls, a: int, power: int) -> int:
        """
        Raise a field element to an integer power.

        Formula:
            a^p = α^(log(a) * p mod 63)

        Special cases:
            - 0^p = 0 for p > 0
            - 0^0 = 1
        """
        if a == 0:
            return 0 if power > 0 else 1
        return cls.EXP_TABLE[(cls.LOG_TABLE[a] * power) % cls.MAX_EXP]

    @classmethod
    def gf_inv(cls, a: int) -> int:
        """
        Multiplicative inverse of a field element.

        Formula:
            a⁻¹ = α^(63 - log(a))
        """
        if a == 0:
            raise ZeroDivisionError()
        return cls.EXP_TABLE[(cls.MAX_EXP - cls.LOG_TABLE[a]) % cls.MAX_EXP]


class _RsMatrix:
    # Generate Matricies

    @staticmethod
    def poly_mul(p: list[int], q: list[int]) -> list[int]:
        """
        Multiply two polynomials over GF(2^6).

        Inputs:
            p, q: polynomials in ascending order

        Output:
            Resulting polynomial (ascending)

        Operation:
            Standard convolution:
                (p * q)[k] = Σ p[i] * q[k - i]

        All multiplications and additions are done in GF(2^6).
        """
        res = [0] * (len(p) + len(q) - 1)
        for i in range(len(p)):
            for j in range(len(q)):
                res[i + j] = GF6.gf_add(res[i + j], GF6.gf_mul(p[i], q[j]))
        return res

    @classmethod
    def rs_generator_poly(cls, nsym: int) -> list[int]:
        """
        Construct the Reed-Solomon generator polynomial.

        g(x) = (x - α¹)(x - α²)...(x - α^nsym)

        This defines the code's error-correcting capability.

        nsym:
            Number of parity symbols (n - k)

        Returns:
            Generator polynomial in ascending form

        Notes:
            Using consecutive roots ensures:
            - Maximum distance separable (MDS)
            - Can correct up to nsym//2 symbol errors
        """
        g = [1]
        for i in range(nsym):
            g = cls.poly_mul(g, [1, GF6.gf_pow(2, i + 1)])
        return g

    # Encode
    @classmethod
    def encode_poly(cls, symbols: HexBits, nsym: int) -> HexBits:
        """
        Encode a message using systematic Reed-Solomon encoding.

        Process:
            1. Append nsym zeros to message
            2. Divide by generator polynomial
            3. Replace zeros with remainder

        Result:
            codeword = message || parity

        Returns:
            Encoded symbol sequence of length k + nsym
        """
        gen = cls.rs_generator_poly(nsym)
        msg = symbols + [0] * nsym

        for i in range(len(symbols)):
            coef = msg[i]
            if coef != 0:
                for j in range(len(gen)):
                    msg[i + j] = GF6.gf_add(msg[i + j], GF6.gf_mul(coef, gen[j]))
        # remainder is last nsym symbols
        return symbols + msg[-nsym:]

    # Decode

    # Error Detection & Error Correction

    @staticmethod
    def horner(poly: Sequence[int], x: int) -> int:
        """
        Evaluate a polynomial at a given field element using Horner's method.

        This implementation assumes a *descending evaluation order*, even though
        polynomials are stored in ascending form. All other formulas expect
        the polynomials of an ascending order.

        Given:
            poly = [a0, a1, a2, ..., an]

        It evaluates the polynomial as:
            a0·x^(n-1) + a1·x^(n-2) + ... + an

        using Horner's scheme:
            y = a0
            y = y·x + a1
            y = y·x + a2
            ...
            y = y·x + an

        Parameters
        ----------
            poly : list[int]
                Polynomial coefficients in ascending storage:
                    [a0, a1, a2, ...] = a0 + a1·x + a2·x² + ...
                BUT interpreted here in reverse (descending) order.

            x : int
                Field element in GF(2^6) at which to evaluate the polynomial.

        Returns
        -------
            int : The value of the polynomial evaluated at x in GF(2^6).

        Notes
        -----
            - All arithmetic is performed in GF(2^6)
            - This function is primarily used in:
                * syndrome calculation
                * Chien search
            - Forney uses explicit ascending evaluation instead to avoid
            ambiguity from this reversed interpretation.

        Warning:
            This function mixes ascending storage with descending evaluation.
            If you change polynomial conventions elsewhere in the code,
            this function must be updated accordingly.
        """
        # Horners' method assumes all symbols are in transmission order
        # Symbol representation | Order
        # MSB-first,            | polynomial
        # LSB-first             | polynomial reverse
        y = poly[0]
        for i in range(1, len(poly)):
            y = GF6.gf_mul(y, x) ^ poly[i]
        return y

    @classmethod
    def calc_syndromes(cls, symbols: HexBits, nsym: int) -> HexBits:
        """
        Compute syndromes of the received codeword.

        S[i] = r(α^(i+1))

        If all syndromes are zero → no errors.

        Interpretation:
            Evaluates received polynomial at consecutive powers
            of the primitive element α.

        Returns:
            List of length nsym
        """
        # S[i] = r(α^(i+1))
        synd = HexBits([0] * nsym)
        for i in range(nsym):
            # α = EXP_TABLE[1]
            alpha_i = GF6.EXP_TABLE[i + 1]
            synd[i] = cls.horner(symbols, alpha_i)
        return synd

    @staticmethod
    def berlekamp_massey(synd: HexBits) -> list[int]:
        """
        Compute the error locator polynomial Λ(x).

        Uses the Berlekamp-Massey algorithm to find the shortest
        linear feedback shift register (LFSR) that generates the syndromes.

        Output:
            Λ(x) = 1 + λ1 x + λ2 x² + ...

        Degree of Λ(x) = number of errors
        """
        # Λ(x) = [1, λ1, λ2, ...]
        C = [1] + [0] * len(synd)
        B = [1] + [0] * len(synd)

        L = 0
        m = 1
        b = 1

        for n in range(len(synd)):
            d = synd[n]
            for i in range(1, L + 1):
                d = GF6.gf_add(d, GF6.gf_mul(C[i], synd[n - i]))
            if d == 0:
                m += 1
            else:
                T = C[:]
                coef = GF6.gf_div(d, b)

                for i in range(len(B)):
                    if B[i] != 0:
                        C[i + m] = GF6.gf_add(C[i + m], GF6.gf_mul(coef, B[i]))

                if n < 2 * L:
                    m += 1
                else:
                    L = n + 1 - L
                    B = T
                    b = d
                    m = 1
        return C[: L + 1]

    @classmethod
    def chien_search(cls, err_loc: list[int], n: int) -> list[int]:
        """
        Find error locations using Chien search.

        Evaluates Λ(x) over all field elements:

            Λ(α^i) = 0  ⇒ error at position corresponding to i

        Efficient root-finding over finite field.

        Returns:
            List of error positions in the symbol array
        """
        err_pos = []
        for i in range(n):
            x = GF6.EXP_TABLE[i % GF6.MAX_EXP]  # α^-i
            if cls.horner(err_loc, x) == 0:
                err_pos.append(n - 1 - i)
        return err_pos

    @classmethod
    def forney(
        cls, symbols: HexBits, synd: HexBits, err_pos: list[int], err_loc: list[int]
    ) -> HexBits:
        """
        Compute and apply error magnitudes using Forney's algorithm.

        Steps:
            1. Construct syndrome polynomial S(x)
            2. Compute error evaluator polynomial:
                Ω(x) = (Λ(x) * S(x)) mod x^nsym
            3. Compute derivative Λ'(x)
            4. For each error location i:
                error = Ω(α^(-i)) / Λ'(α^(-i))
            5. Apply correction via XOR

        Notes:
            - Evaluation uses α^(-i) (reciprocal root)
            - Works for multiple errors
            - Requires consistent polynomial ordering

        Returns:
            Corrected symbol sequence
        """
        n = len(symbols)
        # ensure consistent polynomial orientation
        synd_poly = [0, *synd]

        # compute Ω(x) = (S * Λ) mod x^nsym
        err_eval = [0] * (len(err_loc) + len(synd_poly))
        for i in range(len(err_loc)):
            for j in range(len(synd_poly)):
                err_eval[i + j] ^= GF6.gf_mul(err_loc[i], synd_poly[j])

        err_eval = err_eval[: len(err_loc)]

        # derivative Λ'(x)
        err_loc_deriv = []
        deg = len(err_loc) - 1
        for i in range(deg):
            power = deg - i  # actual power of this coefficient
            if power % 2 == 1:
                err_loc_deriv.append(err_loc[i])
            else:
                err_loc_deriv.append(0)

        corrected = symbols[:]
        for pos in err_pos:
            x_inv = GF6.EXP_TABLE[(63 - (n - 1 - pos)) % 63]
            # evaluate Ω(x_inv)
            num = 0
            for i in range(len(err_eval)):
                num ^= GF6.gf_mul(err_eval[i], GF6.gf_pow(x_inv, i))
            # evaluate Λ'(x_inv)
            den = 0
            for i in range(len(err_loc_deriv)):
                den ^= GF6.gf_mul(err_loc_deriv[i], GF6.gf_pow(x_inv, i))
            err = GF6.gf_div(num, den)
            corrected[pos] ^= err
        return corrected

    @classmethod
    def decode_poly(cls, symbols: HexBits, nsym: int, k: int) -> HexBits:
        """
        Decode a Reed-Solomon codeword.

        Pipeline:
            1. Compute syndromes
            2. If all zero → no errors
            3. Compute error locator polynomial (Berlekamp-Massey)
            4. Locate errors (Chien search)
            5. Correct errors (Forney)
            6. Verify correction

        Returns:
            Original message (first k symbols)

        """
        synd = cls.calc_syndromes(symbols, nsym)
        if all(s == 0 for s in synd):
            return symbols[:k]

        err_loc = cls.berlekamp_massey(synd)
        err_pos = cls.chien_search(err_loc, len(symbols))
        if len(err_pos) == 0:
            raise ValueError("Could not locate errors")
        # 1 Correct any errors found by the Berlekamp and Chien
        corrected = cls.forney(symbols, synd, err_pos, err_loc)
        # Verify the corrected syndromes
        check = cls.calc_syndromes(corrected, nsym)
        if not all(s == 0 for s in check):
            raise ValueError("Correction failed")
        return corrected[:k]


class RS_36_20_17(ReedSolomonScheme):
    K = 20
    nsym = 16


class RS_24_16_9(ReedSolomonScheme):
    K = 16
    nsym = 8


class RS_24_12_13(ReedSolomonScheme):
    K = 12
    nsym = 12


class ReedSolomon:
    """
    High-level interface for Reed-Solomon encoding/decoding.

    Handles:
        - Symbol-level encoding/decoding
        - Stream-to-symbol conversion via HexBits

    Designed for easy integration into systems using 6-bit symbols.
    """

    # Encoders

    @staticmethod
    def encode(symbols: Iterable[int], scheme: type[RsSchemes]) -> HexBits:
        """
        Encoder for a input symbol list.

        Encode a symbol list with a Reed-Solomon(n, k, d) code to yield 'k' hex bits.
        Where:
        - 'n' is the size of the code word
        - 'k' is the number of information bits
        - 'd' is the minimum hamming distance
        """
        symbols_ = HexBits(symbols)
        return scheme.encode(symbols_)

    @classmethod
    def encode_stream(cls, stream: BitStream, scheme: type[RsSchemes]) -> HexBits:
        """
        Encoder for a input BitStream object.

        Encode a bitstream into a symbol list of length 'n', with symbols in 6-bit formats.
        Then encode with a Reed-Solomon(n, k, d) code to yield 'k' hex bits.
        Where:
        - 'n' is the size of the code word
        - 'k' is the number of information bits
        - 'd' is the minimum hamming distance
        """
        symbols = HexBits.from_stream(stream)
        return cls.encode(symbols, scheme)

    # Decoders

    @staticmethod
    def decode(symbols: Iterable[int], scheme: type[RsSchemes]) -> HexBits:
        """
        Decoder for a input symbol list.

        Decode a symbol list with a Reed-Solomon(n, k, d) code to yield 'k' hex bits.
        Where:
        - 'n' is the size of the code word
        - 'k' is the number of information bits
        - 'd' is the minimum hamming distance
        """
        symbols_ = HexBits(symbols)
        return scheme.decode(symbols_)

    @classmethod
    def decode_stream(cls, stream: BitStream, scheme: type[RsSchemes]) -> HexBits:
        """
        Decoder for a input BitStream object.

        Encode a bitstream into a symbol list of length 'n', with symbols in 6-bit formats.
        Then decode with a Reed-Solomon(n, k, d) code to yield 'k' hex bits.
        Where:
        - 'n' is the size of the code word
        - 'k' is the number of information bits
        - 'd' is the minimum hamming distance
        """
        symbols = HexBits.from_stream(stream)
        return cls.decode(symbols, scheme)
