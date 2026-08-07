from collections import UserList
from collections.abc import Iterable
from itertools import combinations
from typing import Protocol, overload

from bitarray import bitarray
from bitarray.util import int2ba
from bitstructures import BitStream


class HexBits(UserList[int]):
    SYMBOL_SIZE: int = 6
    __slots__ = ("data",)

    def __repr__(self) -> str:
        return f"{[f'⧼{i:#08b}⧽' for i in self]}"

    @classmethod
    def from_stream(cls, bitstream: BitStream) -> "HexBits":
        _bitstream = bitstream.copy()
        hex_bits = []
        while _bitstream:
            # They are then separated into 20 symbols of 6 bits each.
            hex_bits.append(int(_bitstream.read(cls.SYMBOL_SIZE)))
        return HexBits(hex_bits)

    def to_stream(self) -> BitStream:
        stream = BitStream()
        for hex_bit in self.data:
            stream.write(hex_bit, self.SYMBOL_SIZE)
        return stream


class BinaryVector(UserList[int]):
    def __repr__(self) -> str:
        string = "0b"
        for i in self.data:
            string += str(i)
        return f"⟪{string}⟫"

    @staticmethod
    def from_int(i: int, size: int) -> "BinaryVector":
        return BinaryVector(list(int2ba(i, size)))

    @staticmethod
    def from_stream(stream: BitStream, size: int) -> "BinaryVector":
        if len(stream) != size:
            raise ValueError(
                f"Cannot convert a stream of size {len(stream)} into a binary vector of size {size}"
            )
        return BinaryVector(list(stream))

    def to_stream(self) -> BitStream:
        if not all(i in [0, 1] for i in self.data):
            raise ValueError(f"All bits must be in {[0, 1]}")
        return BitStream(bitarray(self.data))

    def __xor__(self, value: Iterable[int]) -> "BinaryVector":
        return BinaryVector([x ^ y for x, y in zip(self.data, value, strict=True)])

    def weight(self) -> int:
        # symbols are binary, so we can just sum the list
        return sum(self.data)


class Matrix:
    def __init__(self, matrix: list[list[int]], *, include_transpose: bool = True) -> None:
        self._matrix = matrix
        self._transpose: Matrix | None = None
        if include_transpose:
            # Cache the transpose for quick access
            # NOTE: Don't store this as a Matrix, this will recurse
            self._transpose = Matrix(
                [[row[i] for row in self._matrix] for i in range(self.n_cols)],
                include_transpose=False,
            )

    # ----------------------------
    # Basic properties
    # ----------------------------
    @property
    def n_rows(self) -> int:
        """Returns the number of rows in the Matrix."""
        return len(self._matrix)

    @property
    def n_cols(self) -> int:
        """Returns the number of columns in the Matrix."""
        return len(self._matrix[0])

    @property
    def matrix(self) -> list[list[int]]:
        """Returns a representation of the Matrix which is a nested list of lists."""
        return self._matrix

    # ----------------------------
    # Identities & Transposes
    # ----------------------------
    @property
    def transpose(self) -> "Matrix":
        """
        Defines a transposed matrix Aᵀ calculated by reversing both rows and columns.

        Every Matrix definition will create a transpose of itself and cache the result allowing
        for quick access and calcuations on the transposed Matrix. This has large performance
        improvements over calcuating the transponse in real-time.
        """
        if self._transpose is None:
            # This matrix got defined explicitly without a transpose in cache
            self._transpose = Matrix(
                [[row[i] for row in self._matrix] for i in range(self.n_cols)],
                include_transpose=False,
            )
        return self._transpose

    @staticmethod
    def identity(n: int) -> "Matrix":
        """Creates an identity matrix of a given size n."""
        return Matrix([[1 if i == j else 0 for j in range(n)] for i in range(n)])

    # ----------------------------
    # Multiply row vector (1×n)
    # ----------------------------
    @overload
    def __rmul__(self, other: "Matrix") -> "Matrix": ...
    @overload
    def __rmul__(self, other: "BinaryVector") -> "BinaryVector": ...
    def __rmul__(self, other: "Matrix|BinaryVector") -> "Matrix|BinaryVector":
        # Vectors must be a first argument in multiplications
        if isinstance(other, Matrix):
            return self.__mul__(other)
        vec = other
        if len(vec) != self.n_rows:
            raise ValueError("Vector size mismatch")

        result = [0] * self.n_cols
        for j in range(self.n_cols):
            val = 0
            for i in range(self.n_rows):
                val ^= vec[i] & self._matrix[i][j]
            result[j] = val
        return BinaryVector(result)

    # ----------------------------
    # Matrix multiplication (GF2)
    # (A: m×n, B: n×p → m×p)
    # ----------------------------
    @overload
    def __mul__(self, other: "Matrix") -> "Matrix": ...
    @overload
    def __mul__(self, other: "BinaryVector") -> "BinaryVector": ...
    def __mul__(self, other: "Matrix|BinaryVector") -> "Matrix|BinaryVector":
        if isinstance(other, BinaryVector):
            raise ValueError("Vectors can only be multiplied from the left")
        if self.n_cols != other.n_rows:
            raise ValueError("Matrix dimension mismatch")
        result = [[0] * other.n_cols for _ in range(self.n_rows)]
        for i in range(self.n_rows):
            for j in range(other.n_cols):
                val = 0
                for k in range(self.n_cols):
                    val ^= self._matrix[i][k] & other._matrix[k][j]
                result[i][j] = val
        return Matrix(result)

    # ----------------------------
    # Conjoin matrices (A | B)
    # ----------------------------
    def __or__(self, other: "Matrix") -> "Matrix":
        """Matrix Conjoin operations, (A | B)."""
        if self.n_rows != other.n_rows:
            raise ValueError("Row mismatch for conjoin")
        return Matrix([self._matrix[i] + other.matrix[i] for i in range(self.n_rows)])

    # ----------------------------
    # Helpers
    # ----------------------------

    def shorten_row(self, i: int) -> "Matrix":
        return Matrix(self._matrix[:i])

    def shorten_col(self, i: int) -> "Matrix":
        return Matrix([row[:i] for row in self._matrix])

    def shorten_all(self, i: int) -> "Matrix":
        # remove last column AND last row
        return Matrix([row[:i] for row in self.matrix[:i]])

    def column(self, idx: int) -> "list[int]":
        return [self._matrix[i][idx] for i in range(self.n_rows)]

    def calculate_syndrome(self, vector: BinaryVector) -> BinaryVector:
        # s=r^H^T
        if not isinstance(vector, BinaryVector):
            raise ValueError("Can only calculate syndromes of vector types")
        return vector * self.transpose

    # ----------------------------
    # Debug
    # ----------------------------
    def __repr__(self) -> str:
        return "\n".join("".join(str(x) for x in row) for row in self._matrix)


class Syndromes:
    """Generic syndrome helpers for binary block codes."""

    @staticmethod
    def build_table(H: Matrix, n: int, t: int) -> dict[tuple[int, ...], BinaryVector]:
        """
        Build syndrome lookup table.

        Parameters
        ----------
        - H: Parity-check matrix.
        - n: Codeword length.
        - t: Maximum correctable weight.
        """
        table: dict[tuple[int, ...], BinaryVector] = {}

        #
        # Weight 0
        #
        zero = BinaryVector([0] * n)
        table[tuple(H.calculate_syndrome(zero))] = zero

        #
        # Weight 1..t
        #
        for weight in range(1, t + 1):
            for bits in combinations(range(n), weight):
                error = [0] * n
                for idx in bits:
                    error[idx] = 1
                error_vec = BinaryVector(error)
                syn = tuple(H.calculate_syndrome(error_vec))
                table[syn] = error_vec
        return table


class Scheme(Protocol):
    """
    Base class which most schemes will inherit from.

    All codes are linear block codes over GF(2) and use systematic
    generator matrices of the form:
        G = [ I | P ]

    Where:
        I = identity matrix
        P = parity matrix
    The corresponding parity-check matrix is:
        H = [ Pᵀ | I ]
    such that
        GHᵀ = 0
    and every valid codeword c satisfies:
        cHᵀ = 0
    """

    N: int  # Code word length
    K: int  # Information Symbols
    D: int  # Hamming Distance
    T: int  # Correcting Power

    I: Matrix  # Identity matrix
    P: Matrix  # 12×12 parity construction matrix
    G: Matrix  # Generator Matrix:    [I | P]
    H: Matrix  # Parity-Check Matrix: [Pᵀ | I]
    SYNDROME_TABLE: dict[tuple[int, ...], BinaryVector]

    @classmethod
    def encode(cls, vector: BinaryVector) -> BinaryVector:
        """
        Encode a message vector using G.

        msg: length K
        returns: length N

        Messages are encoded as:
            c = mG
        where:
            m ∈ GF(2)^k
            c ∈ GF(2)^n
        and all arithmetic is performed modulo two.
        """
        assert isinstance(vector, BinaryVector)
        if len(vector) != cls.K:
            raise ValueError("Codeword length mismatch")
        return vector * cls.G

    @classmethod
    def decode(cls, vector: BinaryVector) -> BinaryVector:
        """
        Decodes a message vector using G.

        msg: length N
        returns: length K

        Messages are decoded as:
            r = c + e
        the syndrome is computed as:
            s = rHᵀ
        Since valid codewords satisfy:
            cHᵀ = 0
        the syndrome depends only on the error vector:
            s = eHᵀ
        """
        syndrome = cls.H.calculate_syndrome(vector)
        error = cls.SYNDROME_TABLE.get(tuple(syndrome))
        if error is None:
            raise ValueError("Uncorrectable error pattern")
        return (vector ^ error)[: cls.K]


class EncoderProtocol[T: Scheme](Protocol):
    """
    High-level interface for encoding/decoding binary linear block codes.

    Handles:
        - Symbol-level encoding/decoding
        - Stream-to-symbol conversion via HexBits
    """

    @staticmethod
    def encode(vector: Iterable[int], scheme: type[T]) -> BinaryVector:
        """
        Encode parity bits into a given symbol list.

        Encode a symbol list with a (n, k, d) code to yield 'k' hex bits where:
        - 'n' is the size of the code word
        - 'k' is the number of information bits
        - 'd' is the minimum hamming distance
        """
        vector_ = BinaryVector(vector)
        return scheme.encode(vector_)

    @classmethod
    def encode_stream(cls, stream: BitStream, scheme: type[T]) -> BinaryVector:
        """
        Encode parity bits into a given BitStream.

        Encode a bitstream into a symbol list of length 'n', with symbols in 6-bit formats.
        Then encode with a (n, k, d) code to yield 'k' hex bits where:
        - 'n' is the size of the code word
        - 'k' is the number of information bits
        - 'd' is the minimum hamming distance
        """
        if len(stream) != scheme.K:
            raise ValueError(f"Input stream must be of length {scheme.K}")
        vector = BinaryVector.from_stream(stream, scheme.K)
        return cls.encode(vector, scheme)

    # Decoders

    @staticmethod
    def decode(vector: Iterable[int], scheme: type[Scheme]) -> BinaryVector:
        """
        Decode and correct errors upon a given symbol list.

        Decode a symbol list with a (n, k, d) code to yield 'k' hex bits where:
        - 'n' is the size of the code word
        - 'k' is the number of information bits
        - 'd' is the minimum hamming distance
        """
        vector_ = BinaryVector(vector)
        return scheme.decode(vector_)

    @classmethod
    def decode_stream(cls, stream: BitStream, scheme: type[Scheme]) -> BinaryVector:
        """
        Decode and correct errors upon a given BitStream.

        Decode a bitstream into a symbol list of length 'n', with symbols in 6-bit formats.
        Then encode with a (n, k, d) code to yield 'k' hex bits where:
        - 'n' is the size of the code word
        - 'k' is the number of information bits
        - 'd' is the minimum hamming distance
        """
        if len(stream) != scheme.N:
            raise ValueError(f"Input stream must be of length {scheme.N}")
        vector = BinaryVector.from_stream(stream, scheme.N)
        return cls.decode(vector, scheme)
