from pytest import mark

from binarycorrectioncodes.common import BinaryVector
from binarycorrectioncodes.golay import Golay_18_6_8, Golay_23_12_7, Golay_24_12_8, GolaySchemes


@mark.parametrize(
    ("message", "scheme"),
    [(0xA55, Golay_24_12_8), (0xA55, Golay_23_12_7), (0b110110, Golay_18_6_8)],
)
def test(message: int, scheme: type[GolaySchemes]) -> None:
    print(f"{scheme.__name__}")
    msg = BinaryVector.from_int(message, scheme.K)  # 1010 0101 0101
    # Encode
    print(f"Original codeword ({len(msg)}):\t", msg)
    cw = scheme.encode(msg)
    print(f"Encoded codeword ({len(cw)}):\t", cw)
    # Introduce errors (flip up to 3 bits)
    corrupted = cw.copy()
    corrupted[0] ^= 1
    corrupted[1] ^= 1
    # corrupted[3]^=1
    # corrupted[5]^=1
    print(f"Corrupted codeword ({len(corrupted)}):", corrupted)
    # Decode
    decoded = scheme.decode(corrupted)
    print(f"Decoded codeword ({len(decoded)}):\t", decoded)
    # Extract original message
    assert decoded == msg
    print("Decode successful")
    print()


@mark.parametrize(
    ("message", "scheme"),
    [(0xA55, Golay_24_12_8), (0xA55, Golay_23_12_7), (0b110110, Golay_18_6_8)],
)
def test_no_error(message: int, scheme: type[GolaySchemes]) -> None:
    print(f"{scheme.__name__}")
    msg = BinaryVector.from_int(message, scheme.K)  # 1010 0101 0101
    # Encode
    print(f"Original codeword ({len(msg)}):\t", msg)
    cw = scheme.encode(msg)
    print(f"Encoded codeword ({len(cw)}):\t", cw)
    # Decode
    decoded = scheme.decode(cw)
    print(f"Decoded codeword ({len(decoded)}):\t", decoded)
    # Extract original message
    assert decoded == msg
    print("Decode successful")
    print()
