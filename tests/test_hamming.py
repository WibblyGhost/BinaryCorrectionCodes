from pytest import mark

from binarycorrectioncodes.common import BinaryVector
from binarycorrectioncodes.hamming import Hamming_10_6_3, Hamming_15_11_3, HammingSchemes


@mark.parametrize(("message", "scheme"), [(0b10101, Hamming_15_11_3), (0b10101, Hamming_10_6_3)])
def test(message: int, scheme: type[HammingSchemes]) -> None:
    print(f"{scheme.__name__}")
    msg = BinaryVector.from_int(message, scheme.K)  # 1010 0101 0101
    # Encode
    print(f"Original codeword ({len(msg)}):\t", msg)
    cw = scheme.encode(msg)
    print(f"Encoded codeword ({len(cw)}):\t", cw)
    # Introduce errors (flip up to 3 bits)
    corrupted = cw.copy()
    corrupted[0] ^= 1
    # corrupted[1]^=1
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


@mark.parametrize(("message", "scheme"), [(0b10101, Hamming_15_11_3), (0b10101, Hamming_10_6_3)])
def test_no_error(message: int, scheme: type[HammingSchemes]) -> None:
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
    print("Decode Successful")
    print()
