from collections.abc import Iterable
from random import randint

from pytest import mark

from binarycorrectioncodes.common import HexBits
from binarycorrectioncodes.reed_solomon import (
    GF6,
    RS_24_12_13,
    RS_24_16_9,
    RS_36_20_17,
    ReedSolomon,
    RsSchemes,
)


def test_gf6_field() -> None:
    # 1 Every operation must stay inside GF(2⁶)
    assert all(GF6.EXP_TABLE[GF6.LOG_TABLE[a]] == a for a in range(1, 64))
    # 2 Log/Exp consistency (critical for table correctness)
    assert all(GF6.LOG_TABLE[GF6.EXP_TABLE[i]] == i for i in range(63))
    # 3 Multiplication correctness via logs
    assert all(GF6.gf_pow(2, i) == GF6.EXP_TABLE[i] for i in range(63))
    # 4 Power function correctness
    assert all(GF6.gf_mul(a, 1) == a for a in range(64))
    # 5 Multiplicative inverses
    assert all(GF6.gf_pow(a, 63) == 1 for a in range(1, 64))
    # Addition
    assert all(GF6.gf_add(a, a) == 0 for a in range(64))
    # Inverse
    for a in range(1, 64):
        assert GF6.gf_mul(a, GF6.gf_inv(a)) == 1
    # Cancellation
    for a in range(1, 64):
        for b in range(64):
            assert GF6.gf_mul(GF6.gf_mul(b, a), GF6.gf_inv(a)) == b
    # Power relation
    for a in range(1, 64):
        assert GF6.gf_inv(a) == GF6.gf_pow(a, 62)
    # core correctness
    for a in range(64):
        for b in range(1, 64):
            assert GF6.gf_mul(GF6.gf_div(a, b), b) == a
    # identity
    assert all(GF6.gf_div(a, 1) == a for a in range(64))
    # inverse relation
    assert all(GF6.gf_div(1, a) == GF6.gf_inv(a) for a in range(1, 64))


@mark.parametrize("scheme", [RS_36_20_17, RS_24_16_9, RS_24_12_13])
def fuzz_test(scheme: type[RsSchemes], trials: int = 100) -> None:
    for _ in range(trials):
        msg = [randint(0, GF6.MAX_EXP) for _ in range(scheme.K)]
        enc = ReedSolomon.encode(msg, scheme)

        test = enc.copy()
        for _ in range(scheme.nsym // 2):
            idx = randint(0, len(enc) - 1)
            test[idx] ^= randint(1, 63)

        dec = ReedSolomon.decode(test, scheme)
        assert dec == msg


@mark.parametrize(
    ("message", "scheme"),
    [
        ([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 1, 0, 1], RS_36_20_17),
        ([0, 0, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 1, 0, 1], RS_24_16_9),
        ([0, 0, 0, 0, 32, 0, 0, 0, 0, 1, 0, 1], RS_24_12_13),
    ],
)
def test(message: Iterable[int], scheme: type[RsSchemes]) -> None:
    print(f"{scheme.__name__}")
    msg = HexBits(message)  # 1010 0101 0101
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
    [
        ([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 1, 0, 1], RS_36_20_17),
        ([0, 0, 0, 0, 0, 0, 0, 0, 32, 0, 0, 0, 0, 1, 0, 1], RS_24_16_9),
        ([0, 0, 0, 0, 32, 0, 0, 0, 0, 1, 0, 1], RS_24_12_13),
    ],
)
def test_no_error(message: Iterable[int], scheme: type[RsSchemes]) -> None:
    print(f"{scheme.__name__}")
    msg = HexBits(message)  # 1010 0101 0101
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
