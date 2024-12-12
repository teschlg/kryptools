import pytest
from kryptools import SDESKeySchedule, SDESCipher, DESKeySchedule, DESCipher, list2int


def test_SDES():
    # Test the key schedule
    key = 34
    test_values = [[0, 0, 1, 0, 0, 0, 0, 0], [0, 0, 0, 1, 0, 0, 1, 0]]

    ks = SDESKeySchedule(key)
    assert ks.gen_keys() == test_values

    # Test the encryption
    x = b'\x01'
    y = b'\xb9'

    sdes = SDESCipher(key)
    sdes(x).encrypt()
    assert bytes(sdes) == y
    assert sdes.decrypt() == x


def test_DES():
    # Test the key schedule
    key = 0x0123456789ABCDEF.to_bytes(8, 'big')
    test_values = [0x0B02679B49A5, 0x69A659256A26, 0x45D48AB428D2, 0x7289D2A58257, 0x3CE80317A6C2, 0x23251E3C8545, 0x6C04950AE4C6, 0x5788386CE581,
                   0xC0C9E926B839, 0x91E307631D72, 0x211F830D893A, 0x7130E5455C54, 0x91C4D04980FC, 0x5443B681DC8D, 0xB691050A16B5, 0xCA3D03B87032]

    ks = DESKeySchedule(key)
    keys = ks.gen_keys()
    for k1, k2 in zip(keys, test_values):
        assert list2int(k1) == k2

    # Test the encryption (Handbook of Applied Cryptography)
    key = 0x0123456789ABCDEF.to_bytes(8, 'big')
    x = 0x4E6F772069732074.to_bytes(8, 'big')
    y = 0x3FA40E8A984D4815.to_bytes(8, 'big')

    des = DESCipher(key)
    des(x).encrypt()
    assert bytes(des) == y
    assert des.decrypt() == x
