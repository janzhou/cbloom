from __future__ import absolute_import
import math
import hashlib
from struct import unpack, pack, calcsize

import sys
try:
    import StringIO
    import cStringIO
except ImportError:
    from io import BytesIO

__version__ = '0.1'
__author__  = "Jian Zhou <i@janzhou.org>,\
              "

running_python_3 = sys.version_info[0] == 3

def range_fn(*args):
    if running_python_3:
        return range(*args)
    else:
        return xrange(*args)

def is_string_io(instance):
    if running_python_3:
       return isinstance(instance, BytesIO)
    else:
        return isinstance(instance, (StringIO.StringIO,
                                     cStringIO.InputType,
                                     cStringIO.OutputType))

def make_hashfuncs(num_slices, num_bits):
    if num_bits >= (1 << 31):
        fmt_code, chunk_size = 'Q', 8
    elif num_bits >= (1 << 15):
        fmt_code, chunk_size = 'I', 4
    else:
        fmt_code, chunk_size = 'H', 2
    total_hash_bits = 8 * num_slices * chunk_size
    if total_hash_bits > 384:
        hashfn = hashlib.sha512
    elif total_hash_bits > 256:
        hashfn = hashlib.sha384
    elif total_hash_bits > 160:
        hashfn = hashlib.sha256
    elif total_hash_bits > 128:
        hashfn = hashlib.sha1
    else:
        hashfn = hashlib.md5
    fmt = fmt_code * (hashfn().digest_size // chunk_size)
    num_salts, extra = divmod(num_slices, len(fmt))
    if extra:
        num_salts += 1
    salts = tuple(hashfn(hashfn(pack('I', i)).digest()) for i in range_fn(num_salts))
    def _make_hashfuncs(key):
        if running_python_3:
            if isinstance(key, str):
                key = key.encode('utf-8')
            else:
                key = str(key).encode('utf-8')
        else:
            if isinstance(key, unicode):
                key = key.encode('utf-8')
            else:
                key = str(key)
        i = 0
        for salt in salts:
            h = salt.copy()
            h.update(key)
            for uint in unpack(fmt, h.digest()):
                yield uint % num_bits
                i += 1
                if i >= num_slices:
                    return

    return _make_hashfuncs


class CBloom(object):
    FILE_FMT = b'<dQQQQ'

    def __init__(self, capacity, error_rate=0.001):
        """Implements a space-efficient probabilistic data structure

        capacity
            this BloomFilter must be able to store at least *capacity* elements
            while maintaining no more than *error_rate* chance of false
            positives
        error_rate
            the error_rate of the filter returning false positives. This
            determines the filters capacity. Inserting more than capacity
            elements greatly increases the chance of false positives.

        >>> b = BloomFilter(capacity=100000, error_rate=0.001)
        >>> b.add("test")
        False
        >>> "test" in b
        True

        """
        if not (0 < error_rate < 1):
            raise ValueError("Error_Rate must be between 0 and 1.")
        if not capacity > 0:
            raise ValueError("Capacity must be > 0")
        # given M = num_bits, k = num_slices, P = error_rate, n = capacity
        #       k = log2(1/P)
        # solving for m = bits_per_slice
        # n ~= M * ((ln(2) ** 2) / abs(ln(P)))
        # n ~= (k * m) * ((ln(2) ** 2) / abs(ln(P)))
        # m ~= n * abs(ln(P)) / (k * (ln(2) ** 2))
        num_slices = int(math.ceil(math.log(1.0 / error_rate, 2)))
        counters_per_slice = int(math.ceil(
            (capacity * abs(math.log(error_rate))) /
            (num_slices * (math.log(2) ** 2))))
        self._setup(error_rate, num_slices, counters_per_slice, capacity, 0)
        self.counters = [0] * self.num_counters

    def _setup(self, error_rate, num_slices, counters_per_slice, capacity, count):
        self.error_rate = error_rate
        self.num_slices = num_slices
        self.counters_per_slice = counters_per_slice
        self.capacity = capacity
        self.num_counters = num_slices * counters_per_slice
        self.capacity_count = count
        self.make_hashes = make_hashfuncs(self.num_slices, self.counters_per_slice)

    def count(self, key):
        counters_per_slice = self.counters_per_slice
        counters = self.counters
        hashes = self.make_hashes(key)
        offset = 0
        count = 0
        for k in hashes:
            if count < counters[offset + k]:
                count = counters[offset + k]
            offset += counters_per_slice
        offset = 0
        for k in hashes:
            if count > counters[offset + k]:
                count = counters[offset + k]
            offset += counters_per_slice
        return count

    def __len__(self):
        """Return the number of keys stored by this bloom filter."""
        return self.count

    def add(self, key):
        counters_per_slice = self.counters_per_slice
        hashes = self.make_hashes(key)
        if self.capacity_count > self.capacity:
            raise IndexError("BloomFilter is at capacity")
        hit = count = 0
        offset = 0
        for k in hashes:
            if self.counters[offset + k] > 0:
                hit += 1
            count += 1
            self.counters[offset + k] += 1
            offset += counters_per_slice
        if hit < count:
            self.capacity_count += 1

    def rm(self, key):
        counters_per_slice = self.counters_per_slice
        hashes = self.make_hashes(key)
        hit = count = 0
        offset = 0
        for k in hashes:
            if self.counters[offset + k] > 0:
                hit += 1
            count += 1
            offset += counters_per_slice
        if hit == count:
            offset = 0
            hashes = self.make_hashes(key)
            for k in hashes:
                self.counters[offset + k] -= 1
                offset += counters_per_slice
        for k in hashes:
            if self.counters[offset + k] > 0:
                hit += 1
            count += 1
            offset += counters_per_slice
        if hit < count:
            self.capacity_count -= 1

