#!/usr/bin/env python
# @file        libeuler.py
# @author      Michael Foukarakis
# @version     0.4
# @date        Created:     Tue Oct 11, 2011 08:51 GTB Daylight Time
#              Last Update: Τετ Ιουν 01, 2016 09:14 GTB Daylight Time
#------------------------------------------------------------------------
# Description: Project Euler helper library
#------------------------------------------------------------------------
# History:     Unimportant.
# TODO:        Nothing yet.
#------------------------------------------------------------------------
# -*- coding: utf-8 -*-
#------------------------------------------------------------------------
from itertools import tee, filterfalse
from functools import reduce
from  operator import mul
from    random import randrange
from      math import ceil

def factorial(n):
    return reduce(lambda x, y : x*y, range(1, n+1), 1)

def is_permutation(a, b):
    return sorted(a) == sorted(b)

def is_palindromic(n):
    return str(n) == str(n)[::-1]

def is_pandigital(n, s=9):
    n = str(n)
    return len(n) == s and not '1234567890'[:s].strip(n)

def is_prime(n):
    """Deterministic primality test based on the P3 prime candidate generator.
    """
    if n in [2, 3, 5]:
        return True
    if n == 1 or n & 1 == 0:
        return False
    if n > 5 and (n % 6 not in [1, 5] or n % 5 == 0):
        return False
    for c in range(7, isqrt(n), 2):
        p1, k, p2 = 5 * c, 6 * c, 7 * c
        if (n - p1) % k == 0 or (n - p2) % k == 0:
            return False
    else:
        return True


def miller_rabin_pass(a, s, d, n):
    """Miller-Rabin single pass.
    """
    a_to_power = pow(a, d, n)
    if a_to_power == 1:
        return True
    for i in range(s-1):
        if a_to_power == n - 1:
            return True
        a_to_power = (a_to_power * a_to_power) % n
    return a_to_power == n - 1

def miller_rabin(n):
    """Miller-Rabin primality test. Returns True if N is prime.
    Reference: http://en.literateprograms.org/Miller-Rabin_primality_test_(Python)?action=history&offset=20101013093632
    """
    d = n - 1
    s = 0
    while d % 2 == 0:
        d >>= 1
        s += 1
    for repeat in range(20):
        a = 0
        while a == 0:
            a = randrange(n)
        if not miller_rabin_pass(a, s, d, n):
            return False
    return True

def lucas_lehmer_test(p):
    """Lucas-Lehmer primality test for Mersenne numbers. Returns True if 2**p-1 is prime.
    """
    s = 4
    M = 2**p - 1
    for _ in range(p - 2):
        s = ((s * s) - 2) % M
    return s == 0


def trial_division(n, bound=None):
    """Tests if N can be divided exactly by any integer in [2, max(N, BOUND)). Returns the divisor, or N.
    """
    if n == 1:
        return 1
    for p in [2, 3, 5]:
        if n%p == 0:
            return p
    if bound == None:
        bound = n
    dif = [6, 4, 2, 4, 2, 4, 6, 2]
    i, m = 1, 7
    while m <= bound and m*m <= n:
        if n % m == 0:
            return m
        m += dif[i%8]
        i += 1
    return n

def factor(n):
    """Returns a dictionary with n = product over p ^ d[p].  Seems somewhat faster than
    prime_factors().
    """
    if n in [-1, 0, 1]: return {}
    n = abs(n)
    F = {}
    while n != 1:
        p = trial_division(n)
        e = 1
        n //= p
        while n % p == 0:
            e += 1
            n //= p
        F[p] = e
    return F

def p6():
    """Generates prime number candidates (6*n + 1, 6*n + 5). Used in prime factorization
    routine prime_factors().
    """
    yield 2
    yield 3
    i = 5
    while True:
        yield i
        i += 2
        yield i
        i += 4

def prime_factors(n):
    """Returns a dictionary with n = product p ^ d[p].
    """
    d = {}
    primes = p6()
    for p in primes:
        while n % p == 0:
            n /= p
            d[p] = d.setdefault(p, 0) + 1
        if n == 1:
            return d

def number_of_divisors(n, square=False):
    """Returns the number of divisors of N, computing it as the product of its prime
    factors.
    When SQUARE is True, it returns the number of divisors of N^2, knowing that:
    N^2 = (p1^a1 * p2^a2 * ... * pk^ak)^2
    """
    d = prime_factors(n)
    if square:
        powers_plus = map(lambda x: 2*x + 1, d.values())
    else:
        powers_plus = map(lambda x:   x + 1, d.values())
    return reduce(mul, powers_plus, 1)

def gcd(a, b):
    """Returns the GCD of A and B.
    """
    a, b = abs(a), abs(b)
    if a == 0:
        return b
    if b == 0:
        return a
    while b != 0:
        (a, b) = (b, a%b)
    return a

def permutation(n, s):
    """Returns a permutation of sequence S.
    """
    if len(s) == 1:
        return s
    q, r = divmod(n, factorial(len(s)-1))
    return s[q] + permutation(r, s[:q] + s[q+1:])

def partition(seq, key=bool):
    """Returns two iterators that filter elements from sequence SEQ. The first one returns
    those elements for which KEY is True, and the second those for which KEY is False.
    """
    s1, s2 = tee(seq)
    return filter(key, s1), filterfalse(key, s2)

def binomial(n, k):
    nt = 1
    for t in range(min(k, n-k)):
        nt = nt * (n - t) // (t + 1)
    return nt

# My awesome prime sieve.
def prime_sieve(n):
    """Input n >= 6,
    Returns a list of primes p, where 2 <= p < n
    """
    import numpy as np
    sieve = np.ones(n/3 + (n % 6 == 2), dtype = np.bool)
    sieve[0] = False
    for i in range(int(n**0.5) // 3 + 1):
        if sieve[i]:
            k = 3 * i + 1|1
            sieve[      ((k*k)//3)      ::2*k] = False
            sieve[(k*k+4*k-2*k*(i&1))//3::2*k] = False
    return np.r_[2, 3, ((3 * np.nonzero(sieve)[0] + 1)|1)].tolist()

def triangle_maximal_sum(t):
    """Returns the maximal path from the root of a tree to a leaf.  t is a list of lists,
    representing the tree top-down.
    """
    for row in range(len(t)-1, 0, -1):
        for col in range(0, row):
            dt[row-1][col] += max(t[row][col], t[row][col+1])
    return t[0][0]

def s(n0, primelist):
    """Returns the aliquot sum of n - sum of its proper divisors
    """
    n, i, p, res = n0, 0, primelist[0], 1
    while p * p <= n and n > 1 and i < len(primelist):
        p = primelist[i]
        i += 1
        if n % p == 0:
            j = p * p
            n = n / p
            while n % p == 0:
                j = j * p
                n = n / p
            res = res * (j - 1) / (p - 1)
    if n > 1:
        res *= n + 1
    return res - n0

def polygonal(dim):
    """Generates a series of polygonal numbers for n = [3, 8].
    """
    if dim < 3 or dim > 8:
        return NotImplemented
    n = 0
    func = {
        3 : lambda x: x * (x + 1) // 2,
        4 : lambda x: x ** 2,
        5 : lambda x: x * (3*x - 1) // 2,
        6 : lambda x: x * (2*x - 1),
        7 : lambda x: x * (5*x - 3) // 2,
        8 : lambda x: x * (3*n - 2)
    }[dim]
    while True:
        n += 1
        yield func(n)

def cf(n):
    """Yields the continued fraction expansion of the square root of integer n, or 0 if n
    is a perfect square.
    """
    from math import sqrt
    a0 = int(sqrt(n))      # floor of sqrt(r)
    if a0*a0 == n: yield 0 # perfect square

    d, m, a = 1, 0, a0

    while True:
        yield a

        m = d * a - m
        d = (n - m * m) / d
        a = (a0 + m) / d

        if a == 2 * a0:
            break

def phi(n):
    """Computes the Euler's totient function φ(n) - number of positive numbers less than
    or equal to n which are relatively prime to n.
    """
    from functools import reduce
    from operator  import mul
    return n * reduce(mul, [(1 - 1 / p) for p in prime_factors(n)], 1)

if __name__ == '__main__':
    import doctest
    doctest.testmod()
