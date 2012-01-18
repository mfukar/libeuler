#!/usr/bin/env python
# @file        libeuler.py
# @author      Michael Foukarakis
# @version     0.1
# @date        Created:     Tue Oct 11, 2011 08:51 GTB Daylight Time
#              Last Update: Wed Jan 18, 2012 11:56 GTB Standard Time
#------------------------------------------------------------------------
# Description: Project Euler helper library
#------------------------------------------------------------------------
# History:     Unimportant.
# TODO:        Nothing yet.
#------------------------------------------------------------------------
# -*- coding: utf-8 -*-
#------------------------------------------------------------------------
import random
from functools import reduce
from  operator import mul
from      math import sqrt, ceil

def factorial(n):
    return reduce(lambda x, y : x*y, range(1, n+1), 1)

def is_permutation(a,b):
    return sorted(str(a)) == sorted(str(b))

def is_palindromic(n):
    return str(n) == str(n)[::-1]

def is_pandigital(n, s = 9):
    n = str(n)
    return len(n) == s and not '1234567890'[:s].strip(n)

def is_prime(n):
    if n == 2 or n == 3:
        return True
    if n < 2 or n % 2 == 0:
        return False
    if n < 9:
        return True
    if n % 3 == 0:
        return False
    r = int(sqrt(n))
    f = 5
    while f <= r:
        if n%f == 0: return False
        if n%(f+2) == 0: return False
        f +=6
    return True


# http://en.literateprograms.org/Miller-Rabin_primality_test_(Python)?action=history&offset=20101013093632
def miller_rabin_pass(a, s, d, n):
    a_to_power = pow(a, d, n)
    if a_to_power == 1:
        return True
    for i in range(s-1):
        if a_to_power == n - 1:
            return True
        a_to_power = (a_to_power * a_to_power) % n
    return a_to_power == n - 1


def miller_rabin(n):
    d = n - 1
    s = 0
    while d % 2 == 0:
        d >>= 1
        s += 1
    for repeat in range(20):
        a = 0
        while a == 0:
            a = random.randrange(n)
        if not miller_rabin_pass(a, s, d, n):
            return False
    return True


def trial_division(n, bound=None):
    if n == 1: return 1
    for p in [2, 3, 5]:
        if n%p == 0:
            return p
    if bound == None:
        bound = n
    dif = [6, 4, 2, 4, 2, 4, 6, 2]
    m = 7; i = 1
    while m <= bound and m*m <= n:
        if n%m == 0:
            return m
        m += dif[i%8]
        i += 1
    return n


def factor(n):
    if n in [-1, 0, 1]: return []
    if n < 0: n = -n
    F = []
    while n != 1:
        p = trial_division(n)
        e = 1
        n /= p
        while n%p == 0:
            e += 1; n /= p
        F.append((p,e))
    F.sort()
    return F


# Generator for a superset of primes
# Used for prime factorization, below:
def primes_plus():
    yield 2
    yield 3
    i = 5
    while True:
        yield i
        if i % 6 == 1:
            i += 2
        i += 2

# Returns a dictionary with n = product p ^ d[p]
def prime_factors(n):
    d = {}
    primes = primes_plus()
    for p in primes:
        while n % p == 0:
            n /= p
            d[p] = d.setdefault(p, 0) + 1
        if n == 1:
            return d

def number_of_divisors(n):
    d = prime_factors(n)
    powers_plus = map(lambda x: x + 1, d.values())
    return reduce(mul, powers_plus, 1)

def gcd(a, b):
    a, b = abs(a), abs(b)
    if a == 0:
        return b
    if b == 0:
        return a
    while b != 0:
        (a, b) = (b, a%b)
    return a

def permutation(n, s):
    if len(s)==1: return s
    q, r = divmod(n, factorial(len(s)-1))
    return s[q] + permutation(r, s[:q] + s[q+1:])

def binomial(n, k):
    nt = 1
    for t in range(min(k, n-k)):
        nt = nt*(n-t)//(t+1)
    return nt

# My awesome prime sieve.
import numpy as np
def prime_sieve(n):
    """ Input n >= 6, Returns a list of primes, 2 <= p < n """
    sieve = np.ones(n/3 + (n%6==2), dtype=np.bool)
    sieve[0] = False
    for i in range(int(n**0.5)//3+1):
        if sieve[i]:
            k=3*i+1|1
            sieve[      ((k*k)//3)      ::2*k] = False
            sieve[(k*k+4*k-2*k*(i&1))//3::2*k] = False
    return np.r_[2, 3, ((3 * np.nonzero(sieve)[0]+1)|1)].tolist()

# Find the maximal path from the root of a tree to a leaf
# t is a list of lists, representing the tree top-down
def triangle_maximal_sum(t):
    for row in range(len(t)-1, 0, -1):
        for col in range(0, row):
            dt[row-1][col] += max(t[row][col], t[row][col+1])
    return t[0][0]

# Returns the aliquot sum of n - sum of its proper divisors
def s(n):
    sum = 1
    t = sqrt(n)
    # only proper divisors, please
    for i in range(2, int(t)+1):
        if n % i == 0: sum += i + n / i
    if t == int(t):
        sum -= t
    return sum
