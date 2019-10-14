########################################################################
# This file contains functions required for primality test using AKS   #
# Needed to implement Karp-Rabin string matching algorithm             #
########################################################################

import math

def binpow(a,b):
  """
  Computes a**b efficiently for a and b integers.
  Arguments:
    a : integer base
    b : integer power
  Returns:
    result : integer, a ** b.
  """
  result = 1
  while b > 0:
    if b & 1: # if lowest bit is 1
      result *= a
    a *= a # square a
    b >>= 1 # drop the lowest bit
  return result

def prime_factors(n):
  """
  Finds prime factors of integer n. Will be used to compute Euler's totient function.
  Argument:
    n : integer
  Returns:
    factors : a list of prime numbers that divides n
  """
  factors = []
  if n %2 == 0: # take care of even-ness
    count = 0    
    while n%2 == 0:
      count += 1
      n = n // 2
    factors.append([2, count])

  m = 3
  while m < n: # do all the odds: m increasing and n decreasing with each factor
    if n%m == 0:
      count = 0
      while n%m==0:
        count += 1
        n = n // m
      factors.append([m, count])
    m += 2

  if n!=1: # no factors found, i.e, n is itself prime
    factors.append([n,1])
  return factors


def euler_totient(n):
  """
  Compute Euler's totient function using Euler theorem: phi(n) = n Product(p|n) (1-1/p)
  """
  ps = [x[0] for x in prime_factors(n)]
  result = n*1.0
  for p in ps:
    result *= (p-1)/p
  return int(result)


############ More functions needed below

if __name__ == "__main__":
  # a = 3
  # b = 13
  # c = binpow(a,b)
  # print( "{} ** {} = {}".format(a, b, c))

  #n = 234567892309567825
  #n = 36
  # n = 2345678923*7*7524
  # print( prime_factors(n) )
  nn = [2, 3, 4, 5, 6, 7, 8, 9, 10, 23456780]
  for n in nn:
    print("phi({}) = {}".format(n, euler_totient(n)) )