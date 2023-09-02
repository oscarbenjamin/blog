from operator import mul
from itertools import starmap
from collections import deque
import random
from time import time
from sympy.external.pythonmpq import PythonMPQ

from gmpy2 import mpz, mpq
import matplotlib.pyplot as plt

def rand_ZZ(bits):
    return random.randint(2**(bits-1), 2**bits-1)

def rand_QQ(bits):
    return pythonmpq(rand_ZZ(bits), rand_ZZ(bits))

def iconsume(iterator):
    deque(iterator, 0)

def time_f_ns(func, args, _time=time, _iconsume=iconsume):
    results = starmap(func, args)
    tstart = _time()
    _iconsume(results)
    tdelta = _time()  - tstart
    return tdelta * 1e9

def time_cpython_gmpy2(size, nrepeats):
    vals_cpython = [(rand_ZZ(size), rand_ZZ(size)) for _ in range(nrepeats)]
    vals_gmpy2 = [(mpz(z1), mpz(z2)) for z1, z2 in vals_cpython]
    t_cpython = time_f_ns(mul, vals_cpython)# / nrepeats
    t_gmpy2 = time_f_ns(mul, vals_gmpy2)# / nrepeats
    return t_cpython, t_gmpy2

bitsizes = list(range(1, 32))
nmax = 10 # go up to 2**(nmax+5) bits
for n in range(1, nmax+1):
    bitsizes += list(range(32*2**(n-1), 32*2**n, 2**n))

xticks = [2**n for n in range(bitsizes[-1].bit_length() + 1)]
yticks = [10, 100, 1000]

times_cpython = []
times_gmpy2 = []

t_cpython, t_gmpy2 = time_cpython_gmpy2(1, 10000)

for size in bitsizes:
    nrepeats = 1000 if size < 512 else int(1000/(25*size**0.5)) + 1
    print(f'{size} bits, {nrepeats} repeats')
    t_cpython, t_gmpy2 = time_cpython_gmpy2(size, nrepeats)
    times_cpython.append(t_cpython)
    times_gmpy2.append(t_gmpy2)

fig = plt.figure()

ax = fig.add_subplot(1, 1, 1)
ax.set_title('ZZ multiplication time vs bit size')

ax.plot(bitsizes, times_cpython, marker='+', label='CPython int')
ax.plot(bitsizes, times_gmpy2, marker='.', label='gmpy2 mpz')
ax.plot(bitsizes, [25000*b**0.5 for b in bitsizes])
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xticks(xticks)
ax.set_xticklabels([str(x) for x in xticks], rotation=90)
ax.set_yticks(yticks)
ax.set_yticklabels([str(y) for y in yticks])
ax.set_xlabel('bit size')
ax.set_ylabel('time (nanoseconds)')

ax.legend()
plt.show()
