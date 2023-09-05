from operator import mul
from itertools import starmap
from collections import deque
import random
from time import time
from sympy.external.pythonmpq import PythonMPQ

from gmpy2 import mpz, mpq
import matplotlib.pyplot as plt
import numpy as np

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
    if size < 2*8*1000**2:
        t_cpython = time_f_ns(mul, vals_cpython) / nrepeats
    else:
        t_cpython = None
    t_gmpy2 = time_f_ns(mul, vals_gmpy2) / nrepeats
    return t_cpython, t_gmpy2

bitsizes = list(range(1, 32))
nmax = 25 # go up to 2**(nmax+5) bits
for n in range(1, nmax+1):
    bitsizes += list(range(32*2**(n-1), 32*2**n, 2**n))

xticks = [2**n for n in range(bitsizes[-1].bit_length() + 1)]
yticks = [10, 100, 1000, 10000, 100000, 1000000]

times_cpython = []
times_gmpy2 = []

t_cpython, t_gmpy2 = time_cpython_gmpy2(1, 10000)

N = 100000

for size in bitsizes:
    nrepeats = N if size < 512 else 1 + int(N*(512/size)**1.5)
    print(f'{size} bits, {nrepeats} repeats')
    t_cpython, t_gmpy2 = time_cpython_gmpy2(size, nrepeats)
    if t_cpython is not None:
        times_cpython.append(t_cpython)
    times_gmpy2.append(t_gmpy2)

fig = plt.figure()

ax = fig.add_subplot(1, 1, 1)
ax.set_title('ZZ multiplication time vs bit size')

ax.plot(bitsizes[:len(times_cpython)], times_cpython, marker='+', label='CPython int')
ax.plot(bitsizes, times_gmpy2, marker='.', label='gmpy2 mpz')
bitsizes2 = [b for b in bitsizes if b >= 1024]
ax.plot(bitsizes2, [0.1*b**1.58 for b in bitsizes2], '--', linewidth=3, label=r'$n^{1.58}$')
bitsizes3 = [b for b in bitsizes if b >= 100_000]
ax.plot(bitsizes3, [0.75*b**1.2 for b in bitsizes3], '--', linewidth=3, label=r'$n^{1.2}$')
ax.set_xscale('log')
ax.set_yscale('log')
#ax.set_xticks(xticks)
#ax.set_xticklabels([str(x) for x in xticks], rotation=90)
#ax.set_yticks(yticks)
#ax.set_yticklabels([str(y) for y in yticks])
ax.set_xlabel('bit size')
ax.set_ylabel('time (nanoseconds)')

ax.legend()
fig.savefig('time_ZZ.png')
fig.savefig('time_ZZ.svg')
fig.savefig('time_ZZ.pdf')
plt.show()
