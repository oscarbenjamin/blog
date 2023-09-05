from operator import mul
from itertools import starmap
from collections import deque
import random
from time import time
import sympy
import flint

x = sympy.symbols('x')

from gmpy2 import mpz, mpq
import matplotlib.pyplot as plt
import numpy as np

def rand_ZZ_poly(degree):
    return [random.randint(1, 10) for _ in range(degree+1)]

def time_f_ns(func, items, _time=time):
    tstart = _time()
    for item in items:
        func(item)
    tdelta = _time()  - tstart
    return tdelta * 1e6

def time_sympy_flint(degree, nrepeats):
    c1 = rand_ZZ_poly(degree)
    c2 = rand_ZZ_poly(degree)
    p3_sympy = [sympy.Poly(c1, x, domain='ZZ') * sympy.Poly(c2, x, domain='ZZ') for _ in range(nrepeats)]
    p3_flint = [flint.fmpz_poly(c1[::-1]) * flint.fmpz_poly(c2[::-1]) for _ in range(nrepeats)]
    if degree < 150:
        t_sympy = time_f_ns(lambda p: p.factor_list(), p3_sympy) / nrepeats
    else:
        t_sympy = None
    t_flint = time_f_ns(lambda p: p.factor(), p3_flint) / nrepeats
    return t_sympy, t_flint

degrees = [2, 3, 4, 5, 7, 10, 20, 30, 40, 50, 70, 80, 90, 100, 110, 120, 130, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
times_sympy = []
times_flint = []

t_sympy, t_flint = time_sympy_flint(1, 1)

for degree in degrees:
    nrepeats = 10000 // degree**2 + 1
    print(f'degree {degree}: nrepeat {nrepeats}')
    t_sympy, t_flint = time_sympy_flint(degree, nrepeats)
    if t_sympy is not None:
        times_sympy.append(t_sympy)
    times_flint.append(t_flint)

fig = plt.figure()

ax = fig.add_subplot(1, 1, 1)
ax.set_title('ZZ[x] factorisation time vs degree')

n = len(times_sympy)

ax.plot(degrees[:len(times_sympy)], times_sympy, marker='+', label='SymPy')
ax.plot(degrees, times_flint, marker='.', label='FLINT')
ax.plot(degrees[:n], [200*d**2 for d in degrees[:n]], label=r'$200d^2$')
ax.plot(degrees, [3*d**2 for d in degrees], label=r'$3d^2$')
ax.set_xscale('log')
ax.set_yscale('log')
#ax.set_xticks(xticks)
#ax.set_xticklabels([str(x) for x in xticks], rotation=90)
#ax.set_yticks(yticks)
#ax.set_yticklabels([str(y) for y in yticks])
ax.set_xlabel('degree')
ax.set_ylabel('time (microseconds)')

ax.legend()
fig.savefig('time_factor.png')
fig.savefig('time_factor.svg')
fig.savefig('time_factor.pdf')
plt.show()
