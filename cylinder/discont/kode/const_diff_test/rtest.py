from scipy.special import jv as bessel
from scipy.special import jn_zeros as _zeros
import math
I=100
J=100
hr=0.01
ht=0.01
d=0.1
zeros=_zeros(1, 100) #excludes 0
m=7
zm=zeros[m-2]
print("zm = %d\n" %zm)
def bes(x):
    return bessel(0, x)
def sol(r, t):
    return bes(zm*r) * math.exp(- d * zm**2 * t)

const_diff_ic = open("rconst_diff_ic.scv", "w+")
uss='\n'.join([';'.join(
                  [','.join(
                            map(str, [sol(i*hr, 0)] * 6)
                                              ) for i in range(I)]
                                                                  ) for j in range(J)])
const_diff_ic.write(uss)
const_diff_ic.close()
#print(bessel(1, zm))
each = 1
time = 0.1
const_diff_ex = open("rconst_diff_ex.csv", "w+")
n_rows = int(time/ht/each+1)
print("n_rows = %d\n" %n_rows)
times = [i*ht*each for i in range(0,n_rows)]
def sol_rs(t):
    return ','.join([str(sol(i*hr, t)) for i in range(0,I+1)])
uss='\n'.join(
        map(sol_rs, times))
const_diff_ex.write(uss)
const_diff_ex.close()
