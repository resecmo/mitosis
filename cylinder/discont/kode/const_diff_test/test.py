from scipy.special import jv as bessel
from scipy.special import jn_zeros as _zeros
import math
I=100
J=100
hr=0.01
ht=0.1
zeros=_zeros(1, 4) #excludes 0
m=3
zm=zeros[m-2]
def bes(x):
    return bessel(0, x)
def sol(r, t):
    return bes(zm*r) * math.exp(- zm**2 * t)

const_diff_ic = open("const_diff_ic.scv", "w+")
uss='\n'.join([';'.join(
                  [','.join(
                            map(str, [sol(i*hr, 0)] * 6)
                                              ) for i in range(I)]
                                                                  ) for j in range(J)])
const_diff_ic.write(uss)
const_diff_ic.close()
#print(bessel(1, zm))
