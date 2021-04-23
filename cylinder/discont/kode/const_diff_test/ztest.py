import math
I=100
J=100
hz=0.01
ht=0.01
d=0.1
a=math.sqrt(d)
def sol(z, t):
    if t > 0:
        res = math.exp(-z**2 / (4 * a**2 * t)) / (2 * a * math.sqrt(math.pi * t))
        return res if res > 1e-200 else 0
    else:
        return (z == 0) * 2 / hz

folder = "pres/zdiff0/"
const_diff_ic = open(folder + "zconst_diff_ic.scv", "w+")
t0=10*ht
uss='\n'.join([';'.join(
                  [','.join(
                            map(str, [sol(j*hz, t0)] * 6)
                                              ) for i in range(I)]
                                                                  ) for j in range(J)])
const_diff_ic.write(uss)
const_diff_ic.close()
#print(bessel(1, zm))
each = 1
time = 0.1
const_diff_ex = open(folder + "zconst_diff_ex.csv", "w+")
n_rows = int(time/ht/each+1)
print("n_rows = %d\n" %n_rows)
times = [t0 + i*ht*each for i in range(0,n_rows)]
def sol_zs(t):
    return ','.join([str(sol(j*hz, t)) for j in range(0,J)])
uss='\n'.join(
        map(sol_zs, times))
const_diff_ex.write(uss)
const_diff_ex.close()
