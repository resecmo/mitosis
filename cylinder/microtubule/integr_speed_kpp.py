import matplotlib.pyplot as plt
import numpy as np


with open("output.txt", "r") as f:
    total_time, total_len = tuple(map(float, f.readline().split()))
    steps, l = tuple(map(int, f.readline().split()))
    conc = np.empty(shape=(steps, l, 4))
    for i in range(steps):
        for j in range(l):
            vec = list(map(float, f.readline().split(",")))
            conc[i, j, :] = vec

h = total_len / l

u = conc[-1, :, 0]

u_prime = np.diff(u) / h

denominator = np.trapz(u_prime**2) * h

def numer_int(u):
    return 2.0 * (u**2 / 2 - u**3 / 3)
numerator = numer_int(1) - numer_int(0)

v = numerator / denominator
print(f"v = {numerator}/{denominator} = {v}")


