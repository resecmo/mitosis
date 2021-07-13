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
w = conc[-1, :, 1]
u_between = 0.5 * (u[:-1] + u[1:])
w_between = 0.5 * (w[:-1] + w[1:])
u_prime = np.diff(u) / h
w_prime = np.diff(w) / h

denominator = np.trapz(u_prime**2) * h + np.trapz(w_prime**2) * h

I1 = np.trapz(w_between * u_prime) * h
I2 = np.trapz(u_between * w_prime) * h
numerator = - I1 + I2

v = numerator / denominator
print(f"v = {numerator}/{denominator} = {v}")

