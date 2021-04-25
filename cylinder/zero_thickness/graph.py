import matplotlib.pyplot as plt
import numpy as np


with open("output.txt", "r") as f:
    steps, l = tuple(map(int, f.readline().split()))
    conc = np.empty(shape=(steps, l, 4))
    for i in range(steps):
        for j in range(l):
            vec = list(map(float, f.readline().split(",")))
            conc[i, j, :] = vec


for component in range(4):
    h = conc[:, :, component]
    fig = plt.figure(figsize=(15, 10))
    x = np.array([[i] * l for i in range(steps)]).reshape(-1)
    y = [j for j in range(l)] * steps

    plt.hist2d(x, y, bins=(steps, l), weights=h.reshape(-1), cmap="plasma")
    plt.colorbar()

    plt.xlabel("time steps")
    plt.ylabel("space steps")
    plt.savefig(f"graph{component+1}.png")
