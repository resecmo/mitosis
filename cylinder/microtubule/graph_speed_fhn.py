import matplotlib.pyplot as plt
import numpy as np


with open("output.txt", "r") as f:
    total_time, total_len = tuple(map(int, f.readline().split()))
    steps, l = tuple(map(int, f.readline().split()))
    conc = np.empty(shape=(steps, l, 4))
    for i in range(steps):
        for j in range(l):
            vec = list(map(float, f.readline().split(",")))
            conc[i, j, :] = vec


for component in range(2):
    h = conc[:, :, component]

    t1 = 400
    x1 = t1 / steps * total_time
    max1 = np.argmax(h[t1, :])
    y1 = max1 / l * total_len

    t2 = 1490
    x2 = t2 / steps * total_time
    max2 = np.argmax(h[t2, :])
    y2 = max2 / l * total_len

    speed = (y2 - y1)/(x2 - x1); 
    print(f"v{component} = {speed}")

    fig = plt.figure(figsize=(15, 10))
    x = np.array([[i] * l for i in range(steps)]).reshape(-1) / steps * total_time
    y = np.array([j for j in range(l)] * steps) / l * total_len

    plt.hist2d(x, y, bins=(steps, l), weights=h.reshape(-1), cmap="plasma")
    plt.colorbar()
    plt.scatter([x1, x2], [y1, y2], c='c')
    plt.plot([x1, x2], [y1, y2], c='c', label=f"v={speed:.5f}")


    plt.xlabel("t")
    plt.ylabel("x")
    plt.legend()
    plt.savefig(f"graph{component+1}.png")
