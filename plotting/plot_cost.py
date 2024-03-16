import sys
import matplotlib.pyplot as plt
import numpy as np

for i in range(0, 12):
    with open(f"test{i}", "r") as f:
        sol = f.read()
        
    lines = sol.splitlines()

    y = []

    for line in lines:
        y.append(float(line))

    best = np.argmin(y)
        
    plt.scatter(range(0, len(y)), y, s=1)

    with open(f"test{i}", "w") as f:
        f.write("")

plt.show()