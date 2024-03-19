import sys
import matplotlib.pyplot as plt
import numpy as np

best_x = 0
best_y = 999999
x_width = 0
y_width = 0

for i in range(0, 16):
    with open(f"test{i}", "r") as f:
        sol = f.read()
        
    lines = sol.splitlines()

    y = []

    for line in lines:
        y.append(float(line))

    if best_y > min(y):
        best_y = min(y)
        best_x = np.argmin(y)

    if len(y) > x_width: x_width = len(y)
    if max(y) > y_width: y_width = max(y)
        
    plt.scatter(range(0, len(y)), y, s=1)

    with open(f"test{i}", "w") as f:
        f.write("")

#best_best = min(best)
plt.plot([0, x_width], [best_y, best_y], color = "red")
plt.plot([best_x, best_x], [best_y, y_width], color = "red")
plt.savefig("plot.png")
#plt.show()