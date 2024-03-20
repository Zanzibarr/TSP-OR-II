import matplotlib.pyplot as plt
import numpy as np
import sys

best_x = 0
best_y = 999999
x_width = 0
y_width = 0

if len(sys.argv) == 1:
    a = 0
    b = 16
elif len(sys.argv) == 2:
    a = int(sys.argv[1])
    b = a + 1
elif len(sys.argv) == 3:
    a = int(sys.argv[1])
    b = int(sys.argv[2])

for i in range(a, b):
    with open(f"int_costs_{i}.txt", "r") as f:
        sol = f.read()
        
    lines = sol.splitlines()
    
    if len(lines) == 0: continue

    y = []

    for line in lines:
        y.append(float(line))

    if best_y > min(y):
        best_y = min(y)
        best_x = np.argmin(y)

    if len(y) > x_width: x_width = len(y)
    if max(y) > y_width: y_width = max(y)
        
    plt.plot(range(0, len(y)), y, label=i, linewidth=.6)

plt.plot([0, x_width], [best_y, best_y], color = "red", linewidth=.5)
plt.plot([best_x, best_x], [best_y, y_width], color = "red", linewidth=.5)
plt.legend()
plt.savefig("plot.png")
plt.show()