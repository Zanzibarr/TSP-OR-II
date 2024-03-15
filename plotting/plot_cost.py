import sys
import matplotlib.pyplot as plt

with open(sys.argv[1], "r") as f:
    sol = f.read()
    
lines = sol.splitlines()

y = []

for line in lines:
    y_ = line
    y.append(float(y_))
    
plt.plot(range(0, len(y)), y)
plt.show()