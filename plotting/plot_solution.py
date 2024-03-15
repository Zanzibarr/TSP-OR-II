import sys
import matplotlib.pyplot as plt

with open(sys.argv[1], "r") as f:
    sol = f.read()
    
lines = sol.splitlines()

alg         = lines[0].partition("Algorithm: ")[2].strip()
cost        = lines[1].partition("Cost: ")[2].strip()
time        = lines[2].partition("Time: ")[2].strip()
total_time  = lines[3].partition("Total execution time: ")[2].strip()
over_time   = "has" in lines[4]

x, y = [], []

for line in lines[5:]:
    x_, y_ = line.split(" ")
    x.append(float(x_))
    y.append(float(y_))
    
plt.plot(x, y)
plt.plot(x[0], y[0], ".", markersize=14, color="blue")
plt.figtext(0.5, 0.92, f"{alg} - nodes: {len(x)-1} - cost: {cost} - total time: {total_time}", wrap=True, horizontalalignment='center', fontsize=12)
plt.show()