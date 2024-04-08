import sys, os
import matplotlib.pyplot as plt

with open(sys.argv[1], "r") as f:
    sol = f.read()
    
lines = sol.splitlines()
x, y = [], []

for line in lines[lines.index("--------------------")+1:]:
    _, x_, y_ = line.split()
    x.append(float(x_))
    y.append(float(y_))
    
plt.plot(x, y)
plt.plot(x[0], y[0], ".", markersize=14, color="blue")

#plt.show()
plt.savefig(f"../plots/{os.path.basename(sys.argv[1]).split('/')[-1].replace("file.txt", "plot.png")}")
