import sys, os
import matplotlib.pyplot as plt

with open(sys.argv[1], "r") as f:
    sol = f.read()
    
lines = sol.splitlines()
x, y = [[]], [[]]
loop = -1

for line in lines[lines.index("--------------------")+1:]:
    if "LOOP" in line:
        loop += 1
        x.append([])
        y.append([])
        continue
    if loop == -1: loop = 0
    _, x_, y_ = line.split()
    x[loop].append(float(x_))
    y[loop].append(float(y_))
    
for i in range(0, loop+1):
    plt.plot(x[i], y[i], color="tab:blue")

if loop == 0:
    plt.plot(x[0][0], y[0][0], ".", markersize=14, color="blue")

#plt.show()
if not os.path.exists("../plots"): os.mkdir("../plots")
plt.savefig(f"../plots/{os.path.basename(sys.argv[1]).split('/')[-1].replace("file.txt", "plot.png")}")
