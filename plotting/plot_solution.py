import sys, os, notify
import matplotlib.pyplot as plt

with open(sys.argv[1], "r") as f:
    sol = f.read()

#TODO: Parse also the info of the solution to show in the plot

lines = sol.splitlines()

for line in lines[lines.index("--------------------")+1:]:
    nfrom, nto = line.split(" -> ")
    fromx, fromy = nfrom.partition("(")[2].partition(")")[0].split(",")
    tox, toy = nto.partition("(")[2].partition(")")[0].split(",")
    xlist = (float(fromx), float(tox))
    ylist = (float(fromy), float(toy))
    plt.plot(xlist, ylist, color="tab:blue")

#plt.show()
if not os.path.exists("../plots"): os.mkdir("../plots")
plt.savefig(f"../plots/{os.path.basename(sys.argv[1]).split('/')[-1].replace('.txt', '_plot.png')}")

path = f"../plots/{os.path.basename(sys.argv[1]).split('/')[-1].replace('.txt', '_plot.png')}"

notify.bot(profile="default").send_photo_by_path(path, caption=sys.argv[1])