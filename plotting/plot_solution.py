import sys, os, notify
import matplotlib.pyplot as plt

sol_file = sys.argv[1]

with open(sol_file, "r") as f:
    sol = f.read()

#Parse the info to plot
title = sol.partition("Algorithm: ")[2].partition("\n")[0]

parameters = sol.partition("\n")[2].partition("Number of components: ")[0].strip()

swap = ""
f2opt = False
tenure = ""
variability = ""
frequency = ""
fvns = False
mipstart = False
bloop = False
patching = ""
cand_cb = False
rel_cb = False

incomplete = False

#FIXME: Parsing the new patching lines

for line in parameters.splitlines():
    if "Swap policy:" in line: swap = line.partition("Swap policy: ")[2].partition(" swap")[0]
    elif "f2opt" in line: f2opt = True
    elif "Tabu tenure:" in line: tenure = line.partition("Tabu tenure: ")[2].partition(".\n")[0]
    elif "Tabu variability:" in line: variability = line.partition("Tabu variability: ")[2].partition(".\n")[0]
    elif "Tabu variability frequency:" in line: frequency = line.partition("Tabu tenuvariability frequencyre: ")[2].partition(".\n")[0]
    elif "Fast vns enabled" in line: fvns = True
    elif "Using a mipstart" in line: mipstart = True
    elif "Using benders loop" in line: bloop = True
    elif "normal patching" in line: patching = "normal"
    elif "greedy patching" in line: patching = "greedy"
    elif "candidate callback" in line: cand_cb = True
    elif "relaxation callback" in line: rel_cb = True

#g2opt
if title == "g2opt":
    if swap == "first": title += " (f)"
    elif swap == "best": title += " (b)"
    elif f2opt: title = "f2opt"

#tabu
if title == "tabu":
    title += f" ({tenure}-{variability}-{frequency})"

#vns
if title == "vns":
    if fvns: title = "fvns"

#cplex
if title == "cplex":
    if bloop: title = "benders loop"
    if mipstart: title += " mipst"
    if cand_cb: title += " ccb"
    if rel_cb: title += " rcb"
    if patching == "normal": title += f" n-patch"
    if patching == "greedy": title += f" g-patch"

results = sol.partition("Number of ")[2].partition("--------------------\n")[0].splitlines()

#incomplete
if any(check in results[-1] for check in ["time limit", "terminated"]):
    title += "*"

title += "\n"

for result in results:
    if "components: " in result:
        ncomp = result.partition("components: ")[2].strip()
        if ncomp != "1":
            title += ncomp + " cycles, "
    elif "Cost: " in result: title += f"cost: {result.partition('Cost: ')[2].strip()}, "
    elif "Time: " in result: title += f"{result.partition('Time: ')[2].strip()}/"
    elif "Total execution time: " in result: title += result.partition("Total execution time: ")[2].strip()

lines = sol.splitlines()
split = lines.index("--------------------")

for line in lines[split+1:]:
    nfrom, nto = line.split(" -> ")
    fromx, fromy = nfrom.partition("(")[2].partition(")")[0].split(",")
    tox, toy = nto.partition("(")[2].partition(")")[0].split(",")
    xlist = (float(fromx), float(tox))
    ylist = (float(fromy), float(toy))
    plt.plot(xlist, ylist, color="tab:blue")

plt.title(title)

#plt.show()

if not os.path.exists("../plots"): os.mkdir("../plots")
path = f"../plots/{os.path.basename(sol_file).split('/')[-1].replace('.txt', '_plot.png')}"
plt.savefig(path)

notify.bot(profile="default").send_photo_by_path(path, caption=sol_file)