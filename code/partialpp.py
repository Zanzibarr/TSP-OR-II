import shlex
import subprocess
import notify

# type = sys.argv[1]

file = "plots/pp_whole_"
cost_file, time_file = file + "cost_result.csv", file + "time_result.csv"
end_location = "plots/"

# with open(cost_file, "r") as f:
#    costs = f.read()

# with open(time_file, "r") as f:
#    times = f.read()

# costs = costs.splitlines()
# times = times.splitlines()

# costs_m = []
# times_m = []

# for line in costs:
#    costs_m.append([x.strip() for x in line.split(",")])
# for line in times:
#    times_m.append([x.strip() for x in line.split(",")])

# if type == "heur":
#    start = 0
#    n_algs = 4
#    name = "perfprof_heur"

# if type == "tabu":
#    start = 4
#    n_algs = 5
#    name = "perfprof_tabu"

# if type == "vns":
#    start = 9
#    n_algs = 2
#    name = "perfprof_vns"

# if type == "benders":
#    start = 11
#    n_algs = 2
#    name = "perfprof_benders"

# if type == "bnc":
#    start = 13
#    n_algs = 4
#    name = "perfprof_bnc"

# if type == "bnc_patch":
#    start = 16
#    n_algs = 4
#    name = "perfprof_bnc_patch"

# if type == "hard":
#    start = 22
#    n_algs = 4
#    name = "perfprof_hard"

# if type == "lbv1":
#    start = 26
#    n_algs = 3
#    name = "perfprof_lbv1"

# if type == "lbv2":
#    start = 29
#    n_algs = 3
#    name = "perfprof_lbv2"

# with open(f"{end_location}{name}_costs.csv", "w") as f:
#    f.write(f"{n_algs}, {", ".join(costs_m[0][start+1 : start+1+n_algs])}\n")
#    for i in range(1, len(costs)):
#        f.write(f"{i}, {", ".join(costs_m[i][start+1 : start+1+n_algs])}\n")

# with open(f"{end_location}{name}_times.csv", "w") as f:
#    f.write(f"{n_algs}, {", ".join(times_m[0][start+1 : start+1+n_algs])}\n")
#    for i in range(1, len(times)):
#        f.write(f"{i}, {", ".join(times_m[i][start+1 : start+1+n_algs])}\n")

bot = notify.bot(profile="silent")

for name in [
    "perfprof_heur",
    "perfprof_tabu",
    "perfprof_vns",
    "perfprof_benders",
    "perfprof_bnc",
    "perfprof_bnc_patch",
    "perfprof_hard",
    "perfprof_lbv1",
    "perfprof_lbv2",
]:
    subprocess.run(
        shlex.split(
            f'python3 plotting/pp.py -X "Cost Ratio" -D , -S 2 {end_location}{name}_costs.csv {end_location}{name}_costs.png -P ""'
        )
    )
    subprocess.run(
        shlex.split(
            f'python3 plotting/pp.py -X "Time Ratio" -D , -S 2 {end_location}{name}_times.csv {end_location}{name}_times.png -P ""'
        )
    )

    bot.send_photo_by_path(f"{end_location}{name}_costs.png")
    bot.send_photo_by_path(f"{end_location}{name}_times.png")
