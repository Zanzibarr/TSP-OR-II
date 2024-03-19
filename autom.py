import subprocess, shlex

seed_list = [123, 492, 15, 99, 2122187, 2000139, 55]
node_list = [300]
alg_list = ["g2opt-best -mt", "tabu -tl 360"]

for seed in seed_list:
    for node in node_list:
        for alg in alg_list:
            subprocess.run(shlex.split(f"./main -seed {seed} -nodes {node} -alg {alg}"))