import sys, subprocess, shlex, notify, os

# FORMAT FOR THE COMMAND FILE:
#
# number of istances : common parameters for all istances and algorthms
# alg1 parameters : name displayed on plot
# alg2 parameters : name displayed on plot
# ...

bot = notify.bot(profile="default")
bot.edit_profile(disable_notification=True)

try:

    file = sys.argv[1]

    with open(file, "r") as f:
        lines = f.read().splitlines()

    algs = lines[1:]

    list_cost = [[] for i in range(len(algs))]
    list_time = [[] for i in range(len(algs))]

    n_instances = int(lines[0].partition(":")[0].strip())

    print(f"Running {n_instances} from the {file} file.")

    bot.send_message_by_text(f"Algorithms: {len(algs)}\nInstances: {n_instances}")

    for alg in range(len(algs)):

        bot.create_progress_bar(n_instances, f"Algorithm: {algs[alg].partition(': ')[2]}")

        for i in range(1, n_instances + 1):
            command = f"./tsp -seed {i} {lines[0].partition(':')[2].strip()} -alg {algs[alg].partition(':')[0].strip()} -verbose 0 -noplot"
            print(f"Running command {command}\n")

            #bot.send_message_by_text(command)

            result = subprocess.run(shlex.split(command), capture_output=True, text = True).stdout

            print(result)

            result_cost = result.partition("BEST SOLUTION:")[2].partition("Cost:")[2].partition("\n")[0].strip()
            result_time = float(result.partition("BEST SOLUTION:")[2].partition("Execution time:")[2].partition("s\n")[0].strip())
            list_cost[alg].append(result_cost)
            list_time[alg].append(result_time)

            bot.update_progress_bar()

        bot.conclude_progress_bar()

    if not os.path.exists("plots"): os.mkdir("plots")
    file = file.replace("plotting", "plots")

    with open(f"{file}_cost_result.csv", "w") as f:

        f.write(f"{len(algs)}")
        for alg in range(len(algs)): f.write(f", {algs[alg].partition(':')[2].strip()}")
        f.write("\n")

        for i in range(1, n_instances + 1):
            f.write(f"{i}")
            for alg in range(len(algs)):
                f.write(f", {list_cost[alg][i-1]}")
            f.write("\n")

    with open(f"{file}_time_result.csv", "w") as f:

        f.write(f"{len(algs)}")
        for alg in range(len(algs)): f.write(f", {algs[alg].partition(':')[2].strip()}")
        f.write("\n")

        for i in range(1, n_instances + 1):
            f.write(f"{i}")
            for alg in range(len(algs)):
                f.write(f", {list_time[alg][i-1]}")
            f.write("\n")

    subprocess.run(shlex.split(f'python3 plotting/pp.py -X "Cost Ratio" -D , -S 2 {file}_cost_result.csv {file}_cost_result.png -P ""'))
    subprocess.run(shlex.split(f'python3 plotting/pp.py -X "Time Ratio" -D , -S 2 {file}_time_result.csv {file}_time_result.png -P ""'))

    bot.on()
    bot.send_photo_by_path(f"{file}_cost_result.png")
    bot.send_photo_by_path(f"{file}_time_result.png")

except Exception as e:

    bot.send_exception(type(e).__name__)
    raise e
