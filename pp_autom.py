import sys, subprocess, shlex, notify, os

# FORMAT FOR THE COMMAND FILE:
#
# number of istances : common parameters for all istances and algorthms
# alg1 parameters : name displayed on plot
# alg2 parameters : name displayed on plot
# ...

time = ("-time" in sys.argv)

bot = notify.bot(profile="default")
bot.edit_profile(disable_notification=True)

try:

    file = sys.argv[1]

    with open(file, "r") as f:
        lines = f.read().splitlines()

    algs = lines[1:]

    list = [[] for i in range(len(algs))]

    n_instances = int(lines[0].partition(":")[0].strip())

    print(f"Running {n_instances} from the {file} file.")

    bot.send_message_by_text(f"Algorithms: {len(algs)}\nInstances: {n_instances}")

    for alg in range(len(algs)):

        bot.create_progress_bar(n_instances, f"Algorithm: {algs[alg].partition(': ')[2]}")

        for i in range(1, n_instances + 1):
            command = f"./tsp -seed {i} {lines[0].partition(':')[2].strip()} -alg {algs[alg].partition(':')[0].strip()} -verbose 0"
            print(f"Running command {command}\n")

            #bot.send_message_by_text(command)

            result = subprocess.run(shlex.split(command), capture_output=True, text = True).stdout

            print(result)

            if not time: result = result.partition("BEST SOLUTION:")[2].partition("Cost:")[2].partition("\n")[0].strip()
            else: result = float(result.partition("BEST SOLUTION:")[2].partition("Execution time:")[2].partition("s\n")[0].strip())
            list[alg].append(result)

            print(result)

            bot.update_progress_bar()

        bot.conclude_progress_bar()

    if not os.path.exists("plots"): os.mkdir("plots")
    file = file.replace("plotting", "plots")

    with open(f"{file}_result.csv", "w") as f:

        f.write(f"{len(algs)}")
        for alg in range(len(algs)): f.write(f", {algs[alg].partition(':')[2].strip()}")
        f.write("\n")

        for i in range(1, n_instances + 1):
            f.write(f"{i}")
            for alg in range(len(algs)):
                f.write(f", {list[alg][i-1]}")
            f.write("\n")

    subprocess.run(shlex.split(f'python3 plotting/pp.py -X "{"Time Ratio" if time else "Cost Ratio"}" -D , -S 2 {file}_result.csv {file}_result.png -P ""'))

    bot.on()
    bot.send_photo_by_path(f"{file}_result.png")

except Exception as e:

    bot.send_exception(type(e).__name__)
    raise e
