import random

genome_filename = input("Enter name of .gb file you will be using: ")
output_prefix = input("Enter prefix for output files: ")

run_num = 0

f = open("commands.sh", "w")

for i in range(96):
    run_num += 1
    seed = random.randint(0,4294967295)
    command = "python3 " + "$WORK/src/igem2020/T7-simulation/phage_model.py" + " -i " + genome_filename + " -o " + output_prefix + "_" + str(run_num) + " -s " + str(seed) + " -r\n"
    f.write(command)

f.close()
