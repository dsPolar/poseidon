import plotting
import math
import numpy as np

bad_files = ["JP36_002","JP37_001","JP45_005"]
files = ["JP37_003","JP37_005","JP38_003","JP38_005","JP41_006","JP41_014","JP43_005"]#, "JP45_005", "JP36_002","JP37_001"]
endings = ["_two_out.txt", "_six_out.txt", "_eight_out.txt"]

if __name__ == "__main__":
    for file in files:
        for end in endings:
            try:
                f = open("jon_error/" + file + end, "r+")
                arr = []
                for line in f.readlines():
                    if((line.find("[") == -1) and (line.find("]") == -1)):
                        if(line.find(".") > -1):
                            arr.append(float(line[:line.find("\\")]))
                plotting.hist(arr,file + " Error Histogram", "Final Error Value", "error_graphs/" + file + end[0:4] + ".png")
                f.close()
            except:
                print(file + end + " did not exist")

    for end in endings:
        count = 0
        arr = []
        for file in files:
            try:
                f = open("jon_error/" + file + end, "r+")
                for line in f.readlines():
                    if((line.find("[") == -1) and (line.find("]") == -1)):
                        if(line.find(".") > -1):
                            arr.append(float(line[:line.find("\\")]))
                            count +=1
                f.close()
            except:
                print(file + end + " did not exist")
        type = end[1:4]
        if type == "eig":
            type = "eight"
        plotting.hist(arr, "All Trace Error Histogram for " + type.capitalize() + " Parameters", "Final Error Value", "error_graphs/" + type + ".png")
        mean = sum(arr) / len(arr)
        print(type + " Mean : " + str(mean))
        variance = sum([((x - mean) ** 2) for x in arr]) / len(arr)
        print(type + " Standard Deviation : " + str(math.sqrt(variance)))
        print("Across " + str(count) + " trials")

    f = open("error_out_error.txt", "r+")
    arr = []
    for line in f.readlines():
        if(line.find(".")>-1):
            arr.append(float(line[:line.find("\\")]))
    plotting.hist(arr, "Simulation Two Parameter Optimisation Error", "Final Error Value", "error_graphs/sim_error_hist_two.png")
    means = sum(arr) / len(arr)
    print("Simulation Mean : " + str(means))
    variance = sum([((x - means) ** 2) for x in arr]) / len(arr)
    print("Simulation Standard Deviation : " + str(math.sqrt(variance)))
    f.close()

    f = open("six_error_out_error.txt", "r+")
    arr = []
    for line in f.readlines():
        if(line.find(".")>-1):
            arr.append(float(line[:line.find("\\")]))
    plotting.hist(arr, "Simulation Six Parameter Optimisation Error", "Final Error Value", "error_graphs/sim_error_hist_six.png")
    means = sum(arr) / len(arr)
    print("Simulation Mean : " + str(means))
    variance = sum([((x - means) ** 2) for x in arr]) / len(arr)
    print("Simulation Standard Deviation : " + str(math.sqrt(variance)))
    f.close()

    f= open("jon_error/repeat_variance.txt", "r+")
    arr = []
    for line in f.readlines():
        if(line.find(".")>-1):
            arr.append(float(line[:line.find("\\")]))
    plotting.hist(arr, "Spread of Error for same parameters", "Final Error Value", "error_graphs/repeat_hist.png")
    means = sum(arr) / len(arr)
    print("Repeat mean : " + str(means))
    variance = sum([((x - means) ** 2) for x in arr]) / len(arr)
    print("Repeat Standard Deviation : " + str(math.sqrt(variance)))
    f.close()
