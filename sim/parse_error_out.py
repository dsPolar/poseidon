path = "error_out.txt"



def parse_error():
    true_array = []
    guess_array = []
    f = open(path, "r")
    for line in f.readlines():
        if(line.find(".") > 0):
            a_end = line.rfind(".") - 2
            b_start = line.rfind(".") - 1
            if(line.find("]") == -1):
                a_true = float(line[0:a_end])
                b_true = float(line[b_start:])
                true_array.append((a_true,b_true))
            else:
                a_guess = float(line[1:a_end])
                b_guess = float(line[b_start:-2])
                guess_array.append((a_guess, b_guess))
    f.close()
    return true_array, guess_array

def mean_relative_error(true_array, guess_array):
    assert len(true_array) == len(guess_array)
    a_error = []
    b_error = []
    for i in range(len(true_array)):
        a_error.append((guess_array[i][0] - true_array[i][0]) / true_array[i][0])
        b_error.append((guess_array[i][1] - true_array[i][1]) / true_array[i][1])
    a_mean = sum(a_error) / len(a_error)
    b_mean = sum(b_error) / len(b_error)
    return (a_mean, b_mean), a_error, b_error

def mean_ratio_error(true_array, guess_array):
    assert len(true_array) == len(guess_array)
    a_ratio_error = []
    for i in range(len(true_array)):
        a_ratio_error.append(((guess_array[i][0] / guess_array[i][1]) - (true_array[i][0] / true_array[i][1])) / (true_array[i][0] / true_array[i][1]))
    a_ratio_mean = sum(a_ratio_error) / len(a_ratio_error)
    return a_ratio_mean, a_ratio_error

def categorise_error(a_error, b_error, ratio_errors):
    low_abs = 0
    low_rat = 0
    for i in range(len(a_error)):
        mean_absolute = (a_error[i] + b_error[i]) / 2
        print(str(i))
        if(abs(mean_absolute) > 0.25):
            print("High Absolute Error")
        else:
            print("Low Absolute Error")
            low_abs += 1

        if(abs(ratio_errors[i]) > 0.25):
            print("High Ratio Error")
        else:
            print("Low Ratio Error")
            low_rat += 1
    print("Percent low absolute " + (str(low_abs / len(a_error))))
    print("Percent low ratio " + (str(low_rat / len(a_error))))



if __name__ == "__main__":
    true_array,guess_array = parse_error()
    means, a_error, b_error = mean_relative_error(true_array, guess_array)
    ratio_mean, ratio_errors = mean_ratio_error(true_array, guess_array)
    print("Absolute Mean Error " + str(means[0]) + str(means[1]))
    print(str(ratio_mean))
    categorise_error(a_error,b_error,ratio_errors)
