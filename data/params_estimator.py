import smr_read
import plotting
import scipy.io
from scipy import signal
import numpy as np
import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt


SMR_PATHS = ["JP36_002.smr","JP37_001.smr","JP37_003.smr","JP37_005.smr","JP38_003.smr","JP38_005.smr","JP41_006.smr","JP41_014.smr","JP43_005.smr","JP45_005.smr"]
THRESHOLD = -50
RESTING   = -70

def smr_describe(path, file):
    channel_list = smr_read.smr_read(path)
    file.write(str(path[10:18]) + "\n")
    i = 1
    for ch in channel_list:
        file.write(str(scipy.stats.describe(ch['signal'])) + "\n")

        signal_plot(ch['signal'], path[10:18], "data_graphs/Channel" + str(i) + "/" + path[10:18] + "Channel" + str(i) + ".svg")

        if i == 1:
            #ARRAY SAVING

            pickle_file = open("../arrays/" + path[10:18] + ".npy", "wb+")
            np.save(pickle_file, ch['signal'])
            pickle_file.close()

            #SPECTRA
            array = np.asarray(ch['signal'])
            trace = []
            for el in ch['signal']:
                trace.append(el[0])

            #freqs, times, spectrogram = signal.spectrogram(np.asarray(trace), fs=20000)

            #plotting.plot(spectrogram, None, "Spectrogram", "Time Window", "Frequency Band", "data_graphs/Spectrograms/spectro_" + path[10:18] + ".png" )

            freqs, psd = signal.welch(np.asarray(trace), scaling='spectrum', fs=20000)
            freqs1 = freqs
            psd1 = psd
            freqs = freqs[freqs <= 1000]
            psd = psd[:freqs.size]


            plotting.log_plot(freqs1, psd1, path[10:18] + " Power Spectrum", "Log Frequency/ Hz", "Log Power/ $mV^2Hz^{-1}$", "data_graphs/PSD/pxx_" + path[10:18] + ".png")
            plotting.plot_width(freqs, psd, path[10:18] + " Power Spectrum", "Frequency", "Power", "data_graphs/PSD/pxx_nonlog_" + path[10:18] + ".png", 1.5)
            #low, high = low_high(np.asarray(trace), THRESHOLD)
            #signal_plot(low, path[10:18] + " Low band ", "data_graphs/low/" + path[10:18] + "lowband.png")
            #signal_plot(high, path[10:18] + " High band ", "data_graphs/high/" + path[10:18] + "highband.png")




        i+= 1

def find_next_low(signal, i, thresh):
    i_start = i
    no_overflow = True
    if i > 0:
        start = signal[i-1]
    else:
        start = RESTING
    while signal[i] >= thresh:
        i+= 1
        if i >= signal.size:
            end = start
            no_overflow = False
            break
    if no_overflow:
        end = signal[i]
    return (start+end)/2.0, (i-i_start)

def low_high(signal, thresh):
    high = np.copy(signal)
    low  = np.copy(signal)
    for i in range(signal.size):
        if signal[i] >= thresh:
            high[i] = signal[i]
            mean, count = find_next_low(signal,i,thresh)
            for c in range(count):
                signal[i+c] = mean
        else:
            low[i] = signal[i]
            high[i] = RESTING
    return low, high



def data_describe():
    file = open("smr_descript.txt", "w+")
    for path in SMR_PATHS:
        smr_describe("wholecell/"+ path, file)
    file.close()

def signal_plot(signal, title, save_path):
    plt.close()
    plt.figure(figsize=(10,8))

    ms = len(signal)*0.05
    s  = ms / 1000

    seconds = np.arange(0,int(s),60)
    positions = (seconds*1000)/0.05

    plt.xticks(positions,seconds)


    plt.plot(signal, linewidth=0.5)
    plt.title(title)




    plt.xlabel("time/s")
    plt.ylabel("Membrane Potential/mV")
    plt.savefig(save_path)

def spectrum_plot(f, pxx, title, save_path):
    plt.close()
    plt.plot(f, pxx, linewidth=0.5)
    plt.title(title)
    plt.xlabel("Frequency/Hz")
    plt.ylabel("Power/mV^2")
    plt.savefig(save_path)


data_describe()
