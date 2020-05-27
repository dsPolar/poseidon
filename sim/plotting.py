import matplotlib as mp
mp.use('Agg')
import numpy as np
import matplotlib.pyplot as plt

def plot(x,y, title, xlabel, ylabel, save_path):
    plt.close()
    plt.figure(figsize=(10,8))
    if(y is not None):
        plt.plot(x,y, linewidth=0.75)
    else:
        plt.plot(x, linewidth=0.75)
    plt.title(title, fontsize=18)
    plt.xlabel(xlabel, fontsize=18)
    plt.ylabel(ylabel, fontsize=18)
    plt.savefig(save_path)
    plt.close()

def plot_ticks(x,y, title, xlabel, ylabel, save_path):
    plt.close()
    plt.figure(figsize=(10,8))
    ms = len(x) * 0.05
    s = ms/1000
    millis = np.arange(0,int(ms)+1,50)
    seconds = np.arange(0,int(s)+1,2)
    positions = millis/0.05
    plt.xticks(positions,millis)

    if(y is not None):
        plt.plot(x,y, linewidth=0.75)
    else:
        plt.plot(x, linewidth=0.75)



    plt.title(title)
    plt.xlabel(xlabel, fontsize=18)
    plt.ylabel(ylabel, fontsize=18)
    plt.savefig(save_path)
    plt.close()

def plot_width(x,y, title, xlabel, ylabel, save_path, width):
    plt.close()
    plt.figure(figsize=(20,16))
    if(y is not None):
        plt.plot(x,y, linewidth=width)
    else:
        plt.plot(x, linewidth=width)
    plt.title(title)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel(xlabel, fontsize=18)
    plt.ylabel(ylabel, fontsize=18)
    plt.savefig(save_path)
    plt.close()

def log_plot_x(x,y, title, xlabel, ylabel, save_path):
    plt.close()
    plt.figure(figsize=(20,16))
    if(y is not None):
        plt.plot(x,y, linewidth=1.5)
    else:
        plt.plot(x, linewidth=1.5)
    plt.title(title)
    plt.xscale("log")
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(save_path)
    plt.close()

def log_plot_y(x,y, title, xlabel, ylabel, save_path):
    plt.close()
    plt.figure(figsize=(20,16))
    if(y is not None):
        plt.plot(x,y, linewidth=1.5)
    else:
        plt.plot(x, linewidth=1.5)
    plt.title(title)
    plt.yscale("log")
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(save_path)
    plt.close()

def paulo_valid_plot(freqs, psd, spec, excit, inhib):
    plt.close()
    plt.figure(figsize=(15,12))
    plt.plot(freqs,psd,linewidth = 1.5)
    plt.plot(range(20, 1000, 20), spec,linewidth = 1.5)
    plt.plot(range(20, 1000, 20), excit,linewidth = 1.5)
    plt.plot(range(20, 1000, 20), inhib,linewidth = 1.5)
    plt.legend(["True", "Theoretical", "Theoretical (exc. only)", "Theoretical (inh. only)"])
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Power density (mV^2/Hz)")
    plt.yscale("log")
    plt.xscale("log")
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.savefig("graphs/paulograph.png")
    plt.close()


def log_plot(x,y, title, xlabel, ylabel, save_path):
    plt.close()
    plt.figure(figsize=(10,8))
    if(y is not None):
        plt.plot(x,y, linewidth=1.5)
    else:
        plt.plot(x, linewidth=1.5)
    plt.title(title)
    plt.xscale("log")
    plt.yscale("log")
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel(xlabel, fontsize=18)
    plt.ylabel(ylabel, fontsize=18)
    plt.savefig(save_path)
    plt.close()

def event(x, title, xlabel,ylabel, save_path):
    plt.close()
    plt.figure(figsize=(20,16))
    fig,ax = plt.subplots(1)

    plt.eventplot(x)

    plt.title(title)
    ax.set_yticklabels([])
    plt.xlabel(xlabel)
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.ylabel(ylabel)
    plt.savefig(save_path)
    plt.close()

def surface(x,y,z , title, xlabel, ylabel, save_path):
    plt.close()
    plt.figure(figsize=(15,12))
    fig,ax = plt.subplots()

    CS = ax.contour(x,y,z)
    ax.clabel(CS, inline=1, fontsize=10)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(save_path)
    plt.close()

def show(x, save_path):
    plt.close()
    fig,ax = plt.subplots()

    points = np.arange(0,31,5)
    positions = np.arange(0,61,10)
    plt.xticks(positions, points)
    plt.yticks(positions, points)

    plt.imshow(x, origin='lower')
    plt.title("Nu Error Space")
    plt.xlabel("inhibitory coded")
    plt.ylabel("excitatory coded")
    plt.colorbar()
    plt.savefig(save_path)
    plt.close()

def hist(x, title, xlabel, save_path):
    plt.close()
    plt.hist(x, bins=50)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel("Frequency")
    plt.savefig(save_path)
    plt.close()


def error_comps(psd_error, v_error, std_error, save_path):
    plt.close()
    fig,ax = plt.subplots()

    plt.plot(v_error)
    plt.plot(psd_error)
    plt.plot(std_error)
    plt.legend(['Mean V', 'PSD', 'Std Dev'])
    plt.title("Error components for estimated ratio of rates")
    plt.xlabel("range(1,501,1)/250")
    plt.ylabel("error value")
    plt.savefig(save_path)
    plt.close()

def overlap_psd(x,y, freqs, title, xlabel, ylabel, save_path):
    plt.close()
    plt.figure(figsize=(10,8))

    plt.plot(freqs, x, linewidth=1.5, label="Target")
    plt.plot(freqs, y, linewidth=1.5, label="Fitted Result")
    plt.legend()

    plt.title(title, fontsize=18)
    plt.xscale("log")
    plt.yscale("log")
    plt.xticks(fontsize=18)
    plt.yticks(fontsize=18)
    plt.xlabel(xlabel, fontsize=18)
    plt.ylabel(ylabel, fontsize=18)
    plt.savefig(save_path)
    plt.close()

def overlap_trace(x,y, title, xlabel, ylabel, save_path):
    plt.close()
    plt.figure(figsize=(10,8))

    plt.plot(x, linewidth=0.75)
    plt.plot(y, linewidth=0.75)

    plt.title(title, fontsize=18)
    plt.xlabel(xlabel, fontsize=18)
    plt.ylabel(ylabel, fontsize=18)
    plt.savefig(save_path)
    plt.close()
