import matplotlib as mp
mp.use('Agg')
import matplotlib.pyplot as plt

def plot(x,y, title, xlabel, ylabel, save_path):
    plt.close()
    plt.figure(figsize=(10,8))
    if(y is not None):
        plt.plot(x,y, linewidth=0.5)
    else:
        plt.plot(x, linewidth=0.5)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(save_path)

def plot_width(x,y, title, xlabel, ylabel, save_path, width):
    plt.close()
    plt.figure(figsize=(10,8))
    if(y is not None):
        plt.plot(x,y, linewidth=width)
    else:
        plt.plot(x, linewidth=width)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(save_path)

def log_plot_x(x,y, title, xlabel, ylabel, save_path):
    plt.close()
    plt.figure(figsize=(10,8))
    if(y is not None):
        plt.plot(x,y, linewidth=0.5)
    else:
        plt.plot(x, linewidth=0.5)
    plt.title(title)
    plt.xscale("log")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(save_path)

def log_plot_y(x,y, title, xlabel, ylabel, save_path):
    plt.close()
    plt.figure(figsize=(10,8))
    if(y is not None):
        plt.plot(x,y, linewidth=0.5)
    else:
        plt.plot(x, linewidth=0.5)
    plt.title(title)
    plt.yscale("log")
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(save_path)

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
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(save_path)

def event(x, title, xlabel,ylabel, save_path):
    plt.close()
    plt.figure(figsize=(10,8))
    fig,ax = plt.subplots(1)

    plt.eventplot(x)

    plt.title(title)
    ax.set_yticklabels([])
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.savefig(save_path)
