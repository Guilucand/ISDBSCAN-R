import matplotlib.pyplot as plt
import numpy as np
from cycler import cycler
from matplotlib.animation import FuncAnimation
import os


def show_plot(file_name, cluster_name, directory, interactive=False):
    data = [[float(x[0].strip()), float(x[1].strip())] for x in
            [y.split() for y in open(file_name, "r").readlines()]]
    data = np.array(list(map(list, zip(*data))))

    clusters = np.array([int(x.strip()) for x in open(cluster_name).readlines()])

    fig = plt.figure()
    ax = fig.add_axes([0, 0, 1, 1])

    clus_num = np.max(clusters)

    NUM_COLORS = (clus_num + 1)
    SIZE = 10
    cm = plt.get_cmap('gist_rainbow')

    def plot_cluster(i):
        if i > clus_num:
            return
        ax.scatter(
            data[0][clusters == i],
            data[1][clusters == i],
            s=SIZE)

    def init_ax():
        ax.set_xlabel('Grades Range')
        ax.set_ylabel('Grades Scored')
        ax.set_title('Blank plot' + file_name.split("/")[-1][:-4])
        cy = cycler('color', [cm(1. * i / NUM_COLORS) for i in range(NUM_COLORS)])
        ax.set_prop_cycle(cy)

    init_ax()
    ax.scatter(data[0][clusters == -2], data[1][clusters == -2], s=SIZE, color='black')

    print(f"Noise points count: {np.sum(clusters == -2)}")

    for i in range(0, clus_num + 1):
        plot_cluster(i)

    try:
        os.mkdir(directory)
    except Exception:
        pass

    fig.savefig(f'{directory}/{cluster_name.split("/")[-1][:-4]}.png')

    if interactive:
        ax.cla()
        init_ax()
        ax.scatter(data[0][clusters == -2], data[1][clusters == -2], s=SIZE, color='black')
        ani = FuncAnimation(fig, plot_cluster, interval=250)
        plt.show()
    # c=clusters[clusters == i],
    # )


if __name__ == '__main__':
    file_name = "../data/can473.txt"
    cluster_name = "../data/can473.clustering.txt"
    show_plot(file_name, cluster_name)
