import numpy as np
from numpy import sqrt
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import ListedColormap
from ast import literal_eval as lit
import sys

colors = ['green', 'red']
cmap = ListedColormap(colors)

def init(file_name):
    states = []

    with open(file_name) as fp:
        # read first line
        line = fp.readline().strip()
        s = np.array(lit(line))

        # get constants
        N = len(s)
        L = int(sqrt(N))

        # display initial state
        grid = np.reshape(s, (L, L))
        states.append(grid)

        while line:
            grid = np.reshape(s, (L, L))
            states.append(grid)

            line = fp.readline().strip()
            if line == "":
                break
            else:
                s = np.array(lit(line))

    return states


def visualize(file_name):
    fig = plt.figure(figsize=(8, 8))
    states = init(file_name)
    ims = [[plt.imshow(s, interpolation='none', cmap=cmap)] for s in states]

    ani = animation.ArtistAnimation(fig, ims, interval=10, blit=True,
                                    repeat_delay=1000)

    ani.save(file_name[:-4] + "_anim.mp4")


if __name__ == "__main__":
    fname = sys.argv[1]
    visualize(fname)
