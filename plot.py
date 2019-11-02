import matplotlib.pyplot as plt
import pandas as pd
import sys

# read file
fname = sys.argv[1]
data = pd.read_csv(fname)

# plot data
plot = data.plot.line(title="Plot of M per site for %s" % fname)
fig = plot.get_figure()
plt.xlabel('Time steps')
plt.ylabel('Magnetization per site')
plt.ylim([-1.2, 1.2])
plt.grid()
fig.savefig(fname[:-4] + ".png")
