import matplotlib.pyplot as plt
import pandas as pd
import sys

# read file
fname = sys.argv[1]
data = pd.read_csv(fname)

# plot data
vals = fname.split('_')
T = (vals[1])[1:]
L = (vals[2])[1:]
teq = int((vals[3])[3:-4])

title = "Equilibration of Magnetisation at T=%s with L=%s" % (T, L)

plot = data.plot(kind = 'line', title=title)
fig = plot.get_figure()

plt.xlabel('Time steps')
plt.ylabel('Magnetization per site')
plt.ylim([-1.2, 1.2])
plt.axvline(x=teq, color='red')
plt.grid()

#plt.show()
fig.savefig(fname[:-4] + ".png")
