import matplotlib.pyplot as plt
import pandas as pd
import sys

# read file
fname = sys.argv[1]
data = pd.read_csv(fname, names=['X(t)/X(0)'])

# plot data
plot = data.plot.line(title="Autocorrelation")
fig = plot.get_figure()
plt.xlabel('Time steps')
plt.ylabel('X(t)/X(0)')
#plt.ylim([-0.5, 1.5])
plt.grid()
plt.show()
fig.savefig(fname[:-4] + ".png")
