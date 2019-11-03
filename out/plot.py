import matplotlib.pyplot as plt
import pandas as pd
import sys

# read file
fname = sys.argv[1]
data = pd.read_csv(fname)

# plot data
plot = data.plot.line(x ='T', y='tau', style = '-o', title="Autocorrelation at different temperatures.")
fig = plot.get_figure()
plt.xlabel('Temperature')
plt.ylabel('Autocorrelation time')
#plt.ylim([-0.2, 1.2])
plt.grid()
plt.show()
fig.savefig(fname[:-4] + ".png")
