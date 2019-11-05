import matplotlib.pyplot as plt
import pandas as pd
import sys

# read file
fname = sys.argv[1]
data = pd.read_csv(fname)

# plot data
plot = data.plot(x ='T', y='Mps', kind = 'scatter', title="Phase diagram")
fig = plot.get_figure()
plt.xlabel('Temperature')
plt.ylabel('Magnetization per site')
plt.ylim([-0.2, 1.2])
plt.grid()
plt.show()
fig.savefig(fname[:-4] + ".png")
