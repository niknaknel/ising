import matplotlib.pyplot as plt
import pandas as pd
import sys

# read file
fname = sys.argv[1]
data = pd.read_csv(fname)
L = (fname.split('_')[2])[1:-4]

# plot data
#plot = data.plot(x ='T', y='mean_Mps', kind = 'scatter', title="Phase diagram")
#fig = plot.get_figure()


fig = plt.figure()

plt.errorbar('T', 'mean_Mps', yerr='std_Mps', fmt='o', data=data,
		ecolor='red', capsize=3, markersize=4)

plt.title('Phase diagram for L = %s' % L)
plt.xlabel('Temperature')
plt.ylabel('Magnetization per site')
plt.ylim([-0.2, 1.2])
plt.grid()
plt.show()
fig.savefig(fname[:-4] + ".png")
