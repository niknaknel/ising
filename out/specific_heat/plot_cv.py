import matplotlib.pyplot as plt
import pandas as pd
import sys

# read file
fname = sys.argv[1]
data = pd.read_csv(fname)
L = (fname.split('_')[1])[1:-4]

T = [0.05 * (i+1) for i in range(100)]
Tc = 2.2
alpha = 0.25
CV = [pow(t - Tc,-alpha) - pow(-Tc,-alpha) if not Tc-t == 0  else 1 for t in T]

# plot data
#plot = data.plot.line(x='T', y='Cv', style = '-o', title="Specific heat at different temperatures for L = %s" % L)
#fig = plot.get_figure()

#plt.plot(T, CV, '-')


fig = plt.figure()
plt.errorbar('T', 'mean_Cv', yerr='std_Cv', fmt='-o', data=data,
		ecolor='red', capsize=3, markersize=4)

plt.title('Specific heat at different temperatures for L = %s' % L)
plt.xlabel('Temperature')
plt.ylabel('Specific heat')
plt.ylim([-0.2, 2])
plt.grid()
plt.show()
fig.savefig(fname[:-4] + ".png")
