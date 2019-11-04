import matplotlib.pyplot as plt
import pandas as pd
import sys

# read file
fname = sys.argv[1]
data = pd.read_csv(fname)

T = [0.05 * (i+1) for i in range(100)]
Tc = 2.2
alpha = 0.25
CV = [pow(t - Tc,-alpha) - pow(-Tc,-alpha) if not Tc-t == 0  else 1 for t in T]

# plot data
#plot = data.plot.line(x='T', y='tau', style = '-o', title="Autocorrelation at different temperatures.")
plot = data.plot.line(x='T', y='Cv', style = '-o', title="Specific heat at different temperatures.")
fig = plot.get_figure()


#plt.plot(T, CV, '-')

plt.xlabel('Temperature')
#plt.ylabel('Autocorrelation time')
plt.ylabel('Specific heat')
#plt.ylim([-0.2, 1.2])
plt.grid()
plt.show()
fig.savefig(fname[:-4] + ".png")
