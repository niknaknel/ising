import matplotlib.pyplot as plt
import pandas as pd
import sys

# read file
fname = sys.argv[1]
data = pd.read_csv(fname, names=['X(t)'])

if "corr" in fname:
    # plot data
    vals = fname.split('_')
    T = (vals[1])[1:]
    L = (vals[2])[1:]
    tau = int((vals[3])[3:-4])

    title = "Autocorrelation at T=%s with L=%s" % (T, L)

    plot = data.plot(kind = 'line', title=title)
    fig = plot.get_figure()

    plt.xlabel('Time steps')
    plt.ylabel('X(t)/X(0)')
    plt.ylim([-0.2, 1])
    plt.axvline(x=tau, color='red')
    plt.grid()

    #plt.show()
    fig.savefig(fname[:-4] + ".png")

else:
    plot = data.plot(x='T', y='tau', style = '-o', title="Correlation time vs Temperature")
    fig = plot.get_figure()

    plt.xlabel('Temperature')
    plt.ylabel('Correlation time')
    plt.ylim([0, 180])
    plt.grid()

    #plt.show()
    fig.savefig(fname[:-4] + ".png")
