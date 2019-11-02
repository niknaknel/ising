import matplotlib.pyplot as plt
import pandas as pd
import glob

# read file
files = glob.glob("*.csv")
names = [["temperature", "t%d" % i] for i in range(1, len(files)+1)]
print(files)

df = pd.concat([pd.read_csv(files[i], names=names[i]) for i in range(len(files))], axis=1, sort=False)
df = df.iloc[:,~df.columns.duplicated()]
df['ave'] = df.iloc[:, 1:len(files)+1].mean(axis=1)
# print(df)

# plot data
fig = plt.figure()
ax = plt.axes()
ax.set_xlabel("Temperature (K)")
ax.set_ylabel("Equilibrium time (steps)")
ax.set_title("Phase diagram for L=100")
ax.plot(df['temperature'], df['ave'])
plt.show()
