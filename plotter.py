import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

p = pd.read_csv("out.txt", dtype=np.float64)

print(p)


print(type(p['H'].loc[0]))

p['H'].plot()
plt.show()

plt.plot(p['X'], p['Y'])
plt.show()

#p.dropna().iloc[1000:].plot.hist(bins=5000, density=True)
#plt.title(p['H'].dropna().iloc[1000:].mean())
#plt.show()