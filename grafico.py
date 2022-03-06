import numpy as np
import matplotlib.pyplot as plt

eq1=np.loadtxt('fit.txt');



plt.figure(figsize=(15,10))
plt.plot(eq1[:,0], eq1[:,1], 'o-', label='Anomalia calculada (nT)')
plt.plot(eq1[:,0], eq1[:,2], 'o-', label='Anomalia observada (nT)')
plt.xlabel('m', fontsize=16)
plt.ylabel('nT', fontsize=16)
plt.legend(loc='best', numpoints=1, fontsize=16)
plt.show()

