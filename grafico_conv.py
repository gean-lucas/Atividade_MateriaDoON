import numpy as np
import matplotlib.pyplot as plt

eq1=np.loadtxt('convergencia.txt');



plt.figure(figsize=(15,10))
plt.plot(eq1[:,0], eq1[:,1], 'o-', label='Menor rms por iteração')
plt.xlabel('iterações', fontsize=16)
plt.ylabel('rms', fontsize=16)
plt.legend(loc='best', numpoints=1, fontsize=16)
plt.show()

