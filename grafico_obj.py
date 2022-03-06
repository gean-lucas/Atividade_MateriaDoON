import numpy as np
import matplotlib.pyplot as plt

eq1=np.loadtxt('grid.txt')
eq2=np.loadtxt('phi.txt')

x_grid, y_grid, phi_grid = eq1[:,0], eq1[:,1], eq1[:,2]
x, y, phi = eq2[:,0], eq2[:,1], eq2[:,2]

xmax = np.max(x_grid)
xmin = np.min(y_grid)
    
ymax = np.max(y_grid)
ymin = np.min(y_grid)
    
plt.figure(figsize=(8,6))
plt.axis('scaled')
plt.contourf(y_grid,x_grid,phi_grid,40, cmap = plt.get_cmap('viridis'))
cbar = plt.colorbar()
cbar.set_label('obj function', fontsize=16)
plt.xlabel('xq (m)', fontsize=14)
plt.ylabel('zq (m)', fontsize=14)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.xlim(ymin,ymax)
plt.ylim(xmin,xmax)
plt.show()
