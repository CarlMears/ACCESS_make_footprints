import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 18})
x = np.arange(0.0,1.0,0.001)
y = np.arange(0.0,1.0,0.001)

x_grid,y_grid = np.meshgrid(x,y)

r = np.sqrt(np.square(x_grid)+np.square(y_grid))

r = r*2.5

print(np.sum(r)/1000000.0)
r = np.reshape(r,(1000000))
h,edges =np.histogram(r,bins=50,range=[0,4])
fig,ax = plt.subplots()
x = np.arange(0,4.0,4.0/50.)
ax.plot(x,h/100000.)
ax.set_ylabel('Probability')
ax.set_xlabel('Distance Error (km)')
fig.tight_layout()
plt.show()


print()
