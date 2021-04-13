import numpy as np 
import matplotlib.pyplot as plt 

data = np.loadtxt('../build/test/ooa.txt')

x = data[:,0]
y = data[:,1]

fig, ax = plt.subplots()

ax.plot(x,y,'bo-',label='actual OOA')
ax.axhline(y=2.0, color = 'k', linestyle='dotted', label='nominal OOA')
ax.set_xlabel('$log(h_j)$')
ax.set_ylabel('$Δ log(E_j)/Δ log(h_j)$')
ax.set_title('Code verfication - actual order of accuracy')
ax.legend(loc='lower left')
plt.ylim(-1.0, plt.ylim()[1] + 0.5)
plt.grid(linestyle='dotted')
plt.savefig('../build/test/OVS.png', dpi=300)
