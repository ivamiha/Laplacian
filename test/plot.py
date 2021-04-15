import numpy as np 
import matplotlib.pyplot as plt 

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.family'] = 'STIXGeneral'

data = np.loadtxt('../build/test/ooa.txt')

x = data[:,0]
y = data[:,1]

fig, ax = plt.subplots()

ax.plot(x,y,'ko-',label='actual OOA')
ax.axhline(y=4.0, color = 'k', linestyle='dotted', label='nominal OOA')
ax.set_xlabel('log($h_j$)')
ax.set_ylabel('Δlog($E_j$) / Δlog($h_j$)')
ax.set_title('Code verification for naive 4th-order flat indexing.')
ax.legend(loc='lower left')
plt.ylim(plt.ylim()[0]-1.0, plt.ylim()[1] + 0.5)
plt.savefig('../build/test/OVS.png', dpi=600)
