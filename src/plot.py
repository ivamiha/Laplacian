import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.family'] = 'STIXGeneral'

data = np.loadtxt('../build/src/roof.txt')

# unpack data "visually" to serve as reference
# ridge1 = data[0]     ceil1 = data[1]     ridge2 = data[2]     ceil2 = data[3]
# ridge3 = data[4]     ceil3 = data[5]     opInty = data[6]     flops = data[7]

# construct roofline ceilings
x1 = [0.0,     data[0], 100.0]
x2 = [data[0], data[2], 100.0]
x3 = [data[2], data[4], 100.0]
y1 = [0.0,     data[1], data[1]]
y2 = [data[1], data[3], data[3]]
y3 = [data[3], data[5], data[5]]
# construct vertical line 
x4 = [data[6], data[6]] 
y4 = [0, data[1]]

fig, ax = plt.subplots()

# plot the data
ax.plot(x3, y3, color = 'darkgreen', linestyle = '-', label='+TLP')
ax.plot(x2, y2, color = 'darkred', linestyle = '-', label='+DLP')
ax.plot(x1, y1, 'k-', label='C++')
ax.plot(x4, y4, color = 'k', linestyle = '--')
ax.scatter(data[6], data[7], s = 50, color = 'k')
ax.set_xlabel('Operational intensity [FLOPs/Byte]')
ax.set_ylabel('Performance [GFLOPs/s]')
ax.set_title('Roofline model for naive 4th-order CDS')
ax.legend(loc='lower right')
plt.yscale('log')
plt.xscale('log')
plt.ylim(0.1, 200.0)
plt.xlim(0.01, 100.0)
plt.savefig('../build/src/roof.png', dpi=600)
