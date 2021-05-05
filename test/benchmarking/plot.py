import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.family'] = 'STIXGeneral'

data = np.loadtxt('roof.txt')

# unpack data "visually" to serve as reference
# ridge1 = data[0]     ceil1 = data[1]     ridge2 = data[2]     ceil2 = data[3]
# ridge3 = data[4]     ceil3 = data[5]     
# opInt_zero = data[6] opInt_inft = data[7] 
# flops_zero = data[8] flops_infty = data[9]                    

# construct roofline ceilings
x1 = [0.0,     data[0], 100.0]
x2 = [data[0], data[2], 100.0]
x3 = [data[2], data[4], 100.0]
y1 = [0.0,     data[1], data[1]]
y2 = [data[1], data[3], data[3]]
y3 = [data[3], data[5], data[5]]
# construct vertical line for zero cache OI 
x4 = [data[6], data[6]] 
y4 = [0, data[1]]
# construct vertical line for infinite cache OI
x5 = [data[7], data[7]]
y5 = [0, data[1]]

fig, ax = plt.subplots()

# plot the data
ax.plot(x3, y3, color = 'darkgreen', linestyle = '-',
                                label= f'+TLP (max {data[5]} GFlops/s)')
ax.plot(x2, y2, color = 'darkred', linestyle = '-', 
                                label= f'+DLP (max {data[3]} GFlops/s)')
ax.plot(x1, y1, 'k-', label= f'C++ (max {data[1]} GFlops/s)')
ax.plot(x4, y4, color = 'k', linestyle = '--', 
                                label = f'Zero cache ({data[8]} GFlops/s)')
ax.plot(x5, y5, color = 'k', linestyle = ':', 
                                label = f'Infinite cache ({data[9]} GFlops/s)')
ax.scatter(data[6], data[8], s = 20, color = 'k')
ax.scatter(data[7], data[9], s = 20, color = 'k') 
ax.set_xlabel('Operational intensity (OI) [Flops/Byte]')
ax.set_ylabel('Performance [GFlops/s]')
ax.set_title(r'Roofline model for naive fourth-order CDS (double-precision)'
             '\n'
             r'Architecture: Intel Xeon E5-2670 v3'
             '\n'
             r'Compiler: GCC 10.1.0 with -O3 -g flags', loc='left')
ax.legend(loc='lower right')
ax.set_axisbelow(True)
plt.grid(axis = 'y', linestyle = ':', which = 'both')
plt.yscale('log')
plt.xscale('log')
plt.ylim(0.1, 200.0)
plt.xlim(0.01, 100.0)
plt.savefig('roof.png', dpi=600)
