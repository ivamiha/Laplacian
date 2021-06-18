import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['mathtext.fontset'] = 'cm'
plt.rcParams['font.family'] = 'STIXGeneral'

# benchmark results for single-precision 
sGOLD = 2.819
sGOLD_AUTO = 4.658
sGOLD_ISPC = 8.987
sMIDX = 0.933
sMIDX_AUTO = 0.938
s3DIDX = 3.840
s3DIDX_AUTO = 6.388
sSLICE = 3.628
sSLICE_AUTO = 5.478
sSLICE_ISPC = 6.944

# benchmark results for double-precision
dGOLD = 2.818
dGOLD_AUTO = 4.651
dGOLD_ISPC = 8.992
dMIDX = 0.944
dMIDX_AUTO = 0.955
d3DIDX = 3.827
d3DIDX_AUTO = 6.386
dSLICE = 3.629
dSLICE_AUTO = 5.478
dSLICE_ISPC = 6.961

# kernel specific values
sOI = 1.750
dOI = 0.875

# roofline characteristic values
maxFreq = 3.2
maxBand = 68.0
FMA = 4
sOI = 1.75
dOI = 0.875
ceil1 = maxFreq
sceil2 = maxFreq*FMA*4
dceil2 = maxFreq*FMA*2
ridge1 = ceil1 / maxBand
sridge2 = sceil2 / maxBand
dridge2 = dceil2 / maxBand

# construct roofline lines
x1 = [0.0, ridge1,  100.0]
x2 = [0.0, dridge2, 100.0]
x3 = [0.0, sridge2, 100.0]
x4 = [sOI, sOI];
x5 = [dOI, dOI];

y1 = [0.0, ceil1,  ceil1]
y2 = [0.0, dceil2, dceil2]
y3 = [0.0, sceil2, sceil2]
y4 = [0.0, sceil2]
y5 = [0.0, dceil2]

fig, ax = plt.subplots()

# plot roofline
ax.plot(x3, y3, color = 'k', linestyle = 'solid', zorder = 1,
                                    label = f'single-precision (max {sceil2} GFlops/s)')
ax.plot(x2, y2, color = 'dimgray', linestyle = 'solid', zorder = 1,
                                    label = f'double-precision (max {dceil2} GFlops/s)')
ax.plot(x1, y1, color = 'darkgray', linestyle = 'solid', zorder = 1,
                                    label = f'no parallelism (max {ceil1} GFlops/s')
ax.plot(x4, y4, color = 'k', linestyle = 'dashed', zorder = 1)
ax.plot(x5, y5, color = 'dimgray', linestyle = 'dashed', zorder = 1)

# plot data for single precision
ax.scatter(sOI, sGOLD, s = 20, color = 'gold', marker = 'o', zorder = 2)
ax.scatter(sOI, sGOLD_AUTO, s = 20, color = 'gold', marker = 's', zorder = 2)
ax.scatter(sOI, sGOLD_ISPC, s = 20, color = 'gold', marker = '^', zorder = 2)
ax.scatter(sOI, sMIDX, s = 20, color = 'darkgreen', marker = 'o', zorder = 2)
ax.scatter(sOI, sMIDX_AUTO, s = 20, color = 'darkgreen', marker = 's', zorder = 2)
ax.scatter(sOI, s3DIDX, s = 20, color = 'midnightblue', marker = 'o', zorder = 2)
ax.scatter(sOI, s3DIDX_AUTO, s = 20, color = 'midnightblue', marker = 's', zorder = 2)
ax.scatter(sOI, sSLICE, s = 20, color = 'darkred', marker = 'o', zorder = 2)
ax.scatter(sOI, sSLICE_AUTO, s = 20, color = 'darkred', marker = 's', zorder = 2)
ax.scatter(sOI, sSLICE_ISPC, s = 20, color = 'darkred', marker = '^', zorder = 2)

# plot data for double precision
ax.scatter(dOI, dGOLD, s = 20, color = 'gold', marker = 'o', zorder = 2)
ax.scatter(dOI, dGOLD_AUTO, s = 20, color = 'gold', marker = 's', zorder = 2)
ax.scatter(dOI, dGOLD_ISPC, s = 20, color = 'gold', marker = '^', zorder = 2)
ax.scatter(dOI, dMIDX, s = 20, color = 'darkgreen', marker = 'o', zorder = 2)
ax.scatter(dOI, dMIDX_AUTO, s = 20, color = 'darkgreen', marker = 's', zorder = 2)
ax.scatter(dOI, d3DIDX, s = 20, color = 'midnightblue', marker = 'o', zorder = 2)
ax.scatter(dOI, d3DIDX_AUTO, s = 20, color = 'midnightblue', marker = 's', zorder = 2)
ax.scatter(dOI, dSLICE, s = 20, color = 'darkred', marker = 'o', zorder = 2)
ax.scatter(dOI, dSLICE_AUTO, s = 20, color = 'darkred', marker = 's', zorder = 2)
ax.scatter(dOI, dSLICE_ISPC, s = 20, color = 'darkred', marker = '^', zorder = 2)

# plot configuration
ax.set_xlabel('Operational intensity [Flops/Byte]')
ax.set_ylabel('Performance [GFlops/s]')
ax.legend(loc='upper left')
ax.set_axisbelow(True)
plt.grid(axis = 'y', linestyle = ':', which = 'both')
plt.yscale('log')
plt.xscale('log')
plt.ylim(0.1, 350)
plt.xlim(0.01, 100.0)
plt.savefig('roof.png', dpi=600)
