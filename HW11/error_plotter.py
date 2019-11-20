import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


plt.rcParams.update({'font.size': 22})
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams["legend.handlelength"] = 1.5
plt.rcParams["legend.columnspacing"] = 0.75

L=10.
dt=1.e-3
t_end = 1.

n = int(t_end/dt)
nxf = [10, 20, 80]

dtf = [1e-1,1e-2, 1e-3]
ftcs = [[3.17680264142e-05, 6.32190002819e-05,6.6178208276e-05],
        [5.80308699103e-06, 2.28072088398e-06, 3.43149371069e-06],
        # [0.000541309688957, 2.32316805015e-09, 2.19041674699e-07],
        [0.251889616814, 0.0840887692146, 6.89042177148e-09]]

cnf = [[6.1517728491e-05, 6.6599378832e-05, 6.65179025045e-05],
        [2.8141503822e-06,3.57383197946e-06,  3.57280360729e-06],
        # [9.15383812848e-08, 2.62015052681e-07, 2.64271100485e-07],
        [7.44443020801e-08, 1.67645224284e-08, 1.75290127582e-08]]


nx = [250, 500, 1000]
dt = [0.1, 0.05, 0.025,0.0125, 0.0050625]
cn = [[1.43199547414e-07,7.65537418218e-09,2.19833679799e-10,8.59403587016e-11,1.68393808397e-10],
        [1.52707214609e-07,9.3951636833e-09, 4.7783334906e-10,1.37482530063e-11,7.12937447741e-12],
        [1.55872331268e-07, 9.91847314621e-09, 5.80768656559e-10, 2.9626785174e-11, 2.94207867541e-13]]


lw = 2
ms = 6

c=['k','r','b', 'g', 'c', 'y']
ls = ['-', '--','-.',':']
fig, ax = plt.subplots(figsize=(8,6))

custom_lines=[Line2D([0], [0], color=c[0], lw=lw, ls=ls[0]),
            Line2D([0], [0], color=c[1], lw=lw, ls=ls[0]),
            Line2D([0], [0], color=c[2], lw=lw, ls=ls[0]),
            Line2D([0], [0], color=c[0], lw=lw, ls="-", marker='s', markersize=ms),
            Line2D([0], [0], color=c[0], lw=lw, marker='o',markersize=ms, fillstyle='none', ls=ls[1])]
for i, n in enumerate(nxf):
    ax.plot(dtf, ftcs[i], lw=lw, marker='s', markersize=ms,label='nx={}'.format(nxf[i]), color=c[i], ls=ls[0])
    ax.plot(dtf, cnf[i], lw=lw, marker='o',markersize=ms, fillstyle='none', color=c[i], ls=ls[1])
ax.set(xlabel='$dt$', ylabel='RMSE', yscale='log', xscale='log', ylim=[1e-9, 1e1])
leg=["nx=10", "nx=20", "nx=80", "FTCS", "C-N"]
leg1=ax.legend(custom_lines[:3], leg[:3], ncol=3, frameon=False, loc='upper center')
ax.legend(custom_lines[-2:], leg[-2:], ncol=2,loc='lower center',prop={'size': 16}, frameon=False, bbox_to_anchor= (0.5, 1.01))
plt.gca().add_artist(leg1)

plt.tight_layout()
plt.savefig('figs/error.eps')

fig, ax = plt.subplots(figsize=(8,6))

custom_lines=[Line2D([0], [0], color=c[0], lw=lw, ls=ls[0]),
            Line2D([0], [0], color=c[1], lw=lw, ls=ls[1]),
            Line2D([0], [0], color=c[2], lw=lw, ls=ls[2]),
            Line2D([0], [0], color=c[0], lw=lw, ls="-", marker='s', markersize=ms),
            Line2D([0], [0], color=c[0], lw=lw, marker='o',markersize=ms, fillstyle='none', ls=ls[0])]
for i, n in enumerate(nx):
    ax.plot(dt, cn[i], lw=lw, marker='o', markersize=ms,fillstyle='none',label='nx={}'.format(nxf[i]), color=c[i], ls=ls[i])
ax.set(xlabel='$dt$', ylabel='RMSE', yscale='log', xscale='log', ylim=[1e-13, 1e-5])
leg=["nx=250", "nx=500", "nx=1000", "C-N"]
leg1=ax.legend(custom_lines[:3], leg[:3], ncol=3, frameon=False, loc='upper center')
ax.legend(custom_lines[-1:], leg[-1:], ncol=2,loc='lower center',prop={'size': 16}, frameon=False, bbox_to_anchor= (0.5, 1.01))
plt.gca().add_artist(leg1)

plt.tight_layout()
plt.savefig('figs/error2.eps')




plt.show()
