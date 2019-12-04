import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator


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

n = [4,8,12,16,20,30,40]

ps=[0.0248559669072,5.76778786747e-05,1.11772580711e-08,1.46707542945e-11,1.45854600688e-13,1.08375298231e-14, 1.00311371689e-14]
ps2 =[0.0248559669072, 5.76778786745e-05, 1.11771262695e-08, 8.22254435542e-13, 4.90947084158e-16, 3.04698745629e-15, 2.67953428207e-15]
cn=[0.0405955752319, 0.00723580158323,0.00295678133128,0.00160135186969,0.00100281719845, 0.000433401035878, 0.000240517930066]

lw = 3
ms = 8

c=['k','r','b', 'g', 'c', 'y']
ls = ['-', '--','-.',':']
fig, ax = plt.subplots(figsize=(8,6))

custom_lines=[Line2D([0], [0], color=c[0], lw=lw, ls=ls[0]),
            Line2D([0], [0], color=c[1], lw=lw, ls=ls[0]),
            Line2D([0], [0], color=c[2], lw=lw, ls=ls[0]),
            Line2D([0], [0], color=c[0], lw=lw, ls="-", marker='s', markersize=ms),
            Line2D([0], [0], color=c[0], lw=lw, marker='o',markersize=ms, fillstyle='none', ls=ls[1])]
# ax.plot(n, cn, lw=lw, marker='^', markersize=ms, color=c[1], ls=ls[1], label='Crank-Nicolson')

ax.plot(n, ps2, lw=lw, marker='s', markersize=ms, color=c[0], ls=ls[0], label='Pseudospectral (dopr5)')
# ax.plot(n, ps, lw=lw, marker='s', markersize=ms, color=c[3], ls=ls[3], label='Pseudospectral (dopr8)')

ax.set(xlabel='Number of grid points, $N$', ylabel='RMSE', yscale='log', xscale='linear')
# plt.xscale('log', basex=2)
# ax.legend(frameon=False)
# leg=["$N_x$=10", "$N_x$=20", "$N_x$=80", "FTCS", "C-N"]
# leg1=ax.legend(custom_lines[:3], leg[:3], ncol=3, frameon=False, loc='upper center')
# ax.legend(custom_lines[-2:], leg[-2:], ncol=2,loc='lower center',prop={'size': 16}, frameon=False, bbox_to_anchor= (0.5, 1.01))
# plt.gca().add_artist(leg1)


# plt.tick_params(axis='y', which='minor')

locmaj = matplotlib.ticker.LogLocator(base=10,numticks=8)
ax.yaxis.set_major_locator(locmaj)
# locmin = matplotlib.ticker.LogLocator(base=10.0,subs=(0.25,0.5),numticks=20)
# ax.yaxis.set_minor_locator(locmin)
# ax.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
# ax.xaxis.set_minor_locator(
#     AutoMinorLocator())
# ax.yaxis.set_minor_locator(
#     AutoMinorLocator())
ax.grid(which='both')

plt.tight_layout()
plt.savefig('figs/error.eps')




plt.show()
