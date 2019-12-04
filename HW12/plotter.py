import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


plt.rcParams.update({'font.size': 22})
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams["legend.handlelength"] = 1.5


u = np.loadtxt('sol.txt')
ua = np.loadtxt('analytical.txt')
ucn = np.loadtxt('CN.txt')
ucn = np.transpose(ucn)
u_int=np.loadtxt('sym.txt')
xcn=np.loadtxt('xcn.txt')

t_array = np.loadtxt('time.txt')
t2_array=np.loadtxt('tcn.txt')
x = np.loadtxt('x_array.txt')


times= [0, 0.1, 0.2, 0.4]

tind=[]
tind2=[]
for i, t in enumerate(t_array):
    if t in times:
        tind.append(i)
for i, t, in enumerate(t2_array):
    if t in times:
        tind2.append(i)

c=['k','r','b', 'g', 'c', 'y']
ls = ['-', '--','-.',':']
lw=2.
ms=7.
fig, ax = plt.subplots(figsize=(8,6))

custom_lines=[Line2D([0], [0], color=c[0], lw=lw, ls=ls[0]),
            Line2D([0], [0], color=c[1], lw=lw, ls=ls[1]),
            Line2D([0], [0], color=c[2], lw=lw, ls=ls[2]),
            Line2D([0], [0], color=c[3], lw=lw, ls=ls[3]),
            Line2D([0], [0], color=c[0], lw=lw, ls="-"),
            Line2D([0], [0], color=c[0], lw=0, marker='o',markersize=ms, fillstyle='none')]
leg=[]
for i, t in enumerate(tind):
    ax.plot(x, ua[:,t], lw=lw, color=c[i], ls=ls[i])
    ax.plot(x, u[:,t], lw=0, marker='o',markersize=ms, fillstyle='none',
                    color=c[i], markevery=1)
    leg.append("t={} s".format(t_array[t]))

    ax.plot(x[:11],u_int[:,t], lw=lw, color='c')
# for i, t in enumerate(tind2):
#     ax.plot(xcn,ucn[:,t], color='c')
ax.set(xlabel='$x$', ylabel='$u$', xlim=[-1,1])
leg2=["Analytic", "Pseudospectral"]
leg1=ax.legend(custom_lines[:4], leg[:4], ncol=1, frameon=False, loc='upper left')
ax.legend(custom_lines[-2:], leg2, ncol=1, frameon=False, loc='lower right')
plt.gca().add_artist(leg1)

plt.tight_layout()
plt.savefig('figs/sol.eps')

#
# fig, ax = plt.subplots(figsize=(8,6))
#
# for i, t in enumerate(time):
#     ax.plot(x, u_an[t], lw=lw, label='t={} s'.format(t*dt), color=c[i], ls=ls[i])
#     ax.plot(x, u_cn[t], lw=0, marker='o',markersize=ms, fillstyle='none',
#                     label='t={} s'.format(t*dt), color=c[i], markevery=2)
# ax.set(xlabel='$x$', ylabel='$u$')
# leg=["t=0 s", "t=0.25 s", "t=0.5 s", "t=0.75 s", "Analytic", "Crank-Nicolson"]
# leg1=ax.legend(custom_lines[:4], leg[:4], ncol=1, frameon=False, loc='lower right')
# ax.legend(custom_lines[-2:], leg[-2:], ncol=1, frameon=False, loc='upper right')
# plt.gca().add_artist(leg1)
#
# plt.tight_layout()
# plt.savefig('figs/cn.eps')
# # ax.plot(a[0,25,'x'], a[0,25,'y'], lw=lw, label=r'a$_{{0}}$({1})'.format(0,5))
#

plt.show()
