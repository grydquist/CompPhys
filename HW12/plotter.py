import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


plt.rcParams.update({'font.size': 22})
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams["legend.handlelength"] = 1.5


u = np.loadtxt('u_plot.txt', delimiter=',')
ucn=np.loadtxt('CN_plot.txt', delimiter=',')
# u = np.loadtxt('u.txt', delimiter=',')
# ucn=np.loadtxt('CN.txt', delimiter=',')
L=10.
dt=1.e-3
# dt=0.0075
t_end = 1.

n = int(t_end/dt)
nx=100

dx = L/(nx-1)

x=[]
xn=[]
for i in range(nx):
    x.append(i*dx)
    xn.append(-i*dx)

xn = np.flip(xn)

u_ftcs = []
u_an = []
u_cn =[]
t = []
for i in range(n):
    t.append(dt*i)
    u_ftcs.append(u[(i)*nx:(i+1)*nx, 0])
    u_an.append(ucn[(i)*nx:(i+1)*nx, 1])
    u_cn.append(ucn[(i)*nx:(i+1)*nx, 0])

lw = 2
ms = 6

time = [0, int(n/4), int(n/2), int(3*n/4)]
c=['k','r','b', 'g', 'c', 'y']
ls = ['-', '--','-.',':']
fig, ax = plt.subplots(figsize=(8,6))

custom_lines=[Line2D([0], [0], color=c[0], lw=lw, ls=ls[0]),
            Line2D([0], [0], color=c[1], lw=lw, ls=ls[1]),
            Line2D([0], [0], color=c[2], lw=lw, ls=ls[2]),
            Line2D([0], [0], color=c[3], lw=lw, ls=ls[3]),
            Line2D([0], [0], color=c[0], lw=lw, ls="-"),
            Line2D([0], [0], color=c[0], lw=0, marker='o',markersize=ms, fillstyle='none'),
            Line2D([0], [0], color=c[0], lw=0, marker='x',markersize=ms*1.5)]
for i, t in enumerate(time):
    ax.plot(x, u_an[t], lw=lw, label='t={} s'.format(t*dt), color=c[i], ls=ls[i])
    ax.plot(x, u_ftcs[t], lw=0, marker='o',markersize=ms, fillstyle='none',
                    label='t={} s'.format(t*dt), color=c[i], markevery=3)
    ax.plot(x, u_cn[t], lw=0, marker='x',markersize=ms*1.5, fillstyle='none',
                    label='t={} s'.format(t*dt), color=c[i], markevery=3)

    ax.plot(xn, np.flip(u_an[t]), lw=lw, label='t={} s'.format(t*dt), color=c[i], ls=ls[i])
    ax.plot(xn, np.flip(u_ftcs[t]), lw=0, marker='o',markersize=ms, fillstyle='none',
                    label='t={} s'.format(t*dt), color=c[i], markevery=3)
    ax.plot(xn, np.flip(u_cn[t]), lw=0, marker='x',markersize=ms*1.5, fillstyle='none',
                    label='t={} s'.format(t*dt), color=c[i], markevery=3)
ax.set(xlabel='$x$', ylabel='$u$', xlim=[-5,5])
leg=["t=0 s", "t=0.25 s", "t=0.5 s", "t=0.75 s", "Analytic", "FTCS", "CN"]
leg1=ax.legend(custom_lines[:4], leg[:4], ncol=1, frameon=False, loc='upper left')
ax.legend(custom_lines[-3:], leg[-3:], ncol=1, frameon=False, loc='upper right')
plt.gca().add_artist(leg1)

plt.tight_layout()
plt.savefig('figs/all.eps')

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
