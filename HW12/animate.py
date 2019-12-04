import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from matplotlib.lines import Line2D
from matplotlib.animation import FuncAnimation
import matplotlib.animation as animation
import sys

plt.rcParams.update({'font.size': 22})
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams["legend.handlelength"] = 1.5

plt.rcParams['animation.convert_path']= "/usr/local/Cellar/imagemagick/7.0.9-7/bin/magick"


u = np.loadtxt('sol.txt')
ua = np.loadtxt('analytical.txt')

t_array = np.loadtxt('time.txt')
x = np.loadtxt('x_array.txt')


c=['k','r','b', 'g', 'c', 'y']
ls = ['-', '--','-.',':']
lw=2.
ms=7.

custom_lines=[Line2D([0], [0], color=c[0], lw=lw, ls=ls[0]),
            # Line2D([0], [0], color=c[1], lw=lw, ls=ls[1]),
            # Line2D([0], [0], color=c[2], lw=lw, ls=ls[2]),
            # Line2D([0], [0], color=c[3], lw=lw, ls=ls[3]),
            # Line2D([0], [0], color=c[0], lw=lw, ls="-"),
            Line2D([0], [0], color=c[0], lw=0, marker='o',markersize=ms, fillstyle='none')]

fig, ax = plt.subplots(figsize=(8,6))
fig.set_tight_layout(True)

print('fig size: {0} DPI, size in inches {1}'.format(
    fig.get_dpi(), fig.get_size_inches()))

an, = ax.plot(x, u[:,0], lw=lw, color=c[0], ls=ls[0])
sol, = ax.plot(x, ua[:,0], lw=0, marker='o',markersize=ms, fillstyle='none',color=c[0])
ax.set(xlabel='$x$', ylabel='$u$', xlim=[-1,1])

leg=["Analytic", "Pseudospectral"]
# leg1=["t=0 s"]
# leg1=ax.legend(custom_lines, leg1, ncol=1, frameon=False, loc='upper left')
ax.legend(custom_lines[-2:], leg, ncol=1, frameon=False, loc='lower right')
# plt.gca().add_artist(leg1)
ax.set_title('t=0 s')


def update(i):
    label = 'timestep {0}'.format(i)
    print(label)
    # Update the line and the axes (with a new xlabel). Return a tuple of
    # "artists" that have to be redrawn for this frame.
    an.set_ydata(ua[:,i])
    sol.set_ydata(u[:,i])
    ax.set_title('t = {0:.2f} s'.format(t_array[i]))

    return an,sol, ax


if __name__ == '__main__':
    # FuncAnimation will call the 'update' function for each frame; here
    # animating over 10 frames, with an interval of 200ms between frames.
    anim = FuncAnimation(fig, update, frames=np.arange(0, 100), interval=50)
    anim.save('test.gif', dpi=100, writer='imagemagick')
    # else:
    #     # plt.show() will just loop the animation forever.
    #     plt.show()
# leg=[]
# for i, t in enumerate(tind):
#     ax.plot(x, ua[:,t], lw=lw, color=c[i], ls=ls[i])
#     ax.plot(x, u[:,t], lw=0, marker='o',markersize=ms, fillstyle='none',color=c[i])
#     leg.append("t={} s".format(t_array[t]))
# for i, t in enumerate(tind2):
#     ax.plot(xcn,ucn[:,t], color='c')



# plt.savefig('figs/sol.eps')

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
