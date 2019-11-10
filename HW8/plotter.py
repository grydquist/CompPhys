import matplotlib.pyplot as plt
import numpy as np
from basic_units import radians, degrees, cos

def format_func(value, tick_number):
    # find number of multiples of pi/2
    N = int(np.round(2 * value / np.pi))
    if N == 0:
        return "0"
    elif N == 1:
        return r"$\pi/2$"
    elif N == 2:
        return r"$\pi$"
    elif N % 2 > 0:
        return r"${0}\pi/2$".format(N)
    else:
        return r"${0}\pi$".format(N // 2)


plt.rcParams.update({'font.size': 20})
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams["legend.handlelength"] = 1.0

ns = [0,1,2,10,15]
qs = [5,25]

a = {}
b = {}
for q in qs:
    for n in ns:
        data = np.loadtxt('an{0}q{1}.txt'.format(n, q))
        a[(n,q,'x')] = data[:,0]
        a[(n,q,'y')] = data[:,1]
        if n != 0:
            data = np.loadtxt('bn{0}q{1}.txt'.format(str(n), q))
            b[(n,q,'x')] = data[:,0]
            b[(n,q,'y')] = data[:,1]


lw = 2

for q in qs:
    for n in ns:
        if n==0:
            fig, ax = plt.subplots(figsize=(8,4))
            ax.plot(a[n,q,'x'], a[n,q,'y'], lw=lw, label='$y_\mathrm{I}$'+' n={0}, q={1}'.format(n,q), color='r')
            ax.set(xlabel='$t$', ylabel=r'Eigenfunction $y$')
            ax.legend(ncol=5, loc='lower center',prop={'size': 16}, frameon=False, bbox_to_anchor= (0.5, 1.01))
            ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
            plt.tight_layout()
            plt.savefig('figs/a_{}_{}.eps'.format(n,q))

        else:
            fig, ax = plt.subplots(figsize=(8,4))
            ax.plot(a[n,q,'x'], a[n,q,'y'], lw=lw, label='$y_{\mathrm{I}}$'+' n={0}, q={1}'.format(n,q), color='r')
            ax.plot(b[n,q,'x'], b[n,q,'y'], lw=lw, label='$y_{\mathrm{II}}$'+' n={0}, q={1}'.format(n,q), color='b', ls='--')
            ax.set(xlabel='$t$', ylabel=r'Eigenfunction $y$')
            ax.legend(ncol=5, loc='lower center',prop={'size': 16}, frameon=False, bbox_to_anchor= (0.5, 1.01))
            ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
            plt.tight_layout()
            plt.savefig('figs/{}_{}.eps'.format(n,q))
        # ax.plot(a[0,25,'x'], a[0,25,'y'], lw=lw, label=r'a$_{{0}}$({1})'.format(0,5))



# q=5
# for n in ns:
#     ax.plot(a[(n,q,'x')], a[n,q,'y'], lw=lw, label=r'a$_{{0}}$({1})'.format(n,q))
#     if n != 0:
#         ax.plot(b[(n,q,'x')], b[n,q,'y'], lw=lw, label=r'b$_{{0}}$({1})'.format(n,q))
#
# ax.set(xlabel='$t$', ylabel='Eigenfunction')
# ax.legend(ncol=5, loc='upper center',prop={'size': 16}, frameon=False)
#
# ax=axs[1]
# q=25
# for n in ns:
#     ax.plot(a[(n,q,'x')], a[n,q,'y'], lw=lw, label=r'a$_{{0}}$({1})'.format(n,q))
#
#     if n != 0:
#         ax.plot(b[(n,q,'x')], b[n,q,'y'], lw=lw, label=r'b$_{{0}}$({1})'.format(n,q))
# ax.set(xlabel='$t$', ylabel='Eigenfunction')
# ax.legend(ncol=5, loc='upper center',prop={'size': 16}, frameon=False)



plt.show()
