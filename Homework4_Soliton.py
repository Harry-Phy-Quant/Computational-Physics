"""
差分法求解 Kruskal-Zalusky 孤立子 （笔记本运行时间约为两分半）
"""
import matplotlib.pyplot as plt
from math import *
import matplotlib.animation as animation

# 设置方程参数 a
delta = 0.022
# 设置时间和空间范围
xini, xfinal = 0.0, 2.0
tini, tfinal = 0.0, 4.0
# 设置差分间隔 （N为空间，M为时间）
N = 250
M = 200000
dx = (xfinal - xini) / N
dt = (tfinal - tini) / M
# 初始化，注意空间上考虑边值条件所以多加了几个节点
U = [[0.0 for _ in range(M + 1)] for _ in range(N + 4)]
for i in range(N + 1):
    U[i][0] = cos(pi * (xini + i * dx))
U[N + 1][0] = U[1][0]
U[N + 2][0] = U[2][0]
U[N + 3][0] = U[2][0]

# 第一次迭代差分求解
for i in range(2, N + 2):
    U[i][1] = U[i][0] - (1 / 6) * (dt / dx) * (U[i + 1][0]+U[i][0]+U[i - 1][0]) * (U[i + 1][0] - U[i - 1][0]) - \
              ((delta ** 2) * dt / (2 * (dx ** 3))) * (U[i + 2][0] - 2 * U[i + 1][0] + 2 * U[i - 1][0] - U[i - 2][0])
# 利用周期性边界条件
U[0][1] = U[N][1]
U[1][1] = U[N + 1][1]
U[N + 2][1] = U[2][1]
U[N + 3][1] = U[3][1]

# 第二次及以后迭代差分求解
for j in range(2, M + 1):
    for i in range(2, N + 2):
        U[i][j] = U[i][j - 2] - (1 / 3) * (dt / dx)\
                  * (U[i + 1][j - 1] + U[i][j - 1] + U[i - 1][j - 1]) * (U[i + 1][j - 1] - U[i - 1][j - 1])\
                  - ((delta ** 2) * dt / (dx ** 3)) *\
                  (U[i + 2][j - 1] - 2 * U[i + 1][j - 1] + 2 * U[i - 1][j - 1] - U[i - 2][j - 1])
    # 利用周期性边界条件
    U[0][j] = U[N][j]
    U[1][j] = U[N + 1][j]
    U[N + 2][j] = U[2][j]
    U[N + 3][j] = U[3][j]


# 动态图
fig, ax = plt.subplots()
x = [xini + k * dx for k in range(N + 1)]
t = [tini + k * dt for k in range(M + 1)]
y = [U[x][0] for x in range(N + 1)]
line, = ax.plot(x, y)
text_pt = plt.text(1.6, 2.7, '', fontsize=18)


def animate(m):
    ax.set_title('Kruskal-Zalusky Soliton')
    ax.grid(ls='-.')
    ax.set_xlim(xini, xfinal)  # 设置x轴坐标范围
    ax.set_ylim(-1.5, 3)  # 设置y轴坐标范围
    line.set_ydata([U[k][100*m] for k in range(N + 1)])
    text_pt.set_text("t=%.3f" % t[100*m])
    # 用来保存静态图
    # if m == 0 or m == 159 or m == 573 or m == 1114 or m == 1556 or m == 1999:
    #    plt.pause(2)
    return line,


ani = animation.FuncAnimation(fig, animate, frames=int((M/100)), interval=0.01)
plt.show()
