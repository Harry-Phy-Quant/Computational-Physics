"""
求解点电荷的一维泊松方程
"""
import matplotlib.pyplot as plt
from Homework2_Thomas import Thomas

# 设置起点与终点
x_begin = 0
x_end = 1
# 设置间隔
N = 1000
dx = (x_end - x_begin) / (N + 1)
# 设置 δ 函数 v(x)
v = [0.0 for _ in range(N)]
for i in range(N):
    xx = x_begin + i * dx
    if abs(xx - 0.4) < 0.5 * dx:
        v[i] = 1 / dx
    else:
        v[i] = 0

# 设置三对角矩阵并求解
mainDiagonal = [-2 for _ in range(N)]
upperDiagonal = [1 for _ in range(N - 1)]
lowerDiagonal = [1 for _ in range(N - 1)]
d = [(dx ** 2) * v[i] for i in range(N)]
u = Thomas(mainDiagonal, upperDiagonal, lowerDiagonal, d)

# 加入 u_0, u_N+1
u.insert(0, 0)
u.append(0)

# 作图
x_list = [x_begin + i * dx for i in range(N + 2)]
plt.plot(x_list, u, 'r')
plt.title('1D Poisson')
plt.legend(['Numerical Method with N = %d' % N], loc=0)
plt.grid(ls='-.')
plt.show()
