"""
RK4 算法求解电磁学常微分方程组
"""
import matplotlib.pyplot as plt
from matplotlib.pyplot import MultipleLocator
from math import *

# 磁场方向为 𝑧 向： 𝜔 = 𝑞𝐵/m
w = 1.0
# 设置时间的起点与终点
t_begin = 0
t_end = 40
# 设置间隔
RK_step = 500
RK_h = (t_end - t_begin) / RK_step
Euler_step = 5000
Euler_h = (t_end - t_begin) / Euler_step
# RK4 法设置初值，由于 dvz/dt = 0，因此 RK_vz 为常数，直接赋值为初始值即可
RK_x, RK_y, RK_z, RK_vx, RK_vy, RK_vz = [0.0 for _ in range(RK_step + 1)], [0.0 for _ in range(RK_step + 1)], \
                                        [0.0 for _ in range(RK_step + 1)], [0.0 for _ in range(RK_step + 1)], \
                                        [0.0 for _ in range(RK_step + 1)], [0.1 for _ in range(RK_step + 1)]
RK_x[0], RK_y[0], RK_z[0] = 0.0, 0.0, 0.0
RK_vx[0], RK_vy[0] = 0.0, 2.0
# 向前差分法设置初值
Euler_x, Euler_y, Euler_z, Euler_vx, Euler_vy, Euler_vz = \
    [0.0 for _ in range(Euler_step + 1)], [0.0 for _ in range(Euler_step + 1)], \
    [0.0 for _ in range(Euler_step + 1)], [0.0 for _ in range(Euler_step + 1)], \
    [0.0 for _ in range(Euler_step + 1)], [0.1 for _ in range(Euler_step + 1)]
Euler_x[0], Euler_y[0], Euler_z[0] = 0.0, 0.0, 0.0
Euler_vx[0], Euler_vy[0] = 0.0, 2.0

# RK 法开始循环迭代求解
for k in range(RK_step):
    # 计算 RK_x, RK_y, RK_vx, RK_vy 的 k1~k4
    # RK_vx & RK_vy，由于 RK_vx, RK_vy 的 k1~k4 是耦合在一起的，因此顺序很重要
    k_vx_1 = w * RK_vy[k]
    k_vy_1 = -w * RK_vx[k]
    k_vx_2 = w * (RK_vy[k] + (RK_h / 2) * k_vy_1)
    k_vy_2 = -w * (RK_vx[k] + (RK_h / 2) * k_vx_1)
    k_vx_3 = w * (RK_vy[k] + (RK_h / 2) * k_vy_2)
    k_vy_3 = -w * (RK_vx[k] + (RK_h / 2) * k_vx_2)
    k_vx_4 = w * (RK_vy[k] + RK_h * k_vy_3)
    k_vy_4 = -w * (RK_vx[k] + RK_h * k_vx_3)

    # RK_x
    k_x_1 = RK_vx[k]
    k_x_2 = RK_vx[k] + (RK_h / 2) * k_vx_1
    k_x_3 = RK_vx[k] + (RK_h / 2) * k_vx_2
    k_x_4 = RK_vx[k] + RK_h * k_vx_3

    # RK_y
    k_y_1 = RK_vy[k]
    k_y_2 = RK_vy[k] + (RK_h / 2) * k_vy_1
    k_y_3 = RK_vy[k] + (RK_h / 2) * k_vy_2
    k_y_4 = RK_vy[k] + RK_h * k_vy_3

    # 下面计算所有自变量下标为 k+1 的值
    RK_x[k + 1] = RK_x[k] + (RK_h / 6) * (k_x_1 + 2 * k_x_2 + 2 * k_x_3 + k_x_4)
    RK_y[k + 1] = RK_y[k] + (RK_h / 6) * (k_y_1 + 2 * k_y_2 + 2 * k_y_3 + k_y_4)
    RK_z[k + 1] = RK_z[k] + RK_h * RK_vz[k]
    RK_vx[k + 1] = RK_vx[k] + (RK_h / 6) * (k_vx_1 + 2 * k_vx_2 + 2 * k_vx_3 + k_vx_4)
    RK_vy[k + 1] = RK_vy[k] + (RK_h / 6) * (k_vy_1 + 2 * k_vy_2 + 2 * k_vy_3 + k_vy_4)

# 向前差分法开始循环迭代求解
for k in range(Euler_step):
    Euler_x[k + 1] = Euler_x[k] + Euler_h * Euler_vx[k]
    Euler_y[k + 1] = Euler_y[k] + Euler_h * Euler_vy[k]
    Euler_z[k + 1] = Euler_z[k] + Euler_h * Euler_vz[k]
    Euler_vx[k + 1] = Euler_vx[k] + w * Euler_h * Euler_vy[k]
    Euler_vy[k + 1] = Euler_vy[k] - w * Euler_h * Euler_vx[k]

# 计算精确解
exact_nstep = 10000
exact_x, exact_y, exact_z, exact_vx, exact_vy, exact_vz = \
    [0.0 for _ in range(exact_nstep+1)], [0.0 for _ in range(exact_nstep+1)], [0.0 for _ in range(exact_nstep+1)], \
    [0.0 for _ in range(exact_nstep+1)], [0.0 for _ in range(exact_nstep+1)], [0.1 for _ in range(exact_nstep+1)]
for k in range(exact_nstep + 1):
    exact_x[k] = (2 / w) * (1 - cos(w * k * (t_end - t_begin) / exact_nstep))
    exact_y[k] = (2 / w) * sin(w * k * (t_end - t_begin) / exact_nstep)
    exact_z[k] = exact_vz[0] * k * (t_end - t_begin) / exact_nstep
    exact_vx[k] = 2 * sin(w * k * (t_end - t_begin) / exact_nstep)
    exact_vy[k] = 2 * cos(w * k * (t_end - t_begin) / exact_nstep)

# 作图
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(RK_x, RK_y, RK_z, color='r', lw=2, label='RK4,step=%d' % RK_step)
ax.plot(Euler_x, Euler_y, Euler_z, color='g', lw=2, label='Euler Explicit Method,step=%d' % Euler_step)
ax.plot(exact_x, exact_y, exact_z, color='b', lw=1, label='Exact Solution')
ax.set_xlabel("RK_x")
ax.set_ylabel("RK_y")
ax.set_zlabel("RK_z")
ax.legend()
locator = MultipleLocator(1)
ax.xaxis.set_major_locator(locator)
ax.yaxis.set_major_locator(locator)
ax.zaxis.set_major_locator(locator)
plt.show()
