"""
迎风差分法求解一维对流方程
"""
import matplotlib.pyplot as plt

# 设置方程参数 a
a = -1
# 设置时间和空间范围
xini, xfinal = 0, -15
tini, tfinal = 0, 10
# 设置差分间隔，起始点和终止点与课件上略有不同
dx = 0.05
dt = 0.04
# dt = 0.05
# dt = 0.06
# 完整的段数
nx = int(abs(xfinal - xini) / dx)
nt = int(abs(tfinal - tini) / dt)
r = a * dt / dx
# 初始化，注意空间上考虑边值条件所以多加了一个节点
U = [[0.0 for _ in range(nt + 1)] for _ in range(nx + 2)]

# 初始条件：右侧注入的长度为10倍空间间隔的方波
for k in range(nx + 1):
    if k <= nx - 10:
        U[k][0] = 0
    else:
        U[k][0] = 1
# 边界条件：空间上最右侧的节点固定为0
for n in range(nt + 1):
    U[nx + 1][n] = 0

# 迭代差分求解
for n in range(nt):
    for k in range(nx + 1):
        U[k][n + 1] = (1 + r) * U[k][n] - r * U[k + 1][n]


# 作图
# 先画四张静态图
for t in range(nt + 1):
    if t == 0 or t == (nt + 1) // 3 or t == 2 * (nt + 1) // 3 or t == nt:
        x = [xfinal + i * dx for i in range(nx + 1)]
        plt.figure()
        plt.ylim(0, 1.2)  # 设置y轴坐标范围
        plt.xlim(xfinal, xini)  # 设置x轴坐标范围
        plt.title('Upwind solution for 1D Convective Equation')
        plt.grid(ls='-.')
        y = [U[x][t] for x in range(nx + 1)]
        plt.plot(x, y, linewidth='3', label='t = %.2f s' % (t * dt))
        plt.legend()
        plt.show()


# 动态图
def anni():
    plt.figure()
    plt.ion()  # 打开交互式绘图
    for t in range(nt + 1):
        plt.cla()  # 清除原有图像
        plt.ylim(0, 1.2)  # 设置y轴坐标范围
        plt.xlim(xfinal, xini)  # 设置x轴坐标范围
        plt.title('Upwind solution for 1D Convective Equation')
        plt.grid(ls='-.')
        x = [xfinal + i * dx for i in range(nx + 1)]
        y = [U[x][t] for x in range(nx + 1)]
        plt.plot(x, y, linewidth='2',  label='t = %.2f s' % (t * dt))
        plt.legend()
        plt.pause(0.1)
    plt.show()


if __name__ == '__main__':
    anni()
