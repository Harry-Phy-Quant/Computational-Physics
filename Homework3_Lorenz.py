"""
向前差分法求解洛伦兹吸引子
"""
import matplotlib.pyplot as plt

# 设置起点与终点
x_begin = 0
x_end = 10
# 设置间隔
nstep = 100000
h = (x_end - x_begin) / nstep
# 设置参数
beta = 8 / 3
rho = 28.0
sigma = 10.0
# 设置初值
x, y, z = [0.0 for _ in range(nstep + 1)], [0.0 for _ in range(nstep + 1)], [0.0 for _ in range(nstep + 1)]
x[0], y[0], z[0] = 12.0, 4.0, 0.0

# 向前差分法开始循环迭代求解
for k in range(nstep):
    x[k + 1] = x[k] + h * (-beta * x[k] + y[k] * z[k])
    y[k + 1] = y[k] + h * (-sigma * y[k] + sigma * z[k])
    z[k + 1] = z[k] + h * (-y[k] * x[k] + rho * y[k] - z[k])

# 作图
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot(x, y, z, color='r', lw=2, label='beta=%.3f, rho=%.3f, sigma=%.3f' % (beta, rho, sigma))
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
ax.set_title("Lorenz Attractor")
ax.legend()
plt.show()
