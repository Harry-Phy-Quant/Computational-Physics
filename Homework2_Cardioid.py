"""
三次样条函数插值绘制心形曲线
"""
import matplotlib.pyplot as plt
from Homework2_Advanced_Thomas import Advanced_Thomas
from math import *
from prettytable import PrettyTable


# 采用周期性边界条件的三次样条函数：输入三个参数，t_list, x_list 为一一对应的支撑点列表，t 为自变量的值
def CubicSpline_T(t_list, x_list, t):
    N = len(t_list) - 1  # 区间数
    h = t_list[1] - t_list[0]  # 两点间隔
    d = [0 for _ in range(N)]
    # 由边界条件确定 d[0]
    d[0] = 3 * (x_list[1] - x_list[0] + x_list[N-1] - x_list[N]) / (h ** 2)
    # 计算 d[1] ~ d[N-1]
    for i in range(1, N):
        d[i] = 3 * (x_list[i - 1] + x_list[i + 1] - 2 * x_list[i]) / (h ** 2)
    # 给出扩展的三对角矩阵
    mainDiagonal = [2 for _ in range(N)]
    upperDiagonal = [0.5 for _ in range(N-1)]
    lowerDiagonal = [0.5 for _ in range(N - 1)]
    # 用 Advanced_Thomas 算法解扩展的三对角矩阵得到 m 向量
    m = Advanced_Thomas(mainDiagonal, upperDiagonal, lowerDiagonal, 0.5, 0.5, d)
    # 添加上 m[N](=m[0]) 构成完整的 N+1 维的 m 向量
    m.append(m[0])
    # 用得到的 m 向量来计算任意给定 RK_x 值的预测值
    i = int((t - t_list[0]) / h)  # 确定 RK_x 所属的区间范围 RK_x[i] ~ RK_x[i+1]
    if i == N:
        return x_list[-1]
    else:
        A = ((x_list[i + 1] - x_list[i]) / h) - (h * (m[i + 1] - m[i]) / 6)
        B = x_list[i] - m[i] * (h ** 2) / 6
        S = - m[i] * ((t - t_list[i + 1]) ** 3) / (6 * h) + m[i + 1] * ((t - t_list[i]) ** 3) / (6 * h) + A * (
                t - t_list[i]) + B
        return S


support_t_list = [i for i in range(9)]
support_x_list = [(1-cos(pi*t/4))*cos(pi*t/4) for t in support_t_list]
support_y_list = [(1-cos(pi*t/4))*sin(pi*t/4) for t in support_t_list]

'''
# 方便将输出结果直接存入 Latex 表格
for i in range(len(support_t_list)):
    print(str(i+1)+"&"+str(support_t_list[i])+"&"+format(support_x_list[i], '0.6f')+"&"+format(support_y_list[i], '0.6f') + "\\" + "\\")
'''

# 真实值列表
t_list = [i / 2 for i in range(17)]
x_list = [(1-cos(pi*t/4))*cos(pi*t/4) for t in t_list]
y_list = [(1-cos(pi*t/4))*sin(pi*t/4) for t in t_list]

# 三次样条插值结果
CubicSpline_xlist = [CubicSpline_T(support_t_list, support_x_list, t) for t in t_list]
CubicSpline_ylist = [CubicSpline_T(support_t_list, support_y_list, t) for t in t_list]

# 作图
plt.plot(t_list, x_list, 'b.-', t_list, CubicSpline_xlist, 'r.:')
plt.title('Cubic Spline for RK_x(t)')
plt.legend(['  RK_x(t)',r'$S_{\Delta}(RK_x;t)$'],loc=4)
plt.show()

plt.plot(t_list, y_list, 'b.-', t_list, CubicSpline_ylist, 'r.:')
plt.title('Cubic Spline for RK_y(t)')
plt.legend(['  RK_y(t)',r'$S_{\Delta}(RK_y;t)$'],loc=4)
plt.show()

tb = PrettyTable()
tb.add_column('No.', [i for i in range(1, 18)])
tb.add_column('t', t_list)
tb.add_column('RK_x(t)', x_list)
tb.add_column('$S(RK_x;t)$', CubicSpline_xlist)
tb.add_column('RK_y(t)', y_list)
tb.add_column('$S(RK_y;t)$', CubicSpline_ylist)
print(tb)

'''
# 方便将输出结果直接存入 Latex 表格
for i in range(len(x_list)):
    print(str(i+1)+"&"+str(t_list[i])+"&"+format(x_list[i], '0.6f')+"&"+format(CubicSpline_xlist[i], '0.6f')+"&"+format(y_list[i], '0.6f')+"&"+format(CubicSpline_ylist[i], '0.6f')+"\\"+"\\")
'''

# 绘制二维心形曲线
h = 100
t = [i / h for i in range(8 * h + 1)]
X = [(1-cos(pi*t/4))*cos(pi*t/4) for t in t]
Y = [(1-cos(pi*t/4))*sin(pi*t/4) for t in t]
New_CubicSpline_xlist = [CubicSpline_T(support_t_list, support_x_list, t) for t in t]
New_CubicSpline_ylist = [CubicSpline_T(support_t_list, support_y_list, t) for t in t]

plt.plot(Y, X, 'b-', New_CubicSpline_ylist, New_CubicSpline_xlist, 'r--')
# 标注节点
for i in range(9):
    plt.plot((1-cos(pi*i/4))*sin(pi*i/4), (1-cos(pi*i/4))*cos(pi*i/4), 'ro')
plt.xlabel('RK_y(t) or $S_{\Delta}(RK_y;t)$')
plt.ylabel('RK_x(t) or $S_{\Delta}(RK_x;t)$')
plt.title('Cubic Spline for Cardioid Curve')
plt.legend(['Exact Curve',r'Cubic Spline Curve'], bbox_to_anchor=(0.82, 0.23))

plt.show()
