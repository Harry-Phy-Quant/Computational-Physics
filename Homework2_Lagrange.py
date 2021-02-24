"""
Lagrange 内插
"""
from prettytable import PrettyTable
import matplotlib.pyplot as plt


# 输入三个参数，x_list, y_list 为一一对应的支撑点列表，RK_x 为自变量的值
def Lagrange(x_list, y_list, x):
    # 阶数为 N
    N = len(x_list) - 1
    P = 0
    for i in range(N + 1):
        L = 1
        for m in range(N + 1):
            if m != i:
                L *= (x - x_list[m]) / (x_list[i] - x_list[m])  # 构建 N 次多项式
        P += y_list[i] * L  # 构建 Lagrange 多项式
    return P


if __name__ == '__main__':
    # 以 Runge 现象为例
    def Runge(x):
        return 1 / (1 + 25 * (x**2))

    support_x_list = [-1 + 0.1 * i for i in range(21)]
    support_y_list = [Runge(i) for i in support_x_list]
    x_list = [-1 + 0.05 * i for i in range(41)]
    y_list = [Runge(i) for i in x_list]
    Lagrange_list = [Lagrange(support_x_list, support_y_list, i) for i in x_list]
    difference = [abs(Lagrange_list[i] - y_list[i]) for i in range(len(y_list))]

    tb = PrettyTable()
    tb.add_column('No.', [i for i in range(1, 42)])
    tb.add_column('RK_x', x_list)
    tb.add_column('f(RK_x)', y_list)
    tb.add_column('P20(RK_x)', Lagrange_list)
    tb.add_column('|P20(RK_x) - f(RK_x)|', difference)
    print(tb)

    plt.plot(x_list, y_list, 'b.-', x_list, Lagrange_list, 'r.:')
    plt.ylim((-0.5, 1.5))
    plt.title('Lagrange Interpolation')
    plt.legend(['  f(RK_x)', r'$P_{20}(RK_x)$'], loc=1)
    plt.show()

    ''' 方便将输出结果直接存入 Latex 表格
    for i in range(len(x_list)):
        print(str(i+1)+"&"+format(x_list[i], '0.4f')+"&"+format(y_list[i], '0.4f')+"&"+format(Lagrange_list[i], '0.4f')+"&"+format(difference[i], '0.4f')+"\\"+"\\")
    '''
