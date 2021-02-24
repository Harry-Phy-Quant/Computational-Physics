"""
Householder 和 Givens 变换求 QR 分解用时比较
"""
from timeit import Timer
import matplotlib.pyplot as plt
from numpy import *
from prettytable import PrettyTable


Householder_Timer = Timer("Householder(A)",
                          "from Homework2_Householder import Householder; "
                          "from numpy import random; A = random.rand(6, 6)")
HouseHolder_result = Householder_Timer.repeat(20, 1)
Givens_Timer = Timer("Givens(A)",
                     "from Homework2_Givens import Givens; from numpy import random; A = random.rand(6, 6)")
Givens_result = Householder_Timer.repeat(20, 1)

tb = PrettyTable()
tb.add_column('No.', [i for i in range(1, 21)])
tb.add_column('HouseHolder_result(s)', HouseHolder_result)
tb.add_column('Givens_result(s)', Givens_result)
print(tb)

'''
# 方便将输出结果直接存入 Latex 表格
for i in range(len(HouseHolder_result)):
    print(str(i+1)+"&"+format(HouseHolder_result[i] * 1000, '0.4f')+"ms&"+format(Givens_result[i] * 1000, '0.4f')+"ms"+"\\"+"\\")


# 探究当矩阵维数更高时，两种算法的差别
HouseHolder_result = []
Givens_result = []
n_list = []
for n in range(1000, 1220, 20):
    n_list.append(n)
    Householder_Timer = Timer("Householder(A)",
                              "from Homework2_Householder import Householder; "
                              "from numpy import random; from __main__ import n; A = random.rand(n, n)")
    Givens_Timer = Timer("Givens(A)",
                         "from Homework2_Givens import Givens; "
                         "from numpy import random; from __main__ import n; A = random.rand(n, n)")
    HouseHolder_result.append(mean(Householder_Timer.repeat(5, 1)))
    Givens_result.append(mean(Givens_Timer.repeat(5, 1)))
plt.plot(n_list, HouseHolder_result, 'b.-', n_list, Givens_result, 'r.-')
plt.title('Running time for Householder and Givens')
plt.legend(['Householder', r'Givens'], loc=0)
plt.show()
'''
