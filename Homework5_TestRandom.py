"""
np.random.uniform 的均匀性卡方检验
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# 每次试验撒点数
N = 10000
# 进行 NN 次试验
NN = 500
# NN = 5000
# NN = 10000
Frequency = []
for i in range(NN):
    # 每次试验随机撒点
    numbers = np.random.uniform(0, 1, N)
    # 统计各个区间的撒点数
    frequency = [0 for _ in range(10)]
    for number in numbers:
        for k in range(10):
            if k / 10 <= number < (k + 1) / 10:
                frequency[k] += 1
    Frequency.append(frequency)

    ''' # 方便将输出结果直接存入 Latex 表格
    print(str(Frequency[i][0]) + "&" + str(Frequency[i][1]) + "&" + str(Frequency[i][2]) + "&" + str(
        Frequency[i][3]) + "&" +
          str(Frequency[i][4]) + "&" + str(Frequency[i][5]) + "&" + str(Frequency[i][6]) + "&" + str(
        Frequency[i][7]) + "&" +
          str(Frequency[i][8]) + "&" + str(Frequency[i][9]) + "\\" + "\\")
    '''

# 计算卡方统计量
Chi2 = []
mk = 0.1 * N
for frequency in Frequency:
    chi2 = 0
    for n in frequency:
        chi2 += (n - mk) ** 2 / mk
    # print("第" + str(Frequency.index(frequency)+1) + "次试验的卡方统计量为：" + str(chi2))
    Chi2.append(chi2)

# 作图
plt.rcParams["font.sans-serif"] = 'SimHei'
plt.rcParams['axes.unicode_minus'] = False
plt.hist(x=Chi2, bins=30, density=True, color="steelblue", edgecolor="black")
x = np.linspace(0, 30, 100000)
y = stats.chi2.pdf(x, df=9)
plt.plot(x, y, c='red')
plt.xlabel("卡方统计量")
plt.ylabel("概率密度(频率)")
plt.title("卡方分布 PDF")
plt.show()
