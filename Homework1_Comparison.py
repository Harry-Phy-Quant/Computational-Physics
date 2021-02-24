from Homework1_Cholesky import cholesky
from Homework1_GEM import GEM
import matplotlib.pyplot as plt
from math import *


x_data = []
y_data = []

for n in range(2, 14):
    x_data.append(n)
    difference_list = [GEM(n)[i] - cholesky(n)[i] for i in range(n)]
    difference = sqrt(sum([i ** 2 for i in difference_list]))
    y_data.append(difference)

#y_data = [log10(i) for i in y_data]  # 取对数坐标

plt.plot(x_data, y_data)
plt.show()
