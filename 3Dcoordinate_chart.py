import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib

matplotlib.rcParams['font.family'] = 'Times New Roman'
matplotlib.rcParams['axes.unicode_minus'] = False  # 确保负号显示正确

file_paths = [
    r'C:\Users\XavierWynne\Desktop\wavelets2017\MR&SDR1007_wf3.xlsx',
    r'C:\Users\XavierWynne\Desktop\wavelets2017\MR&SDR1007_wf4.xlsx',
    r'C:\Users\XavierWynne\Desktop\wavelets2017\MR&SDR1007_wf5.xlsx',
    r'C:\Users\XavierWynne\Desktop\wavelets2017\MR&SDR1007_wf6.xlsx'
]


scales = [3, 4, 5, 6]
wavelengths = np.arange(350, 2501)

# 创建一个新的图和3D子图
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(111, projection='3d')

X, Y, Z = [], [], []

# 加载每个尺度的数据
for scale, file_path in zip(scales, file_paths):
    data = pd.read_excel(file_path, usecols=['Wavelength', 'MR'])
    X.extend(data['Wavelength'].values)
    Y.extend([scale] * len(data['Wavelength']))  # 为每个波长重复尺度值
    Z.extend(data['MR'].values)

sc = ax.scatter(X, Y, Z, c=Z, cmap='viridis', depthshade=True, label=scales)

# 设置坐标轴范围和标签
ax.set_xlabel('Wavelength (nm)', fontsize=22,labelpad=15)
ax.set_ylabel('Scale', fontsize=22,labelpad=10) #labelpad坐标轴间隔
ax.set_zlabel('Wavelet Power', fontsize=22,labelpad=10)

# # 设置坐标轴刻度
# ax.set_xlim(350, 2500)
# ax.set_ylim(min(scales), max(scales))
# # 根据小波功率数据调整Z轴范围
# ax.set_zlim(min(Z), max(Z))
# 设置坐标轴刻度
ax.set_xlim(350, 2500)  # 将X轴的起点设置为400
ax.set_ylim(min(scales), max(scales))  # Y轴范围已经是3到6，无需更改
ax.set_xticks([500, 750, 1000, 1250,1500,1750,2000,2250,2500])
ax.set_yticks([3, 4, 5, 6])  # 只显示Y轴的刻度3, 4, 5, 6
ax.set_zticks([-1, -0.5, 0, 0.5, 1])  # 设置Y轴的刻度，可以根据需要调整
# 根据小波功率数据调整Z轴范围
#ax.set_zlim(min(Z), max(Z))
ax.set_zlim(-1.5, 1.5)
# 设置刻度粗细和坐标轴的粗细
ax.tick_params(axis='both', which='major', width=2,labelsize=14)
ax.spines['bottom'].set_linewidth(2)
ax.spines['left'].set_linewidth(2)
ax.spines['right'].set_linewidth(2)
ax.spines['top'].set_linewidth(2)

# # 添加颜色条
# cbar = fig.colorbar(sc, ax=ax, pad=0.1)
# cbar.set_label('Wavelet Power')

# # 添加图例
# ax.legend()

# 显示图形
plt.show()