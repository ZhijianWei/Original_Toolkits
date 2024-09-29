import glob
import os
import cv2
from osgeo import gdal
import numpy as np
from PIL import Image

class Topng():
    def __init__(self, in_img_file):
        """
        :param in_img_file: 输入的4波段遥感图像
        @ created by Sun haoran
        """
        self.in_img_file = in_img_file

    def readTif(self):  # 读取红绿蓝（RGB）三波段数据
        bandsOrder = [1, 2, 3]
        dataset = gdal.Open(self.in_img_file, gdal.GA_ReadOnly)  # 返回一个gdal.Dataset类型的对象
        cols = dataset.RasterXSize  # tif图像的宽度
        rows = dataset.RasterYSize  # tif图像的高度
        data = np.empty([rows, cols, 3], dtype=float)  # 定义结果数组，将RGB三波段的矩阵存储
        for i in range(3):
            band = dataset.GetRasterBand(bandsOrder[i])  # 读取波段数值
            oneband_data = band.ReadAsArray()  # 读取波段数值读为numpy数组
            data[:, :, i] = oneband_data  # 将读取的结果存放在三维数组的一页三
        return data


    def tig_to_jpg(self, bandData, lower_percent=0.5, higher_percent=99.5):
        # banddata为读取的3、2、1波段数据
        band_Num = bandData.shape[2]           # 数组第三维度的大小，在这里是图像的通道数
        JPG_Array = np.zeros_like(bandData, dtype=np.uint8)
        for i in range(band_Num):
            minValue = 0
            maxValue = 255
            #获取数组RGB_Array某个百分比分位上的值
            low_value = np.percentile(bandData[:, :,i], lower_percent)
            high_value = np.percentile(bandData[:, :,i], higher_percent)
            temp_value =(bandData[:, :,i] - low_value) * maxValue / (high_value - low_value)
            temp_value[temp_value < minValue] = minValue
            temp_value[temp_value > maxValue] = maxValue
            JPG_Array[:, :, i] = temp_value
        outputImg = Image.fromarray(np.uint8(bandData))
        return JPG_Array


if __name__ == '__main__':
    in_img_file = r"C:\Users\A\Desktop\output\YH0323\HA_wv_sun.tif"
    out_img_file = r"C:\Users\A\Desktop\output\YH0323\HA_wv_SunHR.png"
    obj_png = Topng(in_img_file)
    data = obj_png.readTif()
    JPG_Array = obj_png.tig_to_jpg(data)
    cv2.imwrite(out_img_file, JPG_Array)




    # imgPath = r"D:\DeepLearning_RS_data\KDXF_VHR\ARS_VHR\im"
    # newPath = r"D:\DeepLearning_RS_data\KDXF_VHR\ARS_VHR\true_im"
    # file_name = glob.glob(os.path.join(imgPath, '*tif'))
    # for name in file_name:
    #     code_name = name.split('\\')[-1].split('.')[0]
    #     new_name = os.path.join(newPath, code_name+'.png')
    #     data = readTif(name)
    #     outputImg = tig_to_jpg(data, lower_percent=0.5, higher_percent=99.5)
    #     cv2.imwrite(new_name, outputImg)