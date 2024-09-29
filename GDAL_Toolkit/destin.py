import os
import math
import scipy.signal as signal
import os
os.environ['PROJ_LIB'] = r'C:\Users\A\anaconda3\envs\Geo\Library\share\proj'
os.environ['GDAL_DATA'] = r'C:\Users\A\anaconda3\envs\Geo\Library\share'
from osgeo import osr, ogr
from osgeo import gdal
import numpy as np
import cv2

class destin:
    def __init__(self, img_path, outImgPath1, intKernelSize1=None, intKernelSize2=None, blocksize=None, closing = None):
        self.img_path = img_path
        self.closing = closing
        self.outImgPath1 = outImgPath1
        self.intKernelSize1 = intKernelSize1
        self.intKernelSize2 = intKernelSize2
        self.blocksize = blocksize
        self.read_file()
        self.creat_tempdir()

    def creat_tempdir(self):
        self.temp_name = os.path.join(self.outImgPath1, 'temp')
        if not os.path.exists(self.temp_name):
            os.makedirs(self.temp_name)

    def read_file(self):
        dataset2 = gdal.Open(self.img_path)
        x_size2 = dataset2.RasterXSize
        y_size2 = dataset2.RasterYSize
        self.data = dataset2.ReadAsArray(0, 0, x_size2, y_size2)


    def get_graphic(self, data, raster_band: int) -> np.ndarray:
        return data[raster_band - 1, :, :]



    def DESTIN_edge_strength(self, data):
        x_size = data.shape[0]
        y_size = data.shape[1]
        data = data.astype("float64")
        strength = np.zeros((x_size, y_size))
        factor = math.sqrt(2) / 2

        strength[0, 0] = abs(data[0, 0] - data[1, 0]) + abs(factor * (data[0, 0] - data[1, 1])) + abs(data[0, 0] - data[0, 1])
        strength[0, y_size - 1] = abs(data[0, y_size - 1] - data[1, y_size - 1]) + abs(factor * (data[0, y_size - 1] - data[1, y_size - 2]))   + abs(data[0, y_size - 1] - data[0, y_size - 2])
        strength[x_size - 1, 0] = abs(data[x_size - 1, 0] - data[x_size - 2, 0]) + abs(factor * (data[x_size - 1, 0] - data[x_size - 2, 1]))   + abs(data[x_size - 1, 0] - data[x_size - 1, 1])
        strength[x_size - 1, y_size - 1] = abs(data[x_size - 1, y_size - 1] - data[x_size - 2, y_size - 1])   + abs(
            factor * (data[x_size - 1, y_size - 1] - data[x_size - 2, y_size - 2])) \
                                           + abs(data[x_size - 1, y_size - 1] - data[x_size - 1, y_size - 2])
        for j in range(1, y_size - 1):
            strength[0, j] = abs(data[0, j] - data[0, j - 1])  + abs(factor * (data[0, j] - data[1, j - 1])) + abs(data[0, j] - data[1, j])  + abs(factor * (data[0, j] - data[1, j + 1]))  + abs(data[0, j] - data[0, j + 1])
            strength[x_size - 1, j] = abs(data[x_size - 1, j] - data[x_size - 1, j - 1])  + abs(factor * (data[x_size - 1, j] - data[x_size - 2, j - 1])) + abs(data[x_size - 1, j] - data[x_size - 2, j])  + abs(factor * (data[x_size - 1, j] - data[x_size - 2, j + 1]))   + abs(data[x_size - 1, j] - data[x_size - 1, j + 1])
        for i in range(1, x_size - 1):
            strength[i, 0] = abs(data[i, 0] - data[i, 1])   + abs(factor * (data[i, 0] - data[i - 1, 1])) + abs(data[i, 0] - data[i - 1, 1])  + abs(factor * (data[i, 0] - data[i + 1, 1]))  + abs(data[i, 0] - data[i + 1, 0])
            strength[i, y_size - 1] = abs(data[i, y_size - 1] - data[i, y_size - 2])  + abs(factor * (data[i, y_size - 1] - data[i - 1, y_size - 2]))  + abs(data[i, y_size - 1] - data[i - 1, y_size - 1]) + abs(factor * (data[i, y_size - 1] - data[i + 1, y_size - 2]))  + abs(data[i, y_size - 1] - data[i + 1, y_size - 1])
        for i in range(1, x_size - 1):
            for j in range(1, y_size - 1):
                # print("123")
                strength[i, j] = abs(data[i, j] - data[i + 1, j])  + abs(data[i, j] - data[i - 1, j]) + abs(data[i, j] - data[i, j + 1]) + abs(data[i, j] - data[i, j - 1]) + abs(factor * (data[i, j] - data[i + 1, j + 1])) + abs(factor * (data[i, j] - data[i - 1, j + 1]))  + abs(factor * (data[i, j] - data[i + 1, j - 1])) + abs(factor * (data[i, j] - data[i - 1, j - 1]))
        return strength

    def readTiff(self, fileName):
        self.dataset = gdal.Open(fileName)
        if self.dataset == None:
            print(fileName + "File cannot open!!")
            return
        self.cols = self.dataset.RasterXSize
        self.rows = self.dataset.RasterYSize
        return self.dataset, self.cols, self.rows

    def writeTiff(self, im_data2, im_width, im_height, im_bands, im_geotrans, im_proj, path):
        if 'int8' in im_data2.dtype.name:
            datatype = gdal.GDT_Byte
        elif 'int16' in im_data2.dtype.name:
            datatype = gdal.GDT_UInt16
        else:
            datatype = gdal.GDT_Float64
        if len(im_data2.shape) == 3:
            im_bands, im_height, im_width = im_data2.shape
        elif len(im_data2.shape) == 2:
            im_data2 = np.array([im_data2])
        else:
            im_bands, (im_height, im_width) = 1, im_data2.shape
        driver = gdal.GetDriverByName("GTiff")
        dataset = driver.Create(path, im_width, im_height, im_bands, datatype)
        if (dataset != None):
            dataset.SetGeoTransform(im_geotrans)
            dataset.SetProjection(im_proj)
        for i in range(im_bands):
            dataset.GetRasterBand(i + 1).WriteArray(im_data2[i])
        del dataset

    def tifInput2jpgShow(self, tifPath):
        parts1 = tifPath.split('/')[:-1]
        pic_space = '/'.join(parts1)+'/'
        parts2 = tifPath.split('/')[-1]
        new_jpg_name = parts2.split('.')[0] + '.png'
        new_jpg_path = new_jpg_name
        picNameWithoutSuffix = parts2.split('.')[0]
        return new_jpg_path, pic_space, picNameWithoutSuffix

    def showSavedJpg(self, img, tifPath):
        jpgPath, _, _ = self.tifInput2jpgShow(tifPath)
        cv2.imwrite(jpgPath, img * 255)
        return jpgPath

    def showSavedJpg2(self, img,tifPath):
        jpgPath = self.tif2png(tifPath)
        cv2.imwrite(jpgPath,img * 255)
        return jpgPath

    def tif2png(self, tigPath):
        from pathlib import Path
        pngpath = Path(tigPath).with_name("step2_png")
        return str(pngpath)

    def fastFilter(self, ct_data, size1):
        ct_data = signal.medfilt2d(ct_data, size1)
        return ct_data

    def genEdgeIntensityMap(self):
        """
        第1步 生成边缘强度图
        """
        # im_data = self.read_file(path1)

        #### 单波段输入####
        NDVI = self.data
        #### 单波段输入####


        # bandVal1, bandVal2 = 1, 2
        # BN = self.get_graphic(self.data, bandVal1)
        # BR = self.get_graphic(self.data, bandVal2)
        # NDVI = (BN-BR)/((BN+BR)+0.000001)
        # NDVI[NDVI > 5] = 0
        # NDVI[NDVI < -5] = 0
        current_data = signal.medfilt2d(NDVI, self.intKernelSize1)
        current_data = self.DESTIN_edge_strength(current_data)
        current_data = signal.medfilt2d(current_data, self.intKernelSize2)
        # current_data[current_data > 65535] = 0
        # Flaten_data = current_data.reshape((1, current_data.shape[0]*current_data.shape[1]))
        # min_index = int(0.10*Flaten_data.shape[1])
        # max_index = int(0.90*Flaten_data.shape[1])
        # Flaten_data = np.sort(Flaten_data)
        # min_data = Flaten_data[:, min_index]
        # max_data = Flaten_data[:, max_index]
        # new_current_data = (current_data-min_data)/(max_data-min_data+0.000001)
        self.ds, self.c, self.r = self.readTiff(self.img_path)
        self.basename1 = os.path.join(self.temp_name, 'step1.tif')
        self.writeTiff(current_data, self.c, self.r, 1, self.ds.GetGeoTransform(), self.ds.GetProjection(), self.basename1)
        self.jpg_Path = self.showSavedJpg(current_data, self.basename1)

    def raster2poly(self, raster, outshp):
        inraster = gdal.Open(raster)  # 读取路径中的栅格数据
        inband = inraster.GetRasterBand(1)  # 这个波段就是最后想要转为矢量的波段，如果是单波段数据的话那就都是1
        prj = osr.SpatialReference()
        prj.ImportFromWkt(inraster.GetProjection())  # 读取栅格数据的投影信息，用来为后面生成的矢量做准备

        drv = ogr.GetDriverByName("ESRI Shapefile")
        if os.path.exists(outshp):  # 若文件已经存在，则删除它继续重新做一遍
            drv.DeleteDataSource(outshp)
        Polygon = drv.CreateDataSource(outshp)  # 创建一个目标文件
        Poly_layer = Polygon.CreateLayer(raster[:-4], srs=prj, geom_type=ogr.wkbMultiPolygon)  # 对shp文件创建一个图层，定义为多个面类
        newField = ogr.FieldDefn('value', ogr.OFTReal)  # 给目标shp文件添加一个字段，用来存储原始栅格的pixel value
        Poly_layer.CreateField(newField)

        gdal.FPolygonize(inband, None, Poly_layer, 0)  # 核心函数，执行的就是栅格转矢量操作
        Polygon.SyncToDisk()
        Polygon = None


    def runStep2(self):
        """
        第2步 自适应阈值分割
        """
        img = cv2.imread(self.jpg_Path)  #   这个png是不是第一步生成的？  然后看一下栅格转Shp里面显示的有个库错了，gdal
        gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
        clahe = cv2.createCLAHE(clipLimit=0.01, tileGridSize=(100, 100))
        dst = clahe.apply(gray)
        self.thresh = cv2.adaptiveThreshold(dst, 255, cv2.ADAPTIVE_THRESH_MEAN_C, cv2.THRESH_BINARY_INV, self.blocksize, 0)
        self.basename2 = os.path.join(self.outImgPath1, 'step2-thresh'+str(self.blocksize)+'.tif')
        self.writeTiff(self.thresh, self.c, self.r, 1, self.ds.GetGeoTransform(), self.ds.GetProjection(), self.basename2)    # 输出路径

    def Closing(self):
        """
        第3步 闭运算
        """
        kernel = np.ones((5, 5), np.uint8)
        closing = cv2.morphologyEx(self.thresh, cv2.MORPH_CLOSE, kernel)
        # self.basename2 = os.path.join(outImgPath1, 'step2-thresh' + str(self.blocksize) + 'closing.tif')
        self.writeTiff(closing, self.c, self.r, 1, self.ds.GetGeoTransform(), self.ds.GetProjection(),
                       self.basename2)


    def destin_run(self):
        self.genEdgeIntensityMap()
        self.runStep2()
        out_shp_file = os.path.join(self.temp_name, 'shp_step1.shp')
        if self.closing == True:
            self.Closing()
            self.raster2poly(self.basename2, out_shp_file)
        else:
            self.raster2poly(self.basename2, out_shp_file)
        return self.temp_name


if __name__ == "__main__":
    path1 = r"D:\SR\Planet\XH\TS\Results\1005\G.tif"# 输入路径
    outImgPath1 = r"D:\SR\Planet\XH\TS\Results\1005\des"    # 输出路径
    obj_destin = destin(path1, outImgPath1, 5, 5, 101, closing=False)
    current_data = obj_destin.destin_run()