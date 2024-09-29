#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
from osgeo import gdal, ogr
from Read_write_data import readTiff,  writeTiff, get_graphic
import geopandas as gpd
import numpy as np
import rasterio as rio
import rasterio.mask


class add_mean_value_field:
    def __init__(self, shp_file, input_raster_file1, input_raster_file2, outpath, NDVI_gap=None):
        """
        :param shp_file: Destin分割边界转矢量
        :param input_raster_file1: NDVI值较大的图像文件
        :param input_raster_file2: NDVI值较小的图像文件
        :param outpath: 输出路径
        """
        self.shp_file = shp_file
        self.raster_file1 = input_raster_file1
        self.raster_file2 = input_raster_file2
        self.outpath = outpath
        self.NDVI_gap = NDVI_gap
        self.mean_name1 = 'NDVI1'
        self.mean_name2 = 'NDVI2'

    def creat_tempdir(self):
        self.temp_name = os.path.join(self.outpath, 'temp1')
        if not os.path.exists(self.temp_name):
            os.makedirs(self.temp_name)

    def basic_data_name(self):
        self.NDVI1_name = os.path.join(self.temp_name, 'NDVI1.tif')
        self.NDVI2_name = os.path.join(self.temp_name, 'NDVI2.tif')
        self.out_shp_file1 = os.path.join(self.outpath, 'shp1.shp')
        self.out_shp_file2 = os.path.join(self.temp_name, 'shp2.shp')

    def Cal_NDVI(self, raster_file):
        """
        计算遥感图像的NDVI
        :param raster_file: 输入的4波段图像
        :return
        """
        dataset, cols, rows = readTiff(raster_file)
        data = dataset.ReadAsArray(0, 0, cols, rows)
        self.geos, self.pros = dataset.GetGeoTransform(), dataset.GetProjection()
        bandVal1, bandVal2 = 4, 3
        BN = get_graphic(data, bandVal1)
        BR = get_graphic(data, bandVal2)
        NDVI = (BN - BR) / (BN + BR + 0.00000000001)
        return NDVI

    def save_NDVI(self):
        self.NDVI1 = self.Cal_NDVI(self.raster_file1)
        self.NDVI2 = self.Cal_NDVI(self.raster_file2)
        self.NDVI_12 = self.NDVI1-self.NDVI2
        self.NDVI_12[self.NDVI_12 < 0.1] = 0.1
        # writeTiff(self.NDVI1, self.NDVI1.shape[1], self.NDVI1.shape[0], 1, self.geos, self.pros, self.NDVI1_name)
        # writeTiff(self.NDVI2, self.NDVI2.shape[1], self.NDVI2.shape[0], 1, self.geos, self.pros, self.NDVI2_name)
        writeTiff(self.NDVI_12, self.NDVI_12.shape[1], self.NDVI_12.shape[0], 1, self.geos, self.pros, self.NDVI1_name)


#计算栅格均值
    def add_mean_value_field(self, raster_file, shp_file, out_shp_file):
        """
        计算矢量图斑范围内的栅格均值并增加字段
        :param raster_file: 待计算均值的栅格底图
        :param shp_file: 与栅格数据叠加的矢量掩膜数据
        :param out_shp_file: 增加均值字段后输出的矢量数据，可与shp_file相同，即重写shp_file文件
        :return
        """

        shp_data = gpd.GeoDataFrame.from_file(shp_file)
        raster_data = rio.open(raster_file)
        profile = raster_data.profile
        # 保存投影信息一致
        shp_data = shp_data.to_crs(raster_data.crs)

        out_shp_data = shp_data.copy()
        mean_value_field = []
        for i in range(0, len(shp_data)):
            # 获取矢量数据的features
            geo = shp_data.geometry[i]
            feature = [geo.__geo_interface__]
            # 通过feature裁剪栅格影像
            out_image, out_transform = rasterio.mask.mask(raster_data, feature, all_touched=True, crop=True,
                                                     nodata=raster_data.nodata)
            # 获取影像Value值，并转化为list
            out_list = out_image.data.tolist()
            # 除去list中的Nodata值
            out_list = out_list[0]
            out_data = []
            for k in range(len(out_list)):
                for j in range(len(out_list[k])):
                    if out_list[k][j] > 0:
                        out_data.append(out_list[k][j])
            # 求数据的平均数
            if len(out_data):
                mean_data = np.mean(out_data)
            else:
                mean_data = None
            mean_value_field.append(mean_data)

        # 增加属性字段，并将GeodataFrame导出为shp文件
        out_shp_data.insert(out_shp_data.shape[1], 'Gap_NDVI', mean_value_field)
        out_shp_data.to_file(out_shp_file)

    def object_class1(self, inshp):
        gdal.AllRegister()
        # 解决中文路径乱码问题
        gdal.SetConfigOption("GDAL_FILENAME_IS_UTF8", "NO")
        driver = ogr.GetDriverByName('ESRI Shapefile')
        pFeatureDataset = driver.Open(inshp, 1)  # 英文路径
        pFeaturelayer = pFeatureDataset.GetLayer(0)

        # 按条件查询空间要素，本例查询字段名为Value，字段值为0的所有要素。
        strValue = 0
        strFilter = "Value = '" + str(strValue) + "'"
        pFeaturelayer.SetAttributeFilter(strFilter)
        #
        # # 删除第二部查询到的矢量要素，注意，此时获取到的Feature皆为选择的Feature.
        pFeatureDef = pFeaturelayer.GetLayerDefn()
        pLayerName = pFeaturelayer.GetName()
        pFieldName = "Value"
        pFieldIndex = pFeatureDef.GetFieldIndex(pFieldName)
        for pFeature in pFeaturelayer:
            pFeatureFID = pFeature.GetFID()
            pFeaturelayer.DeleteFeature(int(pFeatureFID))
        strSQL = "REPACK " + str(pFeaturelayer.GetName())
        pFeatureDataset.ExecuteSQL(strSQL, None, "")
        pFeatureLayer = None
        pFeatureDataset = None

    def object_class2(self, inshp, NDVI_gap):
        gdal.AllRegister()
        # 解决中文路径乱码问题
        gdal.SetConfigOption("GDAL_FILENAME_IS_UTF8", "NO")
        driver = ogr.GetDriverByName('ESRI Shapefile')
        pFeatureDataset = driver.Open(inshp, 1)  # 英文路径
        pFeaturelayer = pFeatureDataset.GetLayer(0)

        # 按条件查询空间要素，本例查询字段名为Value，字段值为0的所有要素。
        strValue = NDVI_gap
        strFilter = "Gap_NDVI <= '" + str(strValue) + "'"
        pFeaturelayer.SetAttributeFilter(strFilter)
        #
        # # 删除第二部查询到的矢量要素，注意，此时获取到的Feature皆为选择的Feature.
        pFeatureDef = pFeaturelayer.GetLayerDefn()
        pLayerName = pFeaturelayer.GetName()
        pFieldName = "Gap_NDVI"
        pFieldIndex = pFeatureDef.GetFieldIndex(pFieldName)
        for pFeature in pFeaturelayer:
            pFeatureFID = pFeature.GetFID()
            pFeaturelayer.DeleteFeature(int(pFeatureFID))
        strSQL = "REPACK " + str(pFeaturelayer.GetName())
        pFeatureDataset.ExecuteSQL(strSQL, None, "")
        pFeatureLayer = None
        pFeatureDataset = None




    def run(self):
        self.creat_tempdir()
        self.basic_data_name()
        self.save_NDVI()
        self.add_mean_value_field(self.NDVI1_name, self.shp_file, self.out_shp_file1)
        self.object_class1(self.out_shp_file1)
        self.object_class2(self.out_shp_file1, self.NDVI_gap)
        # if os.path.exists(self.temp_name):
        #     shutil.rmtree(self.temp_name)
        # self.add_mean_value_field(self.NDVI2_name, self.out_shp_file1, self.mean_name2, self.out_shp_file2)




if __name__ == '__main__':
    shp_file = r"D:\SR\Planet\XH\TS\Results\temp\32\D__SR_Planet_XH_TS_Results_20230812_pre32.shp"
    raster_file1 = r"D:\SR\Planet\XH\TS\NDVI\20230812_023028_22_249c_3B_AnalyticMS_clip.tif"
    raster_file2 = r"D:\SR\Planet\XH\TS\NDVI\20230609_022556_44_2473_3B_AnalyticMS_clip.tif"
    outpath = r"C:\Users\A\Desktop\output\output"
    obj = add_mean_value_field(shp_file, raster_file1, raster_file2, outpath, NDVI_gap=0.3)
    obj.run()
