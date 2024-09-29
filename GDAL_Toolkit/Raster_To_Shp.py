#矢量转栅格shp
#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os

#import gdal

import geopandas as gpd
import numpy as np
import rasterio as rio
from osgeo import osr, ogr
from osgeo import gdal


def raster2poly(raster, outshp):
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

if __name__ == "__main__":
    raser_file = r"D:\SR\Planet\XH\TS\Results\0220\ImageXH_SR_pre.tif"
    shp_file = r"D:\SR\Planet\XH\TS\Results\0220\shp"
    raster2poly(raser_file, shp_file)