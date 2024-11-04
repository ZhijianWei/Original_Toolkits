import numpy as np
import pandas as pd
from osgeo import gdal, ogr, osr

def getWorksheet(path,index):
    workbook =pd.read_excel(path)
    worksheet = workbook.iloc[:,index]
    return worksheet



def getShp(lon,lat,worksheet):
    nrows = worksheet.nrows
    ncols = worksheet.ncols
    gdal.SetConfigOption('GDAl_FILENAME_IS_UTF8','NO')
    gdal.SetConfigOption('SHAPE_ENCODING','GB2312')
    ogr.RegisterAll() #注册所有驱动

    #生成坐标系
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326) #坐标参考

    #生成图层
    driver=ogr.GetDriverByName("ESRI Shapefile")
    dataSource = driver.CreatedataSource('test.shp')
    layer= dataSource.CreateLayer('test',srs,ogr.wkbPoint) #wkbpoint为文件类型

    field_list=[]

    for i in range(0, ncols):
        print(worksheet.cell_value(0,i))
        field_name=worksheet.cell_value(0,i)
        # 生成layer中的feature
        field = ogr.FieldDefn('name', ogr.OFTString)  # 生成属性表的字段
        field.SetWidth(100)
        layer.CreateField(field)

    for i in range(1,nrows):
        feature = ogr.Feature(layer.GetLayerDefn())
        geom = ogr.Geometry(ogr.wkbPoint)
        geom.AddPoint(float(worksheet.cell_value(i,lon)), float(worksheet.cell_value(i,lat)))
        feature.SetGeometry(geom)
        for j in range(0,ncols):
            feature.SetField(field_list[j],worksheet.cell_value(i,j)) #给要素生成几何体
        layer.CreateFeature(feature)

worksheet=getWorksheet('ChinaCity.xlsx',0)
getShp(4,3,worksheet)




