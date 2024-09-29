from osgeo import gdal, ogr


def readShp(shp_name):
    gdal.SetConfigOption("GDAL_FILENAME_IS_UTF8", "YES")
    gdal.SetConfigOption("SHAPE_ENCODING", "UTF-8")
    ogr.RegisterAll()
    ds = ogr.Open(shp_name, 0)
    if ds == None:
        return "打开文件失败！"
    iLayerCount = ds.GetLayerCount()
    oLayer = ds.GetLayerByIndex(0)
    if oLayer == None:
        return "获取图层失败！"
    oLayer.ResetReading()
    num = oLayer.GetFeatureCount(0)
    result_list = []
    # 获取要素
    for i in range(0, num):
        ofeature = oLayer.GetFeature(i)
        id = ofeature.GetFieldAsString("id")
        name = ofeature.GetFieldAsString('name')
        geom = str(ofeature.GetGeometryRef())
        result_list.append([id, name, geom])
    return ds
