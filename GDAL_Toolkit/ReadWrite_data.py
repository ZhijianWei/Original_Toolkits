from osgeo import gdal
import numpy as np



def readTiff(fileName):
    dataset = gdal.Open(fileName)
    if dataset == None:
        print(fileName + "File cannot open!!")
        return
    cols = dataset.RasterXSize
    rows = dataset.RasterYSize
    return dataset, cols, rows


def writeTiff(im_data2, im_width, im_height, im_bands, im_geotrans, im_proj, path):
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

def get_graphic(data, raster_band: int) -> np.ndarray:
    return data[raster_band - 1, :, :]