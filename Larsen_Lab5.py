import numpy as np
import rasterio
import numpy.ma as ma
import pandas as pd
import copy
from rasterstats import zonal_stats

path = '/Users/hannahlarsen/Desktop/UC Denver/FA 2022/Python and GIS/lab5_data/'
path_rasters = '/Users/hannahlarsen/Desktop/UC Denver/FA 2022/Python and GIS/lab5_data/L5_big_elk/'

def slopeAspect(dem, cs):
    """Calculates slope and aspect using the 3rd-order finite difference method
    Parameters
    ----------
    dem : numpy array
        A numpy array of a DEM
    cs : float
        The cell size of the original DEM
    Returns
    -------
    numpy arrays
        Slope and Aspect arrays
    """
    from math import pi
    from scipy import ndimage
    kernel = np.array([[-1, 0, 1], [-2, 0, 2], [-1, 0, 1]])
    dzdx = ndimage.convolve(dem, kernel, mode='mirror') / (8 * cs)
    dzdy = ndimage.convolve(dem, kernel.T, mode='mirror') / (8 * cs)
    slp = np.arctan((dzdx ** 2 + dzdy ** 2) ** 0.5) * 180 / pi
    ang = np.arctan2(-dzdy, dzdx) * 180 / pi
    aspect = np.where(ang > 90, 450 - ang, 90 - ang)
    return slp, aspect

def reclassAspect(npArray):
    """Reclassify aspect array to 8 cardinal directions (N,NE,E,SE,S,SW,W,NW),
    encoded 1 to 8, respectively (same as ArcGIS aspect classes).
    Parameters
    ----------
    npArray : numpy array
        numpy array with aspect values 0 to 360
    Returns
    -------
    numpy array
        numpy array with cardinal directions
    """
    return np.where((npArray > 22.5) & (npArray <= 67.5), 2,
    np.where((npArray > 67.5) & (npArray <= 112.5), 3,
    np.where((npArray > 112.5) & (npArray <= 157.5), 4,
    np.where((npArray > 157.5) & (npArray <= 202.5), 5,
    np.where((npArray > 202.5) & (npArray <= 247.5), 6,
    np.where((npArray > 247.5) & (npArray <= 292.5), 7,
    np.where((npArray > 292.5) & (npArray <= 337.5), 8, 1)))))))

def reclassByHisto(npArray, bins):
    """Reclassify np array based on a histogram approach using a specified
    number of bins. Returns the reclassified numpy array and the classes from
    the histogram.
    Parameters
    ----------
    npArray : numpy array
        Array to be reclassified
    bins : int
        Number of bins
    Returns
    -------
    numpy array
        numpy array with reclassified values
    """
    # array = np.where(np.isnan(npArray), 0, npArray)
    histo = np.histogram(~np.isnan(npArray), bins)[1]
    rClss = np.zeros_like(npArray)
    for i in range(bins):
        print(i + 1, histo[i], histo[i + 1])
        print(np.where((npArray > histo[i]) & (npArray <= histo[i + 1])))
        rClss = np.where((npArray >= histo[i]) & (npArray <= histo[i + 1]),
                         i + 1, rClss)
    return rClss


DEM = rasterio.open(path + 'bigElk_dem.tif')
DEM_read = DEM.read(1)

DEM_arr = np.array(DEM_read)

affine = DEM.transform
cell_size = affine[0]
bins = 10

slope_aspect = slopeAspect(DEM_arr, cell_size)

slope = slope_aspect[0]
aspect = slope_aspect[1]

reclassify_aspect = reclassAspect(aspect)

slope_reclassify = reclassByHisto(slope, bins)


fire_perim_read = rasterio.open(path + 'fire_perimeter.tif')
fire_perim = fire_perim_read.read(1)
    
b3_2002_raster = rasterio.open(path_rasters + 'L5034032_2002_B3.tif')
b3_2002 = b3_2002_raster.read(1)
mean3_2002 = np.mean(b3_2002[b3_2002 != 255])
b4_2002 = rasterio.open(path_rasters + 'L5034032_2002_B4.tif')
b4_2002 = b4_2002.read(1)
mean4_2002 = np.mean(b4_2002[b4_2002 != 255])

b3_2003 = rasterio.open(path_rasters + 'L5034032_2003_B3.tif')
b3_2003 = b3_2003.read(1)
mean3_2003 = np.mean(b3_2003[b3_2003 != 255])
b4_2003 = rasterio.open(path_rasters + 'L5034032_2003_B4.tif')
b4_2003 = b4_2003.read(1)
mean4_2003 = np.mean(b4_2003[b4_2003 != 255])

b3_2004 = rasterio.open(path_rasters + 'L5034032_2004_B3.tif')
b3_2004 = b3_2004.read(1)
mean3_2004 = np.mean(b3_2002[b3_2002 != 255])
b4_2004 = rasterio.open(path_rasters + 'L5034032_2004_B4.tif')
b4_2004 = b4_2004.read(1)
mean4_2004 = np.mean(b4_2004[b4_2004 != 255])

b3_2005 = rasterio.open(path_rasters + 'L5034032_2005_B3.tif')
b3_2005 = b3_2005.read(1)
mean3_2005 = np.mean(b3_2005[b3_2005 != 255])
b4_2005 = rasterio.open(path_rasters + 'L5034032_2005_B4.tif')
b4_2005 = b4_2005.read(1)
mean4_2005 = np.mean(b4_2005[b4_2005 != 255])

b3_2006 = rasterio.open(path_rasters + 'L5034032_2006_B3.tif')
b3_2006 = b3_2006.read(1)
mean3_2006 = np.mean(b3_2006[b3_2006 != 255])
b4_2006 = rasterio.open(path_rasters + 'L5034032_2006_B4.tif')
b4_2006 = b4_2006.read(1)
mean4_2006 = np.mean(b4_2006[b4_2006 != 255])

b3_2007 = rasterio.open(path_rasters + 'L5034032_2007_B3.tif')
b3_2007 = b3_2007.read(1)
mean3_2007 = np.mean(b3_2007[b3_2007 != 255])
b4_2007 = rasterio.open(path_rasters + 'L5034032_2007_B4.tif')
b4_2007 = b4_2007.read(1)
mean4_2007 = np.mean(b4_2007[b4_2007 != 255])

b3_2008 = rasterio.open(path_rasters + 'L5034032_2008_B3.tif')
b3_2008 = b3_2008.read(1)
mean3_2008 = np.mean(b3_2008[b3_2008 != 255])
b4_2008 = rasterio.open(path_rasters + 'L5034032_2008_B4.tif')
b4_2008 = b4_2008.read(1)
mean4_2008 = np.mean(b4_2008[b4_2008 != 255])

b3_2009 = rasterio.open(path_rasters + 'L5034032_2009_B3.tif')
b3_2009 = b3_2009.read(1)
mean3_2009 = np.mean(b3_2009[b3_2009 != 255])
b4_2009 = rasterio.open(path_rasters + 'L5034032_2009_B4.tif')
b4_2009 = b4_2009.read(1)
mean4_2009 = np.mean(b4_2009[b4_2009 != 255])

b3_2010 = rasterio.open(path_rasters + 'L5034032_2010_B3.tif')
b3_2010 = b3_2010.read(1)
mean3_2010 = np.mean(b3_2010[b3_2010 != 255])
b4_2010 = rasterio.open(path_rasters + 'L5034032_2010_B4.tif')
b4_2010 = b4_2010.read(1)
mean4_2010 = np.mean(b4_2010[b4_2010 != 255])

b3_2011 = rasterio.open(path_rasters + 'L5034032_2011_B3.tif')
b3_2011 = b3_2011.read(1)
mean3_2011 = np.mean(b3_2011[b3_2011 != 255])
b4_2011 = rasterio.open(path_rasters + 'L5034032_2011_B4.tif')
b4_2011 = b4_2011.read(1)
mean4_2011 = np.mean(b4_2011[b4_2011 != 255])

band3_files = [mean3_2002, mean3_2003, mean3_2004, mean3_2005, mean3_2006, mean3_2007, mean3_2008, mean3_2009, mean3_2010, mean3_2011]
band4_files = [mean4_2002, mean4_2003, mean4_2004, mean4_2005, mean4_2006, mean4_2007, mean4_2008, mean4_2009, mean4_2010, mean4_2011]

list_NDVI = []
for band3, band4 in zip(band3_files, band4_files):
    NDVI = (band4 - band3) / (band4 + band3)
    list_NDVI.append(NDVI)


b3_2002_m = copy.copy(b3_2002)
b4_2002_m = copy.copy(b4_2002)

b3_2003_m = copy.copy(b3_2003)
b4_2003_m = copy.copy(b4_2003)

b3_2004_m = copy.copy(b3_2004)
b4_2004_m = copy.copy(b4_2004)

b3_2005_m = copy.copy(b3_2005)
b4_2005_m = copy.copy(b4_2005)

b3_2006_m = copy.copy(b3_2006)
b4_2006_m = copy.copy(b4_2006)

b3_2007_m = copy.copy(b3_2007)
b4_2007_m = copy.copy(b4_2007)

b3_2008_m = copy.copy(b3_2008)
b4_2008_m = copy.copy(b4_2008)

b3_2009_m = copy.copy(b3_2009)
b4_2009_m = copy.copy(b4_2009)

b3_2010_m = copy.copy(b3_2010)
b4_2010_m = copy.copy(b4_2010)

b3_2011_m = copy.copy(b3_2011)
b4_2011_m = copy.copy(b4_2011)

    
mask = ma.masked_where(fire_perim != 2, fire_perim)


b3_2002_masked = np.ma.masked_where(ma.getmask(mask), b3_2002_m)
b4_2002_masked = np.ma.masked_where(ma.getmask(mask), b4_2002_m)

b3_2003_masked = np.ma.masked_where(ma.getmask(mask), b3_2003_m)
b4_2003_masked = np.ma.masked_where(ma.getmask(mask), b4_2003_m)

b3_2004_masked = np.ma.masked_where(ma.getmask(mask), b3_2004_m)
b4_2004_masked = np.ma.masked_where(ma.getmask(mask), b4_2004_m)

b3_2005_masked = np.ma.masked_where(ma.getmask(mask), b3_2005_m)
b4_2005_masked = np.ma.masked_where(ma.getmask(mask), b4_2005_m)

b3_2006_masked = np.ma.masked_where(ma.getmask(mask), b3_2006_m)
b4_2006_masked = np.ma.masked_where(ma.getmask(mask), b4_2006_m)

b3_2007_masked = np.ma.masked_where(ma.getmask(mask), b3_2007_m)
b4_2007_masked = np.ma.masked_where(ma.getmask(mask), b4_2007_m)

b3_2008_masked = np.ma.masked_where(ma.getmask(mask), b3_2008_m)
b4_2008_masked = np.ma.masked_where(ma.getmask(mask), b4_2008_m)

b3_2009_masked = np.ma.masked_where(ma.getmask(mask), b3_2009_m)
b4_2009_masked = np.ma.masked_where(ma.getmask(mask), b4_2009_m)

b3_2010_masked = np.ma.masked_where(ma.getmask(mask), b3_2010_m)
b4_2010_masked = np.ma.masked_where(ma.getmask(mask), b4_2010_m)

b3_2011_masked = np.ma.masked_where(ma.getmask(mask), b3_2011_m)
b4_2011_masked = np.ma.masked_where(ma.getmask(mask), b4_2011_m)
    

b3_2002_data_final = b3_2002_masked.data
b4_2002_data_final = b4_2002_masked.data

b3_2003_data_final = b3_2003_masked.data
b4_2003_data_final = b4_2003_masked.data

b3_2004_data_final = b3_2004_masked.data
b4_2004_data_final = b4_2004_masked.data

b3_2005_data_final = b3_2005_masked.data
b4_2005_data_final = b4_2005_masked.data

b3_2006_data_final = b3_2006_masked.data
b4_2006_data_final = b4_2006_masked.data

b3_2007_data_final = b3_2007_masked.data
b4_2007_data_final = b4_2007_masked.data

b3_2008_data_final = b3_2008_masked.data
b4_2008_data_final = b4_2008_masked.data

b3_2009_data_final = b3_2009_masked.data
b4_2009_data_final = b4_2009_masked.data

b3_2010_data_final = b3_2010_masked.data
b4_2010_data_final = b4_2010_masked.data

b3_2011_data_final = b3_2011_masked.data
b4_2011_data_final = b4_2011_masked.data



pixels_2002 = []
for band3, band4 in np.nditer((b3_2002_data_final, b4_2002_data_final)):
    NDVI_pixel = (band4 - band3) / (band4 + band3)
    NDVI_final = (NDVI_pixel) / (list_NDVI[0])
    pixels_2002.append(NDVI_final)

pixels_2003 = []
for band3, band4 in np.nditer((b3_2003_data_final, b4_2003_data_final)):
    NDVI_pixel = (band4 - band3) / (band4 + band3)
    NDVI_final = (NDVI_pixel) / (list_NDVI[1])
    pixels_2003.append(NDVI_final)
    
pixels_2004 = []
for band3, band4 in np.nditer((b3_2004_data_final, b4_2004_data_final)):
    NDVI_pixel = (band4 - band3) / (band4 + band3)
    NDVI_final = (NDVI_pixel) / (list_NDVI[2])
    pixels_2004.append(NDVI_final)
    
pixels_2005 = []
for band3, band4 in np.nditer((b3_2005_data_final, b4_2005_data_final)):
    NDVI_pixel = (band4 - band3) / (band4 + band3)
    NDVI_final = (NDVI_pixel) / (list_NDVI[3])
    pixels_2005.append(NDVI_final)
    
pixels_2006 = []
for band3, band4 in np.nditer((b3_2006_data_final, b4_2006_data_final)):
    NDVI_pixel = (band4 - band3) / (band4 + band3)
    NDVI_final = (NDVI_pixel) / (list_NDVI[4])
    pixels_2006.append(NDVI_final)
    
pixels_2007 = []
for band3, band4 in np.nditer((b3_2007_data_final, b4_2007_data_final)):
    NDVI_pixel = (band4 - band3) / (band4 + band3)
    NDVI_final = (NDVI_pixel) / (list_NDVI[5])
    pixels_2007.append(NDVI_final)
    
pixels_2008 = []
for band3, band4 in np.nditer((b3_2008_data_final, b4_2008_data_final)):
    NDVI_pixel = (band4 - band3) / (band4 + band3)
    NDVI_final = (NDVI_pixel) / (list_NDVI[6])
    pixels_2008.append(NDVI_final)
    
pixels_2009 = []
for band3, band4 in np.nditer((b3_2009_data_final, b4_2009_data_final)):
    NDVI_pixel = (band4 - band3) / (band4 + band3)
    NDVI_final = (NDVI_pixel) / (list_NDVI[7])
    pixels_2009.append(NDVI_final)

pixels_2010 = []
for band3, band4 in np.nditer((b3_2010_data_final, b4_2010_data_final)):
    NDVI_pixel = (band4 - band3) / (band4 + band3)
    NDVI_final = (NDVI_pixel) / (list_NDVI[8])
    pixels_2010.append(NDVI_final)
    
pixels_2011 = []
for band3, band4 in np.nditer((b3_2011_data_final, b4_2011_data_final)):
    NDVI_pixel = (band4 - band3) / (band4 + band3)
    NDVI_final = (NDVI_pixel) / (list_NDVI[9])
    pixels_2011.append(NDVI_final)
    
pixels_2002 = np.array(pixels_2002)
pixels_2003 = np.array(pixels_2003)
pixels_2004 = np.array(pixels_2004)
pixels_2005 = np.array(pixels_2005)
pixels_2006 = np.array(pixels_2006)
pixels_2007 = np.array(pixels_2007)
pixels_2008 = np.array(pixels_2008)
pixels_2009 = np.array(pixels_2009)
pixels_2010 = np.array(pixels_2010)
pixels_2011 = np.array(pixels_2011)



list_2002 = np.full((280, 459), 2003)
list_2003 = np.full((280, 459), 2003)
list_2004 = np.full((280, 459), 2004)
list_2005 = np.full((280, 459), 2005)
list_2006 = np.full((280, 459), 2006)
list_2007 = np.full((280, 459), 2007)
list_2008 = np.full((280, 459), 2008)
list_2009 = np.full((280, 459), 2009)
list_2010 = np.full((280, 459), 2010)
list_2011 = np.full((280, 459), 2011)

x_2002 = list_2002.flatten()
x_2003 = list_2003.flatten()
x_2004 = list_2004.flatten()
x_2005 = list_2005.flatten()
x_2006 = list_2006.flatten()
x_2007 = list_2007.flatten()
x_2008 = list_2008.flatten()
x_2009 = list_2009.flatten()
x_2010 = list_2010.flatten()
x_2011 = list_2011.flatten()

polyfit = []
for a, b, c, d, e, f, g, h, i, j, q, r, s, t, u, v, w, x, y, z in zip(pixels_2002, 
        pixels_2003, pixels_2004, pixels_2005, pixels_2006, pixels_2007, pixels_2008, 
        pixels_2009, pixels_2010, pixels_2011, x_2002, x_2003, x_2004, x_2005, 
        x_2006, x_2007, x_2008, x_2009, x_2010, x_2011):
    temp_list = [a, b, c, d, e, f, g, h, i, j]
    x_list = [q, r, s, t, u, v, w, x, y, z]
    coeff = np.polyfit(x_list, temp_list, deg = 1)
    polyfit.append(coeff)


list_of_means = []
for values in [pixels_2002, pixels_2003, pixels_2004, pixels_2005, pixels_2006, 
              pixels_2007, pixels_2008, pixels_2009, pixels_2010, pixels_2011]:
    mean = sum(values) / len(values)
    list_of_means.append(mean)
    
print("The mean recovery ratio for the years 2002-2011 are: "+ str(list_of_means[0]) + 
      ", " + str(list_of_means[1]) + "," + str(list_of_means[2]) + ", " + 
      str(list_of_means[3]) + ", " + str(list_of_means[4]) + ", " + str(list_of_means[5])
      + ", " + str(list_of_means[6]) + ", " + str(list_of_means[7]) + ", " + 
      str(list_of_means[8]) + ", and " + str(list_of_means[9]) + ".")


polyfit_slope = [item[0] for item in polyfit]    
mean = sum(polyfit_slope) / len(polyfit_slope)

print("The mean coefficient of recovery from 2002-2011 is: "  + str(round((mean), 6)))

print("My conclusions regarding vegetation recovery is that shallower slopes" +
      " and aspects of south, southwest, and southeast have greater vegetation " + 
      "recovery after wildfires.")

polyfit_slope_array = np.array(polyfit_slope)

polyfit_slope_array.shape
polyfit_slope_2d = polyfit_slope_array.reshape(280, 459)

with rasterio.open(path + 'fire_perimeter.tif') as src:
    ras_meta = src.profile

with rasterio.open(path + 'coeff_recovery.tif', 'w', **ras_meta) as dst:
    dst.write(polyfit_slope_2d, indexes = 1)


affine_fire = fire_perim_read.transform


def zonal_stat(zone_raster, value_raster):
    """
    This function calculates the count, minimum, mean, max, and standard
    deviation for the value raster that is bounded by a specific zone from
    the zone raster
    """
    zonal_stats(zone_raster, value_raster, stats="count min mean max std", affine = affine_fire)



mask = ma.masked_where(fire_perim != 2, fire_perim)

reclassify_aspect_mask = np.ma.masked_where(ma.getmask(mask), reclassify_aspect)
reclassify_flatten = reclassify_aspect_mask.flatten()


slope_flatten = slope_reclassify.flatten() 
slope_int = list(map(int, slope_flatten))
slope_int_arr = np.array(slope_int)
slope_2d = slope_int_arr.reshape(280, 459)

slope_reclassify_mask = np.ma.masked_where(ma.getmask(mask), slope_2d)     

    
reclass_1 = zonal_stat((np.where(reclassify_flatten == 1)), polyfit_slope_2d)
reclass_2 = zonal_stat((np.where(reclassify_flatten == 2)), polyfit_slope_2d)
reclass_3 = zonal_stat((np.where(reclassify_flatten == 3)), polyfit_slope_2d)
reclass_4 = zonal_stat((np.where(reclassify_flatten == 4)), polyfit_slope_2d)
reclass_5 = zonal_stat((np.where(reclassify_flatten == 5)), polyfit_slope_2d)
reclass_6 = zonal_stat((np.where(reclassify_flatten == 6)), polyfit_slope_2d)
reclass_7 = zonal_stat((np.where(reclassify_flatten == 7)), polyfit_slope_2d)
reclass_8 = zonal_stat((np.where(reclassify_flatten == 8)), polyfit_slope_2d)

final_reclass = np.append(reclass_1, reclass_2, reclass_3, reclass_4, 
                          reclass_5, reclass_6, reclass_7, reclass_8)

pd.DataFrame(final_reclass).to_csv('reclass.csv', index_label = "Index", 
                                    header  = ['Count','Minimum','Mean', 
                                              'Maximum', 'Standard Deviation'])


slope_0 = zonal_stat((np.where(slope_reclassify_mask == 0)), polyfit_slope_2d)
slope_1 = zonal_stat((np.where(slope_reclassify_mask == 1)), polyfit_slope_2d)
slope_2 = zonal_stat((np.where(slope_reclassify_mask == 2)), polyfit_slope_2d)
slope_3 = zonal_stat((np.where(slope_reclassify_mask == 3)), polyfit_slope_2d)
slope_4 = zonal_stat((np.where(slope_reclassify_mask == 4)), polyfit_slope_2d)
slope_5 = zonal_stat((np.where(slope_reclassify_mask == 5)), polyfit_slope_2d)
slope_6 = zonal_stat((np.where(slope_reclassify_mask == 6)), polyfit_slope_2d)
slope_7 = zonal_stat((np.where(slope_reclassify_mask == 7)), polyfit_slope_2d)
slope_8 = zonal_stat((np.where(slope_reclassify_mask == 8)), polyfit_slope_2d)
slope_9 = zonal_stat((np.where(slope_reclassify_mask == 9)), polyfit_slope_2d)
slope_10 = zonal_stat((np.where(slope_reclassify_mask == 10)), polyfit_slope_2d)

final_slope = np.append(slope_0, slope_1, slope_2, slope_3, slope_4, slope_5, 
                        slope_6, slope_7, slope_8, slope_9, slope_10)

pd.DataFrame(final_slope).to_csv('slope.csv', index_label = "Index", 
                                  header  = ['Count','Minimum','Mean', 
                                            'Maximum', 'Standard Deviation'])
    
    
    


