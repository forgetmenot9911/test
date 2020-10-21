# -*- coding: utf-8 -*-
import cmaps
import netCDF4 as nc
import numpy as np
import pandas as pd
import cartopy
import cartopy.crs as ccrs
import cartopy.feature as cfeat
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from cartopy.io.shapereader import Reader
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import warnings


plt.rcParams['font.sans-serif'] = ['Times New Roman'] 
plt.rcParams['axes.unicode_minus'] = False 


# 数据提取
filename='NC_H08_20200915_1530_L2ARP021_FLDK.02401_02401.nc'
f = nc.Dataset(filename)
print(f.variables.keys(),'\n')#显示所有变量
print(f.variables['AOT'])


lons = f.variables['longitude'][:]
lats = f.variables['latitude'][:]
AOT_arr = f.variables["AOT"]

# 这个循环将所有Nodata的值（即-32768）全部改为0
for i in range(len(AOT_arr)):
    for j in range(len(AOT_arr[0])):
        if AOT_arr[i][j] == -32768:
            AOT_arr[i][j] = 0.0

# AOT_arr = np.array([AOT_arr])
print(AOT_arr)

# 绘制底图
def create_map(extent,xstep,ystep):
    # 加载shp
    shp_path_p = '/Users/wangyucheng/Desktop/kuihua/Province_9/'
    shp_path_c = '/Users/wangyucheng/Desktop/kuihua/diqujie/'
    # --创建画图空间
    proj = ccrs.PlateCarree()  # 创建坐标系
    fig = plt.figure(figsize=(8,10))  # 创建页面
    ax = plt.axes(projection=ccrs.PlateCarree())

    # --设置地图属性
    reader_p = Reader(shp_path_p  +  'Province_9.shp')
    province = cfeat.ShapelyFeature(reader_p.geometries(),proj, edgecolor='k',facecolor='none')
    reader_c = Reader(shp_path_c  +  'diquJie_polyline.shp')
    city = cfeat.ShapelyFeature(reader_c.geometries(),proj, edgecolor='k',facecolor='none')
    # 加载省界线
    ax.add_feature(province)
    # 加载市界
    ax.add_feature(city)
    # 加载经纬范围
    ax.set_extent(extent, crs=proj)
    # --设置网格点属性
    gl = ax.gridlines(
        crs=ccrs.PlateCarree(),
        draw_labels=True,
        linewidth=1.2,
        color='k',
        alpha=0.5,
        linestyle='--'
    )
    gl.xlabels_top = False  # 关闭顶端的经纬度标签
    gl.ylabels_right = False  # 关闭右侧的经纬度标签
    gl.xformatter = LONGITUDE_FORMATTER  # x轴设为经度的格式
    gl.yformatter = LATITUDE_FORMATTER  # y轴设为纬度的格式
    gl.xlocator = mticker.FixedLocator(np.arange(extent[0], extent[1] + xstep, xstep))
    gl.ylocator = mticker.FixedLocator(np.arange(extent[2], extent[3] + ystep, ystep))
    return fig,ax
# 加载经纬范围：南京
extent=[118, 119, 31.6, 32.6]
fig,ax = create_map(extent=extent,xstep=0.2,ystep=0.2)

# 标记nuist
nuist_lon,nuist_lat = 118.72154820009,32.209721542567
plt.plot(nuist_lon,nuist_lat,
        color='blue', linewidth=2, marker='o',
         transform=ccrs.Geodetic())
plt.text(nuist_lon + 0.006, nuist_lat, u"南信大",
         horizontalalignment='left',
         transform=ccrs.Geodetic(),
         fontdict={'family':'Heiti TC'})

# 绘制等值线图
con = plt.contourf(lons,lats,AOT_arr,cmap=cmaps.BlueWhiteOrangeRed)
plt.contour(con,colors='k')
plt.colorbar(con,orientation='horizontal')


plt.show()
