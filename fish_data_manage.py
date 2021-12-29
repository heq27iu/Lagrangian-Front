# -*- coding: utf-8 -*-
"""
Created on Wed Nov 24 21:00:04 2021

@author: 6
"""

from mpl_toolkits.basemap import Basemap
map = Basemap(llcrnrlon = 0, llcrnrlat = -90, urcrnrlon = 360, urcrnrlat = 90,resolution = 'l')
#map = Basemap(llcrnrlon = 100, llcrnrlat = 0, urcrnrlon = 180, urcrnrlat = 60,resolution = 'l')
import matplotlib as mpl
import matplotlib.pyplot as plt 
from scipy import interpolate
import math
from scipy import misc
import netCDF4 as nckl
import time
import sys
import numpy as np
import os
import cmocean
from osgeo import gdal
import pandas as pd
from scipy.interpolate import spline
from multiprocessing import Pool
from scipy.ndimage import filters
from scipy.ndimage import morphology



if __name__=="__main__": 
    
    print ("开始",time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))

    a = []
    b = []
    c = []
    d = []
    e = []
    for year in range (2018,2019,1):
        fish_map = [[0 for x in range(360)]for x in range(180)]
        fish_map = np.array(fish_map, dtype=np.float)
        fish_map_fre = [[0 for x in range(360)]for x in range(180)]
        fish_map_fre = np.array(fish_map, dtype=np.float)
        for i in range (0,125362,1):
            fish_line = fish[i]
            if(str(fish_line[0])==str(year)):
                lat  = fish_line[2]
                lat_symbel = lat[len(lat)-1]
                if (lat_symbel == "N"):
                    s= str(lat)
                    s = s.rstrip('N')
                    lat = int(s)
                else:
                    s= str(lat)
                    s = s.rstrip('S')
                    lat = -int(s)
                    
                lon  = fish_line[3]
                lon_symbel = lon[len(lon)-1]
                if (lon_symbel == "E"):
                    s= str(lon)
                    s = s.rstrip('E')
                    lon = int(s)
                else:
                    s= str(lon)
                    s = s.rstrip('W')
                    lon = -int(s)+360
                specie_data_pos = 12 
                fish_map[lat+90][lon] = fish_map[lat+90][lon] + float(fish_line[specie_data_pos])+float(fish_line[specie_data_pos+4])+float(fish_line[specie_data_pos+8])+float(fish_line[specie_data_pos+12])+float(fish_line[specie_data_pos+16])
                fish_map_fre[lat+90][lon] = fish_map_fre[lat+90][lon] + 1

        frequence_draw_0 = eval(str('Fre_'+str(year)))
        inten_draw_0 = eval(str('draw_'+str(year)))
        frequence_draw = [[0 for x in range(360)]for x in range(180)]
        frequence_draw = np.array(frequence_draw, dtype=np.float)
        inten_draw = [[0 for x in range(360)]for x in range(180)]
        inten_draw = np.array(inten_draw, dtype=np.float)
        for i in range (0,180):
            for j in range (0,360):
                frequence_draw[i][j] = frequence_draw_0[i*4][j*4]+frequence_draw_0[i*4+1][j*4]+frequence_draw_0[i*4+2][j*4]+frequence_draw_0[i*4+3][j*4]+frequence_draw_0[i*4][j*4+1]+frequence_draw_0[i*4+1][j*4+1]+frequence_draw_0[i*4+2][j*4+1]+frequence_draw_0[i*4+3][j*4+1]
                frequence_draw[i][j] = frequence_draw[i][j]+frequence_draw_0[i*4][j*4+2]+frequence_draw_0[i*4+1][j*4+2]+frequence_draw_0[i*4+2][j*4+2]+frequence_draw_0[i*4+3][j*4+2]+frequence_draw_0[i*4][j*4+3]+frequence_draw_0[i*4+1][j*4+3]+frequence_draw_0[i*4+2][j*4+3]+frequence_draw_0[i*4+3][j*4+3]
                inten_draw[i][j] = inten_draw_0[i*4][j*4]+inten_draw_0[i*4+1][j*4]+inten_draw_0[i*4][j*4+1]+inten_draw_0[i*4+1][j*4+1]
                
        for i in range (60,120):
            for j in range (120,200):
                if(fish_map[i][j]>1500 and fish_map[i][j]<10000):
                    a.append(fish_map[i][j])
                    b.append(frequence_draw[i][j])
                    c.append(inten_draw[i][j])
                    
#                    aa = []
#                    bb = []
#                    cc = []
#                    
#                    for ii in range (-1,2,1):
#                        for jj in range (-1,2,1):
#                            if(fish_map[i+ii][j+11]>0):
#                                aa.append(fish_map[i+ii][j+jj])
#                                bb.append(frequence_draw[i+ii][j+jj])
#                                cc.append(inten_draw[i+ii][j+jj])                                
#
#                    if (len(aa)>4):
#                        co_ab = np.corrcoef(aa,bb)
#                        co_ac = np.corrcoef(aa,cc)
#                        d.append(co_ab[0][1])
#                        e.append(co_ac[0][1])
    

    a = np.array(a, dtype=np.float)
    b = np.array(b, dtype=np.float)
    c = np.array(c, dtype=np.float)        
    ab = np.corrcoef(a,b)
    ac = np.corrcoef(a,c)
#    print(np.mean(d))
#    print(np.mean(e))
    print(ab)
    print(ac)
#    fish_corr[1][0] = ab[0][1]
#    fish_corr[9][8] = ac[0][1]