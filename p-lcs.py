# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 10:12:19 2021

@author: 6
"""

# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 11:01:14 2020

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
import netCDF4 as nc
import time
import sys
import numpy as np
import os
from osgeo import gdal
from multiprocessing import Pool
from scipy.ndimage import filters
from datetime import datetime
from scipy.ndimage import morphology

#from scipy.ndimage.morphology import generate_binary_structure, binary_erosion

def load_value(total,time_0):
#        final = []
        time_str = ('H:/g1999/result/'+str(time_0)+'_total')
        f = open(time_str+".txt", "r")  
        while True: 
            line = f.readline()  
            if line:  
                lines = line.split(' ')
                long = len(lines)
                C = []
                for i in range(0,long-1,1):           
                    C.append(float(lines[i]))
                total.append(C)        
            else:  
                break
        f.close()
        
def load_height():
    str_day = ('E:/ETOPO1_Bed_gdal')
    file = nc.Dataset(str_day+".grd")
    data = file['z'][:]
    height = [[0 for x in range(3600)]for x in range(1800)]
    height = np.array(height, dtype=np.float)
    for i in range (0,1800,1):
        for j in range (0,1800,1):
            height[i][j] = data[(i*6)*21601+(j+1800)*6]
    for i in range (0,1800,1):
        for j in range (1800,3600,1):
            height[i][j] = data[(i*6)*21601+(j-1800)*6]
    return height

def EigenCompute(a,b,c):
     mEntry = [0 for x in range(0,3)];
     eigVec = [0 for x in range(0,4)];
     eigVal = [0 for x in range(0,2)];
     ##待求2*2矩阵为mEntry[0]，mEntry[1]，mEntry[1]，mEntry[2]
     mEntry[0] = a
     mEntry[1] = b
     mEntry[2] = c
     sum = math.fabs(mEntry[0]) + math.fabs(mEntry[2])
     ##The matrix M is diagonal (within numerical round-off).
     if (math.fabs(mEntry[1]) + sum == sum):
         eigVec[0] = 1;
         eigVec[1] = 0;
         eigVec[2] = 0;
         eigVec[3] = 1;
         eigVal[0] = mEntry[0];
         eigVal[1] = mEntry[2];
     else:
         trace = mEntry[0] + mEntry[2];
         diff = mEntry[0] - mEntry[2];
         discr = math.sqrt(diff*diff + 4*mEntry[1]*mEntry[1]);
         eigVal0 = 0.5*(trace - discr);
         eigVal1 = 0.5*(trace + discr);
         eigVal[0] = eigVal0;
         eigVal[1] = eigVal1;
         if (diff >= 0):
             cs = mEntry[1];
             sn = eigVal0 - mEntry[0];
         else:
             cs = eigVal0 - mEntry[2];
             sn = mEntry[1];
  
         tempval=cs*cs + sn*sn;
         invLength = 0;
         if (tempval > 0.0):
             invLength= 1.0/(math.sqrt(tempval));
         else:
             print('Division by zero in InvSqr\n')

         cs *= invLength;
         sn *= invLength;

         eigVec[0] = cs;
         eigVec[1] = sn;
         eigVec[2] = -sn;
         eigVec[3] = cs;
         
     return (eigVec[0] , eigVec[1] , eigVec[2] , eigVec[3] , eigVal[0] , eigVal[1])

def cpoint_type(cpoint_a,cpoint_b,new1,new2):
    cpoint_1 = new1(cpoint_a-0.1,cpoint_b)
    cpoint_2 = new2(cpoint_a-0.1,cpoint_b)
    cpoint_3 = new1(cpoint_a+0.1,cpoint_b)
    cpoint_4 = new2(cpoint_a+0.1,cpoint_b)
    cpoint_5 = new1(cpoint_a,cpoint_b-0.1)
    cpoint_6 = new2(cpoint_a,cpoint_b-0.1)
    cpoint_7 = new1(cpoint_a,cpoint_b+0.1)
    cpoint_8 = new2(cpoint_a,cpoint_b+0.1)
    drgree = []
    drgree.append(math.atan2(cpoint_1,cpoint_2))
    drgree.append(math.atan2(cpoint_3,cpoint_4))
    drgree.append(math.atan2(cpoint_5,cpoint_6))
    drgree.append(math.atan2(cpoint_7,cpoint_8))
    drgree.sort()
    edge = 4.19 #4.19
    if((drgree[1]-drgree[0])>edge or (drgree[2]-drgree[1])>edge or (drgree[3]-drgree[2])>edge or (6.28-drgree[3]+drgree[0])>edge):
        point_type = 2#楔形
    else:
        point_type = 1#三角形
    return point_type

def draw ():        
    print ("开始绘制hlcs",time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))          
    fig = plt.figure(figsize=(15,10))
    map.fillcontinents(color=[181/255,180/255,181/255],zorder = 3)
#    cpointlon = cpoint[:,0]
#    cpointlat = cpoint[:,1]
#    plt.scatter(cpointlon, cpointlat, s=2, label = '$dislike$', c = 'red', marker='.', alpha = None)    
    for i in range(0,len(line_all),1):
        line_all_a = np.array(line_all[i], dtype=np.float)
        plt.scatter(line_all_a[0,0], line_all_a[0,1], s=2, label = '$dislike$', c = 'red', marker='.', alpha = None)    
        x = line_all_a[:,0]
        y = line_all_a[:,1]
        if(cpoint[i][3]==1):
            plt.plot(x,y, c = (0,0,0), lw=0.2, zorder = 2)
        else:
            plt.plot(x,y, c = (0,0,1), lw=0.2, zorder = 2)
            
    parallels = np.arange(0.,81,10.)
    map.drawparallels(parallels,labels=[False,True,True,False],zorder = 8)
    meridians = np.arange(20.,351.,20.)
    map.drawmeridians(meridians,labels=[True,False,False,True],zorder = 8)
    plt.savefig('p2', dpi=512)   
    

if __name__=="__main__": 
    for time_0 in range (1,13,1): 
        print (time_0,"开始读取数据",time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
        total = []
        load_value(total,time_0)
        height = load_height()
        c11 = [[0 for x in range(5760)]for x in range(2880)]
        c12 = [[0 for x in range(5760)]for x in range(2880)]
        c22 = [[0 for x in range(5760)]for x in range(2880)]
        c11 = np.array(c11, dtype=np.float)
        c12 = np.array(c12, dtype=np.float)
        c22 = np.array(c22, dtype=np.float)
        value = [[[0 for x in range(6)] for x in range(5760)]for x in range(2880)]
        value_mark = [[0 for x in range(18000)]for x in range(71999)]
        value = np.array(value, dtype=np.float)
        value_mark = np.array(value_mark, dtype=np.float)
        print (time_0,"开始计算特征值",time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
        for i in range (0,2880,1):
            for j in range (0,5760,1):
                c11[i][j] = total[i*5760+j][0]
                c12[i][j] = total[i*5760+j][1]
                c22[i][j] = total[i*5760+j][2]
                (value[i,j]) = EigenCompute (c11[i,j],c12[i,j],c22[i,j])
    #    for i in range(0,800,1):
    #        for j in range(0,1280,1):
    #            (value[i,j]) = EigenCompute (c11[i,j],c12[i,j],c22[i,j])
        print (time_0,"开始计算奇点",time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))        
    #    mark
        cpoint = []
        X = np.linspace(0,360,5760) 
        Y = np.linspace(-90,90,2880)    
        c11_new = interpolate.interp2d(X,Y,c11, kind='linear')
        c12_new = interpolate.interp2d(X,Y,c12, kind='linear')
        c22_new = interpolate.interp2d(X,Y,c22, kind='linear')
        new1 = interpolate.interp2d(X,Y,value[:,:,0], kind='linear')
        new2 = interpolate.interp2d(X,Y,value[:,:,1], kind='linear')
        new3 = interpolate.interp2d(X,Y,value[:,:,2], kind='linear')
        new4 = interpolate.interp2d(X,Y,value[:,:,3], kind='linear')
        new5 = interpolate.interp2d(X,Y,value[:,:,4], kind='linear')
        height_lon = np.linspace(0,360,3600) 
        height_lat = np.linspace(90,-90,1800) 
        height_new = interpolate.interp2d(height_lon,height_lat,height, kind='linear')
        arr = value[:,:,4]
        
        
        
#        start_time=datetime.now()       
#        threadsperblock = 32
#        blockspergrid = (xx.size + (threadsperblock - 1)) // threadsperblock
#        GPU_C[blockspergrid, threadsperblock](xt,yt,zt,deltax,deltay,deltaz,sVec,Norma,k1,k2,k3,xx,yy,zz,phiPrGr_0101,phiPrGr_0202,phiPrGr_0303,phiPrGr_0404,phiPrGr_0505,phiPrGr_0606,C22mC11Gr_00,C12Gr_00,C_all_dict,X,Y,dis_x,dis_y)
#        print(datetime.now()-start_time) # 用时
        
        #求奇点
        for i in range(0,2870,1):
            for j in range(0,5750,1):
    #    for i in range(1440,2400,1):
    #        for j in range(1600,2880,1):
                if(c12[i][j]*c12[i+1][j]<0 or c12[i][j]*c12[i][j+1]<0 or c12[i+1][j]*c12[i+1][j+1]<0 or c12[i][j+1]*c12[i+1][j+1]<0):
                    if((c11[i][j]-c22[i][j])*(c11[i][j+1]-c22[i][j+1])<0 or (c11[i][j]-c22[i][j])*(c11[i+1][j]-c22[i+1][j])<0 or (c11[i+1][j]-c22[i+1][j])*(c11[i+1][j+1]-c22[i+1][j+1])<0 or (c11[i][j+1]-c22[i][j+1])*(c11[i+1][j+1]-c22[i+1][j+1])<0):
                        heigh = height_new(j/16,i/16-90)
                        if (heigh<-50):
                            flag = 0
                            for ii in range (0,10,1):
                                for jj in range (0,10,1):
                                    aabs = abs(c12_new(j/16+jj/160,i/16+ii/160-90))+abs(c11_new(j/16+jj/160,i/16+ii/160-90)-c22_new(j/16+jj/160,i/16+ii/160-90))
                                    if (aabs<0.1):
                                        pointt = []
                                        pointt.append(j/16+jj/160)
                                        pointt.append(i/16+ii/160-90)
                                        valuee = new5(j/16+jj/160,i/16+ii/160-90)
                                        pointt.append(valuee)
                                        flag = 1
                                        cpp_lon = int((j/16+jj/160)*50)
                                        cpp_lat = int((i/16+ii/160)*50)
                                        point_type_0 = cpoint_type(pointt[0],pointt[1],new1,new2)
                                        pointt.append(point_type_0)
                                        cpoint.append(pointt)
                                        for buffa in range (-5,5,1):
                                            for buffb in range (-5,5,1):
                                                value_mark[cpp_lat+buffa][cpp_lon+buffb] = point_type_0
                                        break
                                if(flag == 1):
                                    break
        cpoint = np.array(cpoint, dtype=np.float)
        print(np.mean(cpoint[:,3]))
        time_str = ('H:/g1999/cpoint_'+str(time_0))
        np.save(time_str, cpoint)  
    #    sys.exit()
        del c11,c12,c22,height
        print (time_0,"开始计算plcs",time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
        line_weight = []
    #type == 2
        step = 0.05    
        lone = 300    
        line_all = []
        for i in range(0,len(cpoint),1):
            if (cpoint[i][3]==2):
                lat = cpoint[i][1]
                lon = cpoint[i][0]
                line = []
                flag = 0
                for h in range(0,lone,1): 
                    point = []
                    point.append(lon)
                    point.append(lat)
                    lon_add = (new1(lon,lat)/(new1(lon,lat)**2+new2(lon,lat)**2)**0.5)*step
                    lat_add = (new2(lon,lat)/(new1(lon,lat)**2+new2(lon,lat)**2)**0.5)*step
    #                lon_add = (new3(lon,lat)/(new3(lon,lat)**2+new4(lon,lat)**2)**0.5)*step
    #                lat_add = (new4(lon,lat)/(new3(lon,lat)**2+new4(lon,lat)**2)**0.5)*step
                    lon = lon + lon_add[0]
                    lat = lat + lat_add[0]
                    line.append(point)
                    if (lon>359 or lon<1):
                        break
                    if (lat>89 or lat<-89):
                        break
                    cpp_lon = int(lon*50)
                    cpp_lat = int((lat+90)*50)
                    if(value_mark[cpp_lat][cpp_lon] !=0 and h>15):
                        flag =1
                        break
                if(len(line)>10 and flag==1):
                    line_all.append(line)
                    line_weight_t = new5(line[0][0],line[0][1])
                    line_weight.append(line_weight_t)
                
                
                step = -0.05  
        lone = 300          
        for i in range(0,len(cpoint),1):
            if (cpoint[i][3]==1):
                lat = cpoint[i][1]
                lon = cpoint[i][0]
                line = []
                flag = 0
                for h in range(0,lone,1): 
                    point = []
                    point.append(lon)
                    point.append(lat)
    #                lon_add = (new1(lon,lat)/(new1(lon,lat)**2+new2(lon,lat)**2)**0.5)*step
    #                lat_add = (new2(lon,lat)/(new1(lon,lat)**2+new2(lon,lat)**2)**0.5)*step
                    lon_add = (new3(lon,lat)/(new3(lon,lat)**2+new4(lon,lat)**2)**0.5)*step
                    lat_add = (new4(lon,lat)/(new3(lon,lat)**2+new4(lon,lat)**2)**0.5)*step
                    lon = lon + lon_add[0]
                    lat = lat + lat_add[0]
                    line.append(point)
                    if (lon>359 or lon<1):
                        break
                    if (lat>89 or lat<-89):
                        break
                    cpp_lon = int(lon*50)
                    cpp_lat = int((lat+90)*50)
                    if(value_mark[cpp_lat][cpp_lon] !=0 and h>15):
                        flag =1
                        break
                if(len(line)>10 and flag==1):
                    line_all.append(line)
                    line_weight_t = new5(line[0][0],line[0][1])
                    line_weight.append(line_weight_t)
        
        time_str = ('H:/g1999/line_'+str(time_0))
        np.save(time_str, line_all)    
        weight_str = ('H:/g1999/Lweight_'+str(time_0))
        np.save(weight_str, line_weight)    
        value_str = ('H:/g1999/value_'+str(time_0))
        np.save(value_str, arr)    
        draw ()
        print (time_0,"结束",time.strftime('%Y-%m-%d %H:%M:%S',time.localtime(time.time())))
    
