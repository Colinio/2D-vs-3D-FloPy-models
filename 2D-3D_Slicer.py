# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 14:03:31 2018

@author: manz_cn
"""

# -*- coding: utf-8 -*-


import pandas as pd
import scipy.interpolate
import os
import utm
import numpy as np
import imod

#choose cross section location that you want:
cs_id = 88284       #cross section name
sp = 10	        #n stress periods
time = 1000			#time per stress period in years
nrow = 266			#3D model rows
ncol = 285			#3d model cols
nlay = 34			#layers
m = sp - 1
csv_files =  r'D:\manz_cn\Desktop\_2D_for_Colin\output_final\riv_sens\baseline\88284'
out_folder = r'D:\manz_cn\Desktop\3D_results\RIV_SENS\3D_model_slices\geo'

#Steps:
#1. import x,y,z,head,conc,vx,vz,vy from 3D model
#2. read tefile 3D and assign coordinates to the values:

#1.
tec_3D = imod.tec.load(r'D:\manz_cn\Desktop\3D_results\Bkk_006_3.2.10_f_new_geology\SP_6-11\CONCVELO.tec', times=[0,-1])	#only read last stress period
#2.
#make 3D model grid with coordinates for x and y position of tecfile units
cell_size = 1000 #in meters
#lower left corner:
X_min = 558342.7
Y_min = 1401600.

#define coordinates of the upper right corner of the model area
X_max = X_min + ncol * cell_size
Y_max = Y_min + nrow * cell_size

#define model area grid:
x_grid = np.arange( X_min,  X_max , cell_size) + 0.5 * cell_size
y_grid = np.arange( Y_max,  Y_min , -cell_size) - 0.5 * cell_size

xx, yy = np.meshgrid( x_grid, y_grid)

#make list of x and y values to be pared up with the concvelo "head" and "conc" outputs.
x_tec = []
for i in range(0,((xx[0,:].size))):
    for j in range(0,((xx[:,0].size))):
        x_tec.append(xx[j,i])

y_tec = []
for i in range(0,((yy[0,:].size))):
    for j in range(0,((yy[:,0].size))):
        y_tec.append(yy[j,i])

#read cross section files to obtain dimensions (amount of columns), needed for writing the tecfile header!  
df = pd.read_csv(os.path.join(csv_files, str(cs_id)+ '_SP_' +str(10) + '_layer_'+str(1) + ".csv"),delimiter=' ',
                     names = ["X", "Y", "Z", "HEAD", "CONC"]) #add vx, vy, vz

#get coordinates from cross section: x,y and z per layer needed
x = df.iloc[:,0]
x = np.array(x[::10])


#insert the stressperiods you want to make slices of to be k (in this case I only want the 5th stress period so I insert 4 which is the 5th count due to zero based indexing)
for k in range(m, sp):
    file = open(os.path.join(out_folder, '3D2D_cs_model_' + str(cs_id)+ "_SP_" +str(k+1) + ".tec"), "w")  #writes tecfile for each SP
    #I = 2370, J = 1, K = 34
    file.write("VARIABLES= \"X\", \"Y\", \"Z\", \"HEAD\", \"CONC\", \"VX\", \"VY\", \"VZ\""        + "\n")
    file.write("ZONE T = \""+str(cs_id))
    file.write(" ")
    file.write(str((k+1)*time))
    file.write(" YEARS\""" I = ")
    file.write(str(x.size))
    file.write(", J = ")
    file.write(str(1))
    file.write(", K = ")
    file.write(str(nlay))
    file.write(", DATAPACKING = POINT"+ "\n")

    for l in range(0,nlay):
        df = pd.read_csv(os.path.join(csv_files, str(cs_id)+ '_SP_' +str(k+1) + '_layer_'+str(l+1) + ".csv"),delimiter=' ',
                             names = ["X", "Y", "Z", "HEAD", "CONC", "VX", "VZ"]) #to do: add vx, vz

        #get coordinates from cross section: x,y and z per layer needed
        x = df.iloc[:,0]
        y = df.iloc[:,1]
        z = df.iloc[:,2]
       
        x = np.array(x[::10])
        y = np.array(y[::10])
        z = np.array(z[::10])
        
        utm_lst = []
        
        x_utm = x
        y_utm = y

        conc_tec_array = np.array(tec_3D.conc[-1,:,:,:])
        conc_tec = []
        for i in range(0,((yy[0,:].size))):
            for j in range(0,((yy[:,0].size))):
                conc_tec.append(conc_tec_array[l,j,i])            
        conc_tec = np.array(conc_tec)

        head_tec_array = np.array(tec_3D.head[-1,:,:,:])
        head_tec = []
        for i in range(0,((yy[0,:].size))):
            for j in range(0,((yy[:,0].size))):
                head_tec.append(head_tec_array[l,j,i])  
        head_tec = np.array(head_tec)
        
        vx_tec_array = np.array(tec_3D.vx[-1,:,:,:])
        vx_tec = []
        for i in range(0,((yy[0,:].size))):
            for j in range(0,((yy[:,0].size))):
                vx_tec.append(vx_tec_array[l,j,i])           
        vx_tec = np.array(vx_tec)       
        
        vy_tec_array = np.array(tec_3D.vy[-1,:,:,:])
        vy_tec = []
        for i in range(0,((yy[0,:].size))):
            for j in range(0,((yy[:,0].size))):
                vy_tec.append(vy_tec_array[l,j,i])  
        vy_tec = np.array(vy_tec)
         
        vz_tec_array = np.array(tec_3D.vz[-1,:,:,:])
        vz_tec = []
        for i in range(0,((yy[0,:].size))):
            for j in range(0,((yy[:,0].size))):
                vz_tec.append(vz_tec_array[l,j,i])  
        vz_tec = np.array(vz_tec)
        
        
        #assign nearest value from 3D grid to the 2D grid locations:
        head_3D_2D = scipy.interpolate.griddata(points=(y_tec, x_tec), values = head_tec, xi=(y_utm, x_utm), method="nearest")
        conc_3D_2D = scipy.interpolate.griddata(points=(y_tec, x_tec), values = conc_tec, xi=(y_utm, x_utm), method="nearest")
        
        vx_3D_2D = scipy.interpolate.griddata(points=(y_tec, x_tec), values = vx_tec, xi=(y_utm, x_utm), method="nearest")
        vy_3D_2D = scipy.interpolate.griddata(points=(y_tec, x_tec), values = vy_tec, xi=(y_utm, x_utm), method="nearest")
        vz_3D_2D = scipy.interpolate.griddata(points=(y_tec, x_tec), values = vz_tec, xi=(y_utm, x_utm), method="nearest")
        
		#for each stress period right results in tecfile
        file = open(os.path.join(out_folder, '3D2D_cs_model_' + str(cs_id)+ "_SP_" +str(k+1) + ".tec"), mode = "a")
        for j in range (0, x_utm.size):
                file.write(str(x_utm[j]))
                file.write(str(","))
                file.write(str(y_utm[j]))
                file.write(str(","))
                file.write(str(z[j]))
                file.write(str(","))
                file.write(str(head_3D_2D[j]))
                file.write(str(","))
                file.write(str(conc_3D_2D[j])) 
                file.write(str(","))
                file.write(str(vx_3D_2D[j]))
                file.write(str(","))
                file.write(str(vy_3D_2D[j]))
                file.write(str(","))
                file.write(str(vz_3D_2D[j])+ "\n")#\n for next line!
        
    file.write("TEXT X=6 Y=83 T=\"Time=")
    file.write(str((k+1)*time))
    file.write("years\", F=TIMES-BOLD, CS=FRAME, H=3, ZN= ")
    file.write(str(k+1)+"\n")
    file.close()
