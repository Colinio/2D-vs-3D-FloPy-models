# -*- coding: utf-8 -*-
"""
Created on Mon July 26 13:02:19 2018

@author: colin
"""

import matplotlib.pyplot as plt
import numpy as np
#from modeltools import cs_model 
import os
from matplotlib.colors import Normalize
from matplotlib import colors
from matplotlib import colorbar
import matplotlib.cm as cmx
import matplotlib.patches as mpatches
import matplotlib 
import collections
from time import time
import flopy
import itertools
import flopy.utils.binaryfile as bf
import pandas as pd
import geopandas as gpd
import csv
import utm

#delete old csv files:

#   plotting at times
plot_profile_yrs_all = [0.1, 10, 15, 30, 50, 75, 80, 100, 150, 200, 250, 300, 350, 400, 450, 500, 750, 1000]#, 2000, 3000, 4000, 5000, 7500, 10000, 100010, 100020, 10030, 10040, 10050, 10060, 10070, 10080, 10090, 10100, 10125, 10150, 10175, 10200]

#   Define all important directories 

model_dir = r'\models'
 
mf_exe_dir = r'\mf2005.exe'
mt3d_exe_dir = r'\MT3D-USGS_64.exe'
swat_exe_dir = r'\swt_v4x64.exe'

model_dir = r'\model_dir'
output_dir = r'\output_dir'

#   loop through the model_dir and go COSCAT by COSCAT region
for subdir, dirs, files in os.walk(model_dir):
    for file in files:
        if file.endswith("_top_bot_elev.npy"):
            full_dir = os.path.join(subdir, file)

            #   split the directory string to get the COSCAT number and type of coast
            cs_id = int(full_dir.split('_model_')[1].split('\\')[0])
            print('Building model for cross-section ID : ' + str(cs_id))
            
            #   define the output directory
            out_dir = os.path.join(model_dir + '\\cs_model_' + str(cs_id), '_model_files')
            modelname = 'cs_%s' % (cs_id)

            #   1) create the models via flopy
            mf = flopy.modflow.Modflow(modelname, exe_name = mf_exe_dir, model_ws = out_dir)
            mt = flopy.mt3d.Mt3dms(modelname, 'nam_mt3dms', mf, model_ws = out_dir, exe_name = mt3d_exe_dir)
            mswt = flopy.seawat.Seawat(modelname, 'nam_swt', mf, mt, model_ws = out_dir, exe_name = swat_exe_dir)

            #   read in the numpy dictionary with the model IBOUND and other files and start building the model
            dict_top_bot_in = np.load(full_dir).item()      
            dict_kh_in = np.load(os.path.join(model_dir + '\\cs_model_' + str(cs_id), '_kh_vals.npy')).item()
            dict_ghb_in = np.load(os.path.join(model_dir + '\\cs_model_' + str(cs_id), '_ghb_vals.npy')).item()
            dict_riv_in = np.load(os.path.join(model_dir + '\\cs_model_' + str(cs_id), '_riv_vals.npy')).item()
            dict_rch_in = np.load(os.path.join(model_dir + '\\cs_model_' + str(cs_id), '_rch_vals.npy')).item()
            dict_bnd_in = np.load(os.path.join(model_dir + '\\cs_model_' + str(cs_id), '_bnd_vals.npy')).item()
            dict_xy_coord_in = np.load(os.path.join(model_dir + '\\cs_model_' + str(cs_id), '_xy_coord.npy')).item()
            
            #   the top elevation list is the one with dictionary key 0
            top_elev = dict_top_bot_in['0']            
            
            #   read in the IBOUND values for the cross-section
            ibound_lst = [int(x) for x in dict_bnd_in['BND_ACTIVE']]
            
            """
            ['HEAD_L10_40M',
             'COND_L10_40M_SEA_10K_LAND_0.1',
             'CONC_L10_40M',
             'DENS_L10_40M',
             'COND_L10_40M_SEA_100_LAND_0.1',
             'COND_L10_40M_SEA_10K_LAND_1',
             'COND_L1_40M_SEA_100_LAND_0.1',
             'COND_L1_40M_SEA_10K_LAND_0.1',
             'HEAD_L1_40M',
             'CONC_L1_40M',
             'DENS_L1_40M',
             'COND_L1_40M_SEA_10K_LAND_1'
             'COND_L1_40M_SEA_100k_LAND_0.1'
             'COND_L1_40M_SEA_1M_LAND_0.1
             'COND_L1_40M_SEA_100K_LAND_0.01']
             
            if you want to change to different boundary condition just change the name below in the ghb_lst = dict_ghb_in['INSERT HERE']             
            """
            
            #   do the same for GHB, RIV, RCH packages
            ghb_lst = dict_ghb_in['COND_L1_40M_SEA_1M_LAND_0.1']
            riv_bot_lst = dict_riv_in['RBOT']
            riv_cond_lst = dict_riv_in['COND']
            riv_stage_lst = dict_riv_in['STAGE']
            rch_lst = dict_rch_in['RCH']
            xy_coord_lst = dict_xy_coord_in['xy_wgs84']
            
            """
            #   create an object of the Model class (see modeltools - dont look there it is too much :D). Basically what this does is that
            #   it creates a model object that has quite some methods/functions/attributes attached to it. 
            model = cs_model(cs_id, top_elev, None, None, None, None, None, None, None,\
                             None, None, None, None, None, None, None, None, None, None, None, None, None,\
                             None, None, cst_bound = False)
            """
            
            #   start building the model itself, starting with the DIS package             
            #       1) Find the start_idx and end_idx based on the BND information provided                
            #       first find the first element in the ibound list with an active cell (=1), then find the index     
            #       of the lowest top elevation in the top_elev list, this will be the offshore limit of the model   
            #       Another issue is to find the position of the inland part of the cross-section - can either be
            #       east or west (left or right) from the coastal point. Based on that the indexes are selected.
            
            #   check the amount of cells below sea level on both sides from the coastal point
            cells_left = sum(1 for i in top_elev[400:420] if i < 0)
            cells_right = sum(1 for i in top_elev[380:400] if i < 0)

            del_col = 100.
            cs_points_dist = 500.            
            
            #   if there are more cells in the left then invert the lists
            if cells_right > cells_left:
                ghb_lst = ghb_lst[::-1]
                riv_bot_lst = riv_bot_lst[::-1]
                riv_cond_lst = riv_cond_lst[::-1]
                riv_stage_lst = riv_stage_lst[::-1]
                rch_lst = rch_lst[::-1]
                top_elev = top_elev[::-1]
                ibound_lst = ibound_lst[::-1]
                xy_coord_lst = xy_coord_lst[::-1]

                ibound_start_idx = ibound_lst.index(1)
                ibound_end_idx = top_elev.index(min(top_elev))
            
                #   the rest are the bottom of each layer, also read in the Kh values for each layer
                bot_elev, kh_lst = [], []            
                for a in range(1, len(dict_top_bot_in.keys())):
                    #   the line below repeats the same value for each column 
                    lst_bot = list(itertools.chain.from_iterable(itertools.repeat(x, int(500 / del_col)) for x in\
                                   dict_top_bot_in[str(a)][::-1][ibound_start_idx : ibound_end_idx]))
                    bot_elev.append(lst_bot)                
                    lst_kh =  list(itertools.chain.from_iterable(itertools.repeat(x, int(500 / del_col)) for x in\
                                   dict_kh_in[str(a)][::-1][ibound_start_idx : ibound_end_idx]))               
                    kh_lst.append(lst_kh)

            else:
                ibound_start_idx = ibound_lst.index(1)
                ibound_end_idx = top_elev.index(min(top_elev))    
                
                #   the rest are the bottom of each layer, also read in the Kh values for each layer
                bot_elev, kh_lst = [], []            
                for a in range(1, len(dict_top_bot_in.keys())):
                    #   the line below repeats the same value for each column 
                    lst_bot = list(itertools.chain.from_iterable(itertools.repeat(x, int(500 / del_col)) for x in\
                                   dict_top_bot_in[str(a)][ibound_start_idx : ibound_end_idx]))
                    bot_elev.append(lst_bot)                
                    lst_kh =  list(itertools.chain.from_iterable(itertools.repeat(x, int(500 / del_col)) for x in\
                                   dict_kh_in[str(a)][ibound_start_idx : ibound_end_idx]))               
                    kh_lst.append(lst_kh)       
                    
            #       2) Calculate the number of columns necessary for the model domain, based on del_col, keep in
            #       mind that the distance between the individual cross-section points is 500m so the column width
            #       should be something like 50m, 100m, 250m, or 500m...
            #       Define the rest of the input parameters for the DIS package
            tot_len = (ibound_end_idx - ibound_start_idx) * 500       #    total model length in meters      
            ncol = int(tot_len / del_col)                             #    total number of columns
            nlay = len(dict_top_bot_in.keys()) - 1                    #    total number of layers
            nrow = 1                                                  #    constant value, only one row
            perlen_ys = 10                                                 
            perlen = [365.25*perlen_ys]  
            nstp = [1]                             # number of time steps in each stress period
            nper = len(perlen)                                               #    number of stress periods
            delc = 1.                                                 #    distance along the columns (= 1m), these delr, delc are a bit confusing, check the MODFLOW manual if you want   
            delr = [del_col] * ncol                                   #    distance along the row..
            laycbd = 0                                                #    0 indicates no confining bed... 
            botm = np.asarray(bot_elev)                               #    create an array with bottom elevations of each layer
            top = np.asarray([list(itertools.chain.from_iterable(itertools.repeat(x, int(500 / del_col)) for x in\
                  top_elev[ibound_start_idx : ibound_end_idx]))])

            #       3) Create the DIS package itself
            dis = flopy.modflow.ModflowDis(mf, nlay, nrow, ncol, nper, delr, delc, laycbd, top, botm, perlen, nstp)            
            
            #   now build the BAS package that includes the starting head profile and the IBOUND array
            #       4) Create the IBOUND array itself, based on the amount of layer and columns
            ibound_arr = np.ones((nlay, 1, ncol), dtype = np.int32)  
            strt_arr = np.zeros((nlay, 1, ncol), dtype = np.int32)  

            #       5) Write the BAS package
            bas = flopy.modflow.ModflowBas(mf, ibound_arr, strt_arr)

            "   create the X and Y lists for plotting and saving the tecplot files later on     "
            "       first get the coordinates of the start and end point of the cross-section   "
            pt_start = xy_coord_lst[ibound_start_idx]
            pt_end = xy_coord_lst[ibound_end_idx]     
            "       create a list of column coordinates for all columns                         "
            pt_x_lst = np.linspace(pt_start[0], pt_end[0], ibound_arr.shape[-1]) # ibound_arr.shape[-1] = number of columns in the model domain
            pt_y_lst = np.linspace(pt_start[1], pt_end[1], ibound_arr.shape[-1])
            
            #   create the LPF package that holds all the hydrogeological information            
            #       6) Create the necassry HK and VK arrays

            anis_fact = 0.1                                          #    anisotropy factor to multiply the HK_arr
            hk_arr = np.asarray(kh_lst)                               #    define the HK array based on the previously created list                        
            hk_arr[hk_arr < 10000] = 10                     #makint it all sand with hk of 10 m/d
            vk_arr = hk_arr * anis_fact                              #    based on the anisotropy factor create the VK array
            
            #       7) Write the LPF package
            lpf = flopy.modflow.ModflowLpf(mf, laytyp = 0, hk = hk_arr, vka = vk_arr) # check the LAYTYP values!!
            
            #       8) Create the GHB, RCH, RIV and DRN packages - all boundary conditions and top system 
            ghb_inland = 0.1                                           #    inland conductance
            ghb_sea = 1000000.                                            #    offshore conductance
            sea_lvl = .0                                              #    define the sea level elevation (relative to present)
                                   
            #   define and initiate the output arrays and lists, one for the GHB and one for the SSM package 
            ghb_input_lst = []
            ghb_ssmdata = []

            #  create the SSM dictionary where the ssm input will be written to
            itype = flopy.mt3d.Mt3dSsm.itype_dict()

            #generate noise array for conc and sea lvl:
            noise = np.random.normal(-0.001, 0.001, ncol)
            
            for c in xrange(nlay):
                ghb_input_lst.append([c, 0, 0, top[0][0], ghb_inland])
                ghb_ssmdata.append([c, 0, 0, 0.0, itype['GHB']])
                ghb_input_lst.append([c, 0, ncol - 1, sea_lvl, ghb_sea])
                ghb_ssmdata.append([c, 0, ncol - 1, 35.0, itype['GHB']])
            
            #   loop through the top elevation list and wherever the top_elev < sea level assign sea GHB 
            for d in xrange(ncol - 1):
                if top[0][d] < sea_lvl:
                    ghb_input_lst.append([0, 0, d, sea_lvl + noise[d], ghb_sea])
                    ghb_ssmdata.append([0, 0, d, 35.0 + noise[d], itype['GHB']])
    
            #   write the final output dictionary, inlcude each stress period
            ghb_arr_in = {}
            for c in xrange(len(perlen)):
                ghb_arr_in[c] = ghb_input_lst

            #       9) Write the GHB package
            ghb = flopy.modflow.ModflowGhb(mf, ipakcb = 1, stress_period_data = ghb_arr_in)  #   ipakcb - write output in cbc file

            #       10) Set up the RCH package, first define the necessary input parameters and lists
            nrchop = 3                                               #  3 - rchrge applied to the highest active cell in each column
            rch_arr = np.array([[0.0] * 1 * ncol], dtype = np.float32) #    the RCH array itself, first fill with 0.0 m/d
            #irch_arr = np.zeros((1, 1, ncol), dtype=np.float)          #    the concentration array of RCH (only fresh so 0.0)

            #   pre process the recharge input list
            rch_in = list(itertools.chain.from_iterable(itertools.repeat(x, int(500 / del_col)) for x in\
                  rch_lst[ibound_start_idx : ibound_end_idx]))

            #   loop through the input recharge list and fill in the respective values in cells that are above sea level
            for d in xrange(len(rch_in)):
                if top[0][d] >= sea_lvl:
                    rch_arr[0][d] = rch_in[d]
                    
            #       11) Write the RCH package
            rch = flopy.modflow.ModflowRch(mf, nrchop = nrchop, ipakcb = 1, rech = rch_arr)
    
            #       12) Create and write the DRN package, define the drainage level as below surface
            drn_lvl = 0.5                                               #   level to be substracted from surface level
            drn_input_lst = []                                          #   drainage is assigned only to cells that receive recharge - cells with elev above sea level
            
            #   loop through the columns and if the elevation is above sea level then assign the drainage 
            for e in xrange(ncol):
                if top[0][e] >= sea_lvl:
                    #   apply the ghb_inland conductance
                    drn_input_lst.append([0, 0, e, top[0][e] - drn_lvl, ghb_inland])

            #   write the final output dictionary, inlcude each stress period
            drn_arr_in = {}
            for f in xrange(len(perlen)):
                drn_arr_in[f] = drn_input_lst    
    
            #       13) Write the DRN package
            drn = flopy.modflow.ModflowDrn(mf, ipakcb = 1, stress_period_data = drn_arr_in)
        
            #       14) Define the input for the RIV package
            riv_input_lst = []

            #   pre process the river input list
            riv_bot_in = list(itertools.chain.from_iterable(itertools.repeat(x, int(500 / del_col)) for x in\
                  riv_bot_lst[ibound_start_idx : ibound_end_idx]))
            riv_cond_in = list(itertools.chain.from_iterable(itertools.repeat(x, int(500 / del_col)) for x in\
                  riv_cond_lst[ibound_start_idx : ibound_end_idx]))
            riv_stage_in = list(itertools.chain.from_iterable(itertools.repeat(x, int(500 / del_col)) for x in\
                  riv_stage_lst[ibound_start_idx : ibound_end_idx]))            
            
            #   loop through the riv_bot_in list and if there are any cells where the value is different from -9999. then assign a RIV cell
            for g in xrange(len(riv_cond_in)):
                if riv_cond_in[g] == -9999.:
                    pass
                else:
                    riv_input_lst.append([0, 0, g, riv_stage_in[g], riv_cond_in[g], riv_bot_in[g]])

            #   write the final output dictionary, inlcude each stress period
            riv_arr_in = {}
            for h in xrange(len(perlen)):
                riv_arr_in[h] = riv_input_lst        
        
            #   only write the package if the input list is not empty
            #if riv_input_lst != []:
            #    riv = flopy.modflow.ModflowRiv(mf, ipakcb = 1, stress_period_data = riv_arr_in)
        
            #       15) Write the rest of the packages..
        
            #           write the PCG package
            hclose = 1e-4           #   criteria for head convergence
            rclose = 1e+1           #   criteria for concentration convergence (check this..)
            pcg = flopy.modflow.ModflowPcg(mf, hclose = hclose, rclose = rclose)

            #           write the OC package
            ihedfm = 1          # a code for the format in which heads will be printed.
            iddnfm = 0          # a code for the format in which drawdowns will be printed.
            extension = ['oc','hds','ddn','cbc']
            unitnumber = [14, 30, 0, 50]
            th_time_step = 10    #   the frequency of time steps to be saved
            #   create the dictionary that defines how to write the output file
            spd = {(0, 0): ['save head', 'save budget']}
            #for t in range(0, nper):
            #per = t #+ 1
            #   xrange allows to iterate through the list with specified step size - 25
            #   to save space on disk, every nth timestep is saved
            for g in xrange(nstp[0]):
                spd[(0, int(g))] = ['save head', 'save budget']
                #spd[(0, int(g) + 1)] = []

            spd[(0, int(g) + 1)] = ['save head', 'save budget']
            #spd[(0, int(g) - 1)] = ['save head', 'save budget']            
            
            oc = flopy.modflow.ModflowOc(mf, ihedfm, iddnfm, cboufm='(20i5)', stress_period_data = spd)
            
            #           write the BTN package
            #   the nprs parameter defines how many transport steps are going to be exported to the UCN file
            th_time_step = [(nstp[0])]           
            dt0 = 0
            porosity = 0.3
            ifmtcn = 0
            nprs = 1
            chkmas = False
            nprmas = 10
            nprobs = 10

            timprs_sp1 = np.linspace(1., perlen[0], th_time_step[0], endpoint = False)
            if nper > 1:
                timprs_sp2 = np.linspace(perlen[0], perlen[0] + perlen[1], th_time_step[1], endpoint = True)
                timprs = np.concatenate((timprs_sp1, timprs_sp2[1:]), axis = 0)
            else:
                timprs = timprs_sp1
            
            #   define the starting concnetration profile (start all saline)
            sconc_arr = ibound_arr * 0.000001
            
            btn = flopy.mt3d.Mt3dBtn(mt, nprs = nprs, timprs = timprs, prsity = porosity, sconc = sconc_arr, ifmtcn = ifmtcn,
                                 chkmas = chkmas, nprobs = nprobs, nprmas = nprmas, dt0 = dt0)
            
            #   12) write the ADV package
            mixelm = 0
            
            adv = flopy.mt3d.Mt3dAdv(mt, mixelm = mixelm)

            #   13) write the DSP package
            dmcoef = 0.0000864    # effective molecular diffusion coefficient [M2/D]
            al = 1.
            trpt = 0.1
            trpv = 0.1
            
            dsp = flopy.mt3d.Mt3dDsp(mt, al = al, trpt = trpt, trpv = trpv, dmcoef = dmcoef)

            #   14) write the GCG package
            iter1 = 500
            mxiter = 1
            isolve = 1
            cclose = 1e-7
            
            gcg = flopy.mt3d.Mt3dGcg(mt, iter1 = iter1, mxiter = mxiter, isolve = isolve, cclose = cclose)

            #   15) write the VDF package
            iwtable = 0
            densemin = 1000.
            densemax = 1025.
            denseref = 1000.
            denseslp = 0.7143
            firstdt = 1e-3

            vdf = flopy.seawat.SeawatVdf(mswt, iwtable = iwtable, densemin = densemin, densemax = densemax,\
                                              denseref = denseref, denseslp = denseslp, firstdt = firstdt)
            #   16) write the SSM package
            ssm_rch_in = np.array([[0.0] * 1 * ncol], dtype = np.float32)
            ssmdata_dict = {0: ghb_ssmdata,\
                            1: ghb_ssmdata}

            ssm = flopy.mt3d.Mt3dSsm(mt, crch = ssm_rch_in, stress_period_data = ssmdata_dict)

            #   17) write all the input files
            mf.write_input()
            mt.write_input()
            mswt.write_input()            
       
            #   18) Run the model and measure the time
            t0 = time()
            #   run the model
            v = mswt.run_model(silent = False, report = True)
            for idx in range(-3, 0):
                print(v[1][idx])
            #   stop measuring time and calculate total run time
            t1 = time()
            run_time = t1 - t0
    
            #  19) Read model output; heads, concentrations and cell budget flow
            ml_results = flopy.modflow.Modflow.load(modelname + ".nam", model_ws = out_dir, verbose = False, check = False, exe_name = "mfnwt")
            hdsobj = bf.HeadFile(os.path.join(out_dir, modelname + '.hds'), model = ml_results)
            head = hdsobj.get_alldata()
            ucnobj = bf.UcnFile(os.path.join(out_dir, 'MT3D001.UCN'), model = ml_results)
            time_steps = ucnobj.get_times()
            conc = ucnobj.get_alldata()
            cbbobj = bf.CellBudgetFile(os.path.join(out_dir, modelname + '.cbc'))
            times_heads = cbbobj.get_times()       
            
            utm_lst = []
            
            #for i in x:
            #transform lat-lon to utm:
            for i in range(0,pt_x_lst.size):
                utm_lst.append(utm.from_latlon(float(pt_y_lst[i]), float(pt_x_lst[i])))
                
                utm_array = np.array(utm_lst)
                utm_array.shape
                
                x_utm = utm_array[:,0]      #check if this is really x or not y
                
                y_utm = utm_array[:,1]
            
                for k in range(0, nstp[0]):

                file = open(os.path.join(output_dir , "cs_model_" + str(cs_id)+ "_SP_" +str(k+1) + ".tec"), "w")  #writes tecfile for each SP

                """Example header:
                VARIABLES= "X", "Y", "Z", "HEAD" , "CONC" , "VC" , "VR", "VL"
                ZONE T = " 168.00 DAYS" I = 51, J = 18, K = 84, DATAPACKING = POINT
                """
                timezn = (perlen[0]/ (nstp[0])) *(k+1)
                file.write("VARIABLES= \"X\", \"Y\", \"Z\", \"HEAD\", \"CONC\"" + "\n")
                file.write("ZONE T = \""+str(cs_id))
                file.write(" ")
                file.write(str(timezn))
                file.write(" DAYS\""" I = ")
                file.write(str(ncol))
                file.write(", J = ")
                file.write(str(nrow))
                file.write(", K = ")
                file.write(str(nlay))
                file.write(", DATAPACKING = POINT"+ "\n")
               
            
                #then write the other stuff for each layer:
                for i in range (0, nlay):                                   #for each layer
                    for j in range(0,(len(pt_x_lst))):                      #for each element in the layer
                        file.write(str(x_utm[j]))
                        file.write(str(","))
                        file.write(str(y_utm[j]))
                        file.write(str(","))
                        file.write(str(botm[i,j]))
                        file.write(str(","))
                        file.write(str(head[k,i,0,j]))
                        file.write(str(","))
                        file.write(str(conc[k+1,i,0,j]) + "\n")   #\n for next line!
                        
                file.write("TEXT X=6 Y=83 T=\"Time=")
                file.write(str((k+1)*perlen_ys))
                file.write("years\", F=TIMES-BOLD, CS=FRAME, H=3, ZN= ")
                file.write(str(k+1)+"\n")
            file.close()
            
            #make csv file
            for k in range(0, nstp[0]):
                for i in range(0,nlay):
                    with open(os.path.join(output_dir, "3_SP_" +str(k+1) + '_layer_' + str(i+1)+ ".csv"), mode='a+') as all_points:
                        all_points_writer = csv.writer(all_points, delimiter=' ')
                        for j in range(0, (len(pt_x_lst))):
                                all_points_writer = csv.writer(all_points, delimiter=' ')
                                all_points_writer.writerow(
                                        [str(x_utm[j]),
                                         str(y_utm[j]),
                                         str(botm[i,j]),
                                         str(head[k,i,0,j]),
                                         str(conc[k+1,i,0,j])])
           

