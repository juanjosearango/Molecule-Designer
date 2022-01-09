#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#            Molecule Designer - Engine            |
#                                                  |
#                           © 2021 Juan José Arango|
#                  Universidad Nacional de Colombia|
#                                                  |
#    Licensed under the Apache License, Version 2.0|
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Libraries and resources to be implemented

import os
original_wd=os.getcwd()
os.chdir('./../MoleculeDesigner')
from USER_CONFIGURED_database_coupling_coefficients import kappa
from USER_CONFIGURED_database_propagation_parameters import p
os.chdir(original_wd)
from pixel_selector import Selector
from obj_arrow import color_arrow
from obj_arrow import great_coloring
from obj_triangle import triangle

#(External)
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import pandas as pd
import numpy as np
import warnings
import datetime
import openpyxl
import time
import copy
import math

I=0+1j

def run_engine(LAMBDA0_START, LAMBDA0_STOP, LAMBDA0_SAMPLES, NU_START, NU_STOP, NU_SAMPLES, VAL_VIS, DIR_DEV_CONFIG, DIR_RESULTS, ADJ, COORD_FROM_IMAGE, COORD, GEN_E_DIST, GEN_P_DIST, CMAP_STYLE_E, CMAP_STYLE_P, INCLUDE_CBAR, INCLUDE_DIR_ARROWS, SAT_LEV_E, WG_GRAPHIC_WIDTH, SAVE_FIGS, IMAGE_DPI, INCLUDE_LOCAL_TAGS, EVAL_TAGS_FONTSIZE, OPEN_INTERACTIVE_W, T_PLOT_YLIMS, PLOT_LOG_T, GENERATE_T_FILE, GENERATE_P_FILE):
    
    chronos=time.perf_counter()
    
    #Study settings unloading
    lambda0_start=LAMBDA0_START
    lambda0_stop=LAMBDA0_STOP
    lambda0_samples=LAMBDA0_SAMPLES
    nu_start=NU_START
    nu_stop=NU_STOP
    nu_samples=NU_SAMPLES
    val_vis=VAL_VIS
    carpeta=DIR_RESULTS
    generate_energy_dist=GEN_E_DIST
    generate_phase_dist=GEN_P_DIST
    colormap_style_energy_pre=CMAP_STYLE_E
    colormap_style_phase_pre=CMAP_STYLE_P
    include_color_bar=INCLUDE_CBAR
    include_arrows=INCLUDE_DIR_ARROWS
    saturation_level_energy=SAT_LEV_E
    waveguide_graphic_width=WG_GRAPHIC_WIDTH
    save_figs=SAVE_FIGS
    image_dpi=IMAGE_DPI
    
    #Molecule Designer header printing
    print('                 ╔═════════════════════════════╗')
    print('                 ║      MOLECULE DESIGNER      ║')
    print('                 ║  Photonic device simulator  ║')
    print('                 ╚═════════════════════════════╝')
    print('')
    print('╭────────────────────────╮              ╭────────────────────────╮')
    print('│      ╭──╮╭──╮╭──╮      │╭────╮  ╭────╮│ ╭────────╮ ╭────────╮  │')
    print('│      ╰╮╭╯╰╮╭╯╰╮╭╯      ││    ╰──╯    ││ ╰────────╯ ╰────────╯  │')
    print('│  ╭────╯╰──╯╰──╯╰────╮  ││    ╭──╮    ││╭─╮╭──╮╭─╮╭──╮╭─╮╭──╮╭─╮│')
    print('│  ╰──┉┉┉┉┉───┉┉┉┉┉───╯  │╰────╯  ╰────╯│╰─╯╰──╯╰─╯╰──╯╰─╯╰──╯╰─╯│')
    print('│╭─╮╭──╮╭─╮╭──╮╭─╮╭──╮╭─╮│╭────╮  ╭────╮│  ╭──┉┉┉┉┉───┉┉┉┉┉───╮  │')
    print('│╰─╯╰──╯╰─╯╰──╯╰─╯╰──╯╰─╯││    ╰──╯    ││  ╰────╮╭──╮╭──╮╭────╯  │')
    print('│ ╭────────╮ ╭────────╮  ││    ╭──╮    ││      ╭╯╰╮╭╯╰╮╭╯╰╮      │')
    print('│ ╰────────╯ ╰────────╯  │╰────╯  ╰────╯│      ╰──╯╰──╯╰──╯      │')
    print('╰────────────────────────╯              ╰────────────────────────╯ ')
    print('')
    print('                                    ©Juan José Arango-Uribe (2021)')
    print('                                  Universidad Nacional de Colombia')
    print('')
    print('')
    
    #Adjacency matrix information import
    adjacency_df= pd.read_csv(DIR_DEV_CONFIG+ADJ, header=None)
    
    #Circuit nodes coordinates information generation/import
    if(generate_energy_dist or generate_phase_dist):
        if(COORD_FROM_IMAGE):
            print('Use the image pop-up window for providing graph nodes coordinates by clicking on nodes locations.\nIMPORTANT: The node numbering convention used for adjacency matrix construction needs to be preserved in coordinates selection.\n\n To restart the selection press r (lower case).\n To finish the selection press n (lower case). You can close the image window after code execution.\n Pressing Esc will stop the program.\n')
            pix_selector = Selector(DIR_DEV_CONFIG)
            pix_selector.run()
            wb_obj = openpyxl.load_workbook(DIR_DEV_CONFIG+str("/coords.xlsx"))
            sheet = wb_obj.active
            string_coords=sheet["C2"].value.replace('], [','];[').replace(', ',',').replace('[','').replace(']','')
            XY_list_coords=string_coords.split(';')
            for i in np.arange(len(XY_list_coords)):
                XY_list_coords[i]=XY_list_coords[i].split(',')
            for i in np.arange(len(XY_list_coords)):
                for j in np.arange(len(XY_list_coords[i])):
                    XY_list_coords[i][j]=int(XY_list_coords[i][j])
            XY=np.array(XY_list_coords)
        else:
            coord_df= pd.read_csv(DIR_DEV_CONFIG+COORD, header=None)
            XY=coord_df.to_numpy()
            
        #Coordinates matrix adjustment
        X_min=min(XY[:,0])
        X_max=max(XY[:,0])
        Y_min=min(XY[:,1])
        Y_max=max(XY[:,1])
        
        height=Y_max-Y_min
        width=X_max-X_min
        
        if(height>width):
            molec_vertical=True
        else:
            molec_vertical=False
        
        if(molec_vertical):
            X_ref=50-50*width/height
            Y_ref=0    
        else:
            X_ref=0
            Y_ref=50-50*height/(width)
        
        if(include_color_bar):
            factor=90
        else:
            factor=100
        
        for i in np.arange(XY.shape[0]):
            XY[i,0]=X_ref+factor*(XY[i,0]-X_min)/max(width,height)
            XY[i,1]=Y_ref+100*(height-XY[i,1]+Y_min)/max(width,height)
        
        #Colormap legend preparation
        if(include_color_bar):
            color_level_coords=[[max(XY[:,0])+10,min(XY[:,1])+5],[max(XY[:,0])+10,max(XY[:,1])-5]]
        else:
            color_level_coords=['X']
    
    #x-axis setting
    displayed_in_freq=int(input('Please select the spectral domain:'+'\n'+'\n'+'[0] Wavelength'+'\n'+'[1] Frequency'+'\n\n'+'>>: '))
    
    print()
    print('🚩 🚩 Running simulation 🚩 🚩')
    print()
    
    #Results directory creation (if it does not exist already)
    if (not os.path.exists(carpeta)):
        os.makedirs(carpeta)
        print('✓  Results directory created.')
        print()
    
    #Visualization data structures preparation
    if(displayed_in_freq==1):
        for i in np.arange(len(val_vis)):
            val_vis[i]=299.792458/val_vis[i]

    VIS=[]
    for LLL in val_vis:
        VIS.append([LLL])

    A=adjacency_df.to_numpy()
    A=A[:,:A.shape[1]-1]
    A=np.array(A, dtype='complex')
    
    
    #Visualization colormaps preparation
    if(type(colormap_style_energy_pre)!=type(cm.get_cmap('Blues'))):
        colormap_style_energy=cm.get_cmap(colormap_style_energy_pre)
    else:
        colormap_style_energy=colormap_style_energy_pre
        
    if(type(colormap_style_phase_pre)!=type(cm.get_cmap('Blues'))):
        colormap_style_phase=cm.get_cmap(colormap_style_phase_pre)
    else:
        colormap_style_phase=colormap_style_phase_pre
    
    
    
    #Inclusion/Hiding of signal direction arrows
    if(include_arrows):
        opacity=1
    else:
        opacity=0
    
    #Connectivity matrix creation (for visualization)
    A_u=[]
    for elem in A:
        A_u.append(elem)
    A_u=np.array(A_u)
    
    for i in np.arange(A_u.shape[0]):
        for j in np.arange(A_u.shape[1]):
            if(A_u[i][j]>299):
                A_u[i][j]=1
            elif((A_u[i][j]<200 and A_u[i][j]>0)):
                A_u[i][j]=-1
            else:
                A_u[i][j]=0
    
    #Adjacency matrix processing (for spectral analysis)
    warnings.filterwarnings(action='ignore', message='Casting')
    
    used_parameters_report=['\n'+str(datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S"))+'\n','Parameters of modelled device:','\n','(according to input graph parameter encoding)','\n','\n']
    for i in np.arange(A.shape[0]):
        for j in np.arange(A.shape[1]):
            if (A[i][j]>=100 and A[i][j]<200):              #Transmission coefficients (freq. independent **approximation)
                used_parameters_report.append(str(int(abs(A[i][j])))+': '+str(np.sqrt(1-(kappa[int(abs(A[i][j])-100)])**2))+'\n')
                A[i][j]=np.sqrt(1-(kappa[int(abs(A[i][j])-100)])**2)
            if (A[i][j]>=200 and A[i][j]<300):              #Coupling coefficients (freq. independent **approximation)
                qqq=int(math.floor(A[i][j])-200)
                used_parameters_report.append(str(int(abs(A[i][j])))+': '+str(I*kappa[qqq])+'\n')
                A[i][j]=I*kappa[qqq]
    
    print('✓  Device configuration and visualization matrices processed.')
    print()
    
    #Horizontal axis preparation
    if(displayed_in_freq):
        FREQ=np.linspace(nu_start,nu_stop,nu_samples)
        LDA0=299.792458/FREQ
    else:
        LDA0=np.linspace(lambda0_start,lambda0_stop,lambda0_samples)
    
    print('⏳ Calculating transmission curve...')
    print('⏳ [ 0.0 % ]')
    
    #Transfer functions calculation
    TRANSMITTANCE=[]
    PHASE=[]

    cont=1
    contt=1
    for lambda0 in LDA0:
        A_temp=np.array(A)
        for i in np.arange(A_temp.shape[0]):
            for j in np.arange(A_temp.shape[1]):
                if (abs(A_temp[i][j])>=300):          #Propagation coefficients (freq. dependent)
                    if(lambda0==LDA0[0]):
                        used_parameters_report.append(str(int(abs(A_temp[i][j])))+' (alpha): '+str(p[int(abs(A_temp[i][j]))-300][0])+' [1/um]\n')
                        used_parameters_report.append(str(int(abs(A_temp[i][j])))+' (n_eff): '+str(p[int(abs(A_temp[i][j]))-300][1])+'\n')
                        used_parameters_report.append(str(int(abs(A_temp[i][j])))+' (L_eff_prop): '+str(p[int(abs(A_temp[i][j]))-300][2])+' [um]\n')
                        
                    alpha_temp=p[int(abs(A_temp[i][j]))-300][0]
                    neff_temp=p[int(abs(A_temp[i][j]))-300][1]
                    L_temp=p[int(abs(A_temp[i][j]))-300][2]
                    
                    A_temp[i][j]=np.exp((-0.5*alpha_temp*L_temp)+I*(2*np.pi*neff_temp*L_temp/lambda0))

        U=np.eye(A_temp.shape[0])
        T=np.linalg.inv(U-A_temp)

        TT=T[0][A_temp.shape[1]-1]*np.conj(T[0][A_temp.shape[1]-1])
        TP=np.angle(T[0][A_temp.shape[1]-1])
        TRANSMITTANCE.append(TT)
        PHASE.append(TP)
        
        if (contt>=0.1*LDA0.shape[0]):
            print('⏳ [ '+str(100*cont/LDA0.shape[0])[0:5]+' % ]')
            contt=0
        cont+=1
        contt+=1
    print('✓  Transmission curve calculated.')
    print()
    
    #Visualization graphics generation
    if(len(VIS)!=0):
        print('⏳ Preparing graphics...')
        for L in VIS:
            if np.size(L)==1:
                lambda0=L[0]
                A_temp=np.array(A)
                for i in np.arange(A_temp.shape[0]):
                    for j in np.arange(A_temp.shape[1]):
                        if (abs(A_temp[i][j])>=300):          #Propagation coefficients (freq. dependent)
                            
                            alpha_temp=p[int(abs(A_temp[i][j]))-300][0]
                            neff_temp=p[int(abs(A_temp[i][j]))-300][1]
                            L_temp=p[int(abs(A_temp[i][j]))-300][2]
                            
                            A_temp[i][j]=np.exp((-0.5*alpha_temp*L_temp)+I*(2*np.pi*neff_temp*L_temp/lambda0))
                
                U=np.eye(A_temp.shape[0])
                T=np.linalg.inv(U-A_temp)
                L.append(T[:])
                
        if(displayed_in_freq==1):
            for i in np.arange(len(VIS)):
                VIS[i][0]=299.792458/VIS[i][0]
        
        VIS_E=copy.deepcopy(VIS)
        VIS_P=copy.deepcopy(VIS)

        #Device visualization
        print('⏳ Plotting and saving...')
        print('⏳ [ 0.0 % ]')
        cont=1
        tot_figs_per_value=0
        if(generate_energy_dist):
            tot_figs_per_value+=2
        if(generate_phase_dist):
            tot_figs_per_value+=2
            
        for i in np.arange(len(VIS)):
            cont2=1
            if(generate_energy_dist):
                #Field modules at matrices preparation calculation
                temp_list=[]
                for j in np.arange(VIS_E[i][1].shape[0]):
                    for k in np.arange(VIS_E[i][1].shape[1]):
                        VIS_E[i][1][j][k]=abs(VIS_E[i][1][j][k])
                        
                        if j==0:
                            temp_list.append(abs(VIS_E[i][1][j][k]))
                            
                VIS_E[i][1]=np.array(VIS_E[i][1],dtype='float')

                for j in np.arange(VIS_E[i][1].shape[0]):
                    for k in np.arange(VIS_E[i][1].shape[1]):
                        VIS_E[i][1][j][k]=VIS_E[i][1][j][k]/(max(temp_list))
                        
                #Plots: Field intensity distribution and spectrum
                if(displayed_in_freq):
                    fig=plt.figure('Energy distribution @'+str(VIS_E[i][0])+'THz_Spectra')        
                else:
                    fig=plt.figure('Energy distribution @'+str(VIS_E[i][0])+'$\mu$m_Spectra')
                    
                grid=plt.GridSpec(20,13,wspace=0.8)
                
                if(molec_vertical):
                    format_adjustment=10
                else:
                    format_adjustment=9
                
                ax1 = fig.add_subplot(grid[1:-2,0:format_adjustment])
                ax1.set_aspect('equal')
                
                if(displayed_in_freq):        
                    ax1.set_title('Energy distribution @'+str(VIS_E[i][0])+'THz',loc='left',fontdict={'fontsize':8})
                else:        
                    ax1.set_title('Energy distribution @'+str(VIS_E[i][0])+'$\mu$m',loc='left',fontdict={'fontsize':8})
                
                ax1.axis("off")
                
                if(color_level_coords!=['X']):
                    color_arrow(0.999999999,1,ax1, [color_level_coords[0][0],color_level_coords[0][1]-0.3], [color_level_coords[1][0],color_level_coords[1][1]+0.3], cmap=cm.get_cmap('gray'), n=500, lw=waveguide_graphic_width+2.7, sat=saturation_level_energy, include_tags=False)
                    color_arrow(0,1,ax1, color_level_coords[0], color_level_coords[1], cmap=colormap_style_energy, n=500, lw=waveguide_graphic_width+2.5, sat=saturation_level_energy, include_tags=False)
                    ax1.text(color_level_coords[0][0]-0.8,color_level_coords[0][1]+1,'0',rotation=90,fontsize=5,color='white')
                    ax1.text(color_level_coords[0][0]-1,color_level_coords[0][1]+1.2,'0',rotation=90,fontsize=5,color='white')
                    ax1.text(color_level_coords[0][0]-0.9,color_level_coords[0][1]+1.1,'0',rotation=90,fontsize=5)
                    ax1.text(color_level_coords[1][0]-0.8,color_level_coords[1][1]-9,str(float(max(temp_list)))[0:5],rotation=90,fontsize=7,color='white')
                    ax1.text(color_level_coords[1][0]-1,color_level_coords[1][1]-8.8,str(float(max(temp_list)))[0:5],rotation=90,fontsize=7,color='white')
                    ax1.text(color_level_coords[1][0]-0.9,color_level_coords[1][1]-8.9,str(float(max(temp_list)))[0:5],rotation=90,fontsize=7)
    
                
                for iii in np.arange(A_u.shape[0]):
                    for jjj in np.arange(A_u.shape[1]):
                        if(abs(A_u[iii][jjj])==1):
                            if(abs(A_u[jjj][iii])==1):
                                if((XY[iii][0]-XY[jjj][0])!=0):
                                    pendiente=-1/((XY[iii][1]-XY[jjj][1])/(XY[iii][0]-XY[jjj][0]))
                                    dx=2*np.cos(np.arctan(pendiente))
                                    dy=2*np.sin(np.arctan(pendiente))
                                else:
                                    dx=2.5
                                    dy=0
                                if(iii<jjj):
                                    triangle(VIS_E[i][1][0][iii],VIS_E[i][1][0][jjj],ax1, (XY[iii]+[dx,dy]).tolist(), (XY[jjj]+[dx,dy]).tolist(),opacity)
                                else:
                                    triangle(VIS_E[i][1][0][iii],VIS_E[i][1][0][jjj],ax1, (XY[iii]-[dx,dy]).tolist(), (XY[jjj]-[dx,dy]).tolist(),opacity)
                            else:
                                triangle(VIS_E[i][1][0][iii],VIS_E[i][1][0][jjj],ax1, XY[iii].tolist(), XY[jjj].tolist(),opacity)
                                      
                for iii in np.arange(A_u.shape[0]):
                    for jjj in np.arange(A_u.shape[1]):
                        if(A_u[iii][jjj]==1):
                            s_t='p'
                        elif(A_u[iii][jjj]==-1):
                            s_t='t'
                        if(abs(A_u[iii][jjj])==1):
                            if(abs(A_u[jjj][iii])==1):
                                if((XY[iii][0]-XY[jjj][0])!=0):
                                    pendiente=-1/((XY[iii][1]-XY[jjj][1])/(XY[iii][0]-XY[jjj][0]))
                                    dx=2*np.cos(np.arctan(pendiente))
                                    dy=2*np.sin(np.arctan(pendiente))
                                else:
                                    dx=2.5
                                    dy=0
                                if(iii<jjj):
                                    color_arrow(VIS_E[i][1][0][iii],VIS_E[i][1][0][jjj],ax1, (XY[iii]+[dx,dy]).tolist(), (XY[jjj]+[dx,dy]).tolist(), cmap=colormap_style_energy, n=500, lw=waveguide_graphic_width, sat=saturation_level_energy, segment_type=s_t, include_tags=INCLUDE_LOCAL_TAGS, eval_tags_fontsize=EVAL_TAGS_FONTSIZE)
                                else:
                                    color_arrow(VIS_E[i][1][0][iii],VIS_E[i][1][0][jjj],ax1, (XY[iii]-[dx,dy]).tolist(), (XY[jjj]-[dx,dy]).tolist(), cmap=colormap_style_energy, n=500, lw=waveguide_graphic_width, sat=saturation_level_energy, segment_type=s_t, include_tags=INCLUDE_LOCAL_TAGS, eval_tags_fontsize=EVAL_TAGS_FONTSIZE)
                            else:
                                color_arrow(VIS_E[i][1][0][iii],VIS_E[i][1][0][jjj],ax1, XY[iii].tolist(), XY[jjj].tolist(), cmap=colormap_style_energy, n=500, lw=waveguide_graphic_width, sat=saturation_level_energy, segment_type=s_t, include_tags=INCLUDE_LOCAL_TAGS, eval_tags_fontsize=EVAL_TAGS_FONTSIZE)
                            
                plt.tick_params(left=False,
                                bottom=False,
                                labelleft=False,
                                labelbottom=False)
                
                fig.add_subplot(grid[1:7,10:])
                
                if (displayed_in_freq):
                    if(PLOT_LOG_T):
                        plt.plot(FREQ,10*np.log10(TRANSMITTANCE))
                        plt.vlines(VIS_E[i][0],min(10*np.log10(TRANSMITTANCE)),max(10*np.log10(TRANSMITTANCE)),colors='crimson',zorder=5,linestyle='dashed')
                        plt.ylabel('Transmission (dB)',fontdict={'fontsize':7})
                    else:
                        plt.plot(FREQ,TRANSMITTANCE)
                        plt.vlines(VIS_E[i][0],min(TRANSMITTANCE),max(TRANSMITTANCE),colors='crimson',zorder=5,linestyle='dashed')
                        plt.ylabel('Transmission (adim.)',fontdict={'fontsize':7})
                    plt.xticks(fontsize=7)
                    plt.xlabel('Frequency (THz)',fontdict={'fontsize':7})
                    
                else:
                    if(PLOT_LOG_T):
                        plt.plot(LDA0,10*np.log10(TRANSMITTANCE),zorder=1)
                        plt.vlines(VIS_E[i][0],min(10*np.log10(TRANSMITTANCE)),max(10*np.log10(TRANSMITTANCE)),colors='crimson',zorder=5,linestyle='dashed')
                        plt.ylabel('Transmission (dB)',fontdict={'fontsize':7})
                    else:
                        plt.plot(LDA0,TRANSMITTANCE,zorder=1)
                        plt.vlines(VIS_E[i][0],min(TRANSMITTANCE),max(TRANSMITTANCE),colors='crimson',zorder=5,linestyle='dashed')
                        plt.ylabel('Transmission (adim.)',fontdict={'fontsize':7})
                    plt.xticks(fontsize=7)
                    plt.xlabel('Wavelength in vacuum ($\mu$m)',fontdict={'fontsize':7})
                plt.ticklabel_format(useOffset=False)  
                plt.yticks(fontsize=7)
                if(len(T_PLOT_YLIMS[0])!=0 and len(T_PLOT_YLIMS[1])!=0):
                    plt.ylim(bottom=T_PLOT_YLIMS[0][0],top=T_PLOT_YLIMS[1][0])
                
                if(save_figs):
                    plt.savefig(carpeta+str('\EnergyDistribution@')+str(VIS_E[i][0])+str('_Spectra.pdf'),dpi=image_dpi, transparent=True, bbox_inches='tight')
                plt.show()
                print('⏳ [ '+str(100*(tot_figs_per_value*(cont-1)+cont2)/(tot_figs_per_value*len(VIS)))[0:5]+' % ]')
                cont2+=1
                
                #Plots: Just field intensity distribution
                
                if(displayed_in_freq):
                    fig=plt.figure('Energy distribution @'+str(VIS_E[i][0])+'THz')        
                else:
                    fig=plt.figure('Energy distribution @'+str(VIS_E[i][0])+'um')
                
                ax3 = fig.add_subplot(1,1,1)
                ax3.set_aspect('equal')
                
                if(displayed_in_freq):        
                    ax3.set_title('Energy distribution @'+str(VIS_E[i][0])+'THz',loc='left',fontdict={'fontsize':10})
                else:        
                    ax3.set_title('Energy distribution @'+str(VIS_E[i][0])+'um',loc='left',fontdict={'fontsize':10})
                
                ax3.axis("off")
                
                if(color_level_coords!=['X']):
                    color_arrow(0.999999999,1,ax3, [color_level_coords[0][0],color_level_coords[0][1]-0.3], [color_level_coords[1][0],color_level_coords[1][1]+0.3], cmap=cm.get_cmap('gray'), n=500, lw=waveguide_graphic_width+2.7, sat=saturation_level_energy, include_tags=False)
                    color_arrow(0, 1, ax3, color_level_coords[0], color_level_coords[1], cmap=colormap_style_energy, n=500, lw=waveguide_graphic_width+2.5, sat=saturation_level_energy, include_tags=False)
                    ax3.text(color_level_coords[0][0]-0.8,color_level_coords[0][1]+1,'0',rotation=90,fontsize=5,color='white')
                    ax3.text(color_level_coords[0][0]-1,color_level_coords[0][1]+1.2,'0',rotation=90,fontsize=5,color='white')
                    ax3.text(color_level_coords[0][0]-0.9,color_level_coords[0][1]+1.1,'0',rotation=90,fontsize=5)
                    ax3.text(color_level_coords[1][0]-0.8,color_level_coords[1][1]-9,str(float(max(temp_list)))[0:5],rotation=90,fontsize=5,color='white')
                    ax3.text(color_level_coords[1][0]-1,color_level_coords[1][1]-8.8,str(float(max(temp_list)))[0:5],rotation=90,fontsize=5,color='white')
                    ax3.text(color_level_coords[1][0]-0.9,color_level_coords[1][1]-8.9,str(float(max(temp_list)))[0:5],rotation=90,fontsize=5)
                
                for iii in np.arange(A_u.shape[0]):
                    for jjj in np.arange(A_u.shape[1]):
                        if(abs(A_u[iii][jjj])==1):
                            if(abs(A_u[jjj][iii])==1):
                                if((XY[iii][0]-XY[jjj][0])!=0):
                                    pendiente=-1/((XY[iii][1]-XY[jjj][1])/(XY[iii][0]-XY[jjj][0]))
                                    dx=2*np.cos(np.arctan(pendiente))
                                    dy=2*np.sin(np.arctan(pendiente))
                                else:
                                    dx=2.5
                                    dy=0
                                if(iii<jjj):
                                    triangle(VIS_E[i][1][0][iii],VIS_E[i][1][0][jjj],ax3, (XY[iii]+[dx,dy]).tolist(), (XY[jjj]+[dx,dy]).tolist(),opacity)
                                else:
                                    triangle(VIS_E[i][1][0][iii],VIS_E[i][1][0][jjj],ax3, (XY[iii]-[dx,dy]).tolist(), (XY[jjj]-[dx,dy]).tolist(),opacity)
                            else:
                                triangle(VIS_E[i][1][0][iii],VIS_E[i][1][0][jjj],ax3, XY[iii].tolist(), XY[jjj].tolist(),opacity)
                                           
                for iii in np.arange(A_u.shape[0]):
                    for jjj in np.arange(A_u.shape[1]):
                        if(A_u[iii][jjj]==1):
                            s_t='p'
                        elif(A_u[iii][jjj]==-1):
                            s_t='t'
                        if(abs(A_u[iii][jjj])==1):
                            if(abs(A_u[jjj][iii])==1):
                                if((XY[iii][0]-XY[jjj][0])!=0):
                                    pendiente=-1/((XY[iii][1]-XY[jjj][1])/(XY[iii][0]-XY[jjj][0]))
                                    dx=2*np.cos(np.arctan(pendiente))
                                    dy=2*np.sin(np.arctan(pendiente))
                                else:
                                    dx=2.5
                                    dy=0
                                if(iii<jjj):
                                    color_arrow(VIS_E[i][1][0][iii],VIS_E[i][1][0][jjj],ax3, (XY[iii]+[dx,dy]).tolist(), (XY[jjj]+[dx,dy]).tolist(), cmap=colormap_style_energy, n=500,lw=waveguide_graphic_width, sat=saturation_level_energy, segment_type=s_t, include_tags=INCLUDE_LOCAL_TAGS, eval_tags_fontsize=EVAL_TAGS_FONTSIZE)
                                else:
                                    color_arrow(VIS_E[i][1][0][iii],VIS_E[i][1][0][jjj],ax3, (XY[iii]-[dx,dy]).tolist(), (XY[jjj]-[dx,dy]).tolist(), cmap=colormap_style_energy, n=500,lw=waveguide_graphic_width, sat=saturation_level_energy, segment_type=s_t, include_tags=INCLUDE_LOCAL_TAGS, eval_tags_fontsize=EVAL_TAGS_FONTSIZE)
                            else:
                                color_arrow(VIS_E[i][1][0][iii],VIS_E[i][1][0][jjj],ax3, XY[iii].tolist(), XY[jjj].tolist(), cmap=colormap_style_energy, n=500,lw=waveguide_graphic_width, sat=saturation_level_energy, segment_type=s_t, include_tags=INCLUDE_LOCAL_TAGS, eval_tags_fontsize=EVAL_TAGS_FONTSIZE)  
                      
                plt.tick_params(left=False,
                                bottom=False,
                                labelleft=False,
                                labelbottom=False)
                
                if(save_figs):
                    plt.savefig(carpeta+str('\EnergyDistribution@')+str(VIS_E[i][0])+str('.pdf'),dpi=image_dpi, transparent=True, bbox_inches='tight')
                plt.show()
                print('⏳ [ '+str(100*(tot_figs_per_value*(cont-1)+cont2)/(tot_figs_per_value*len(VIS)))[0:5]+' % ]')
                cont2+=1
            
            if(generate_phase_dist):
                
                #Field phases at matrices preparation calculation
                temp_list=[]
                for j in np.arange(VIS_P[i][1].shape[0]):
                    for k in np.arange(VIS_P[i][1].shape[1]):
                        VIS_P[i][1][j][k]=np.angle(VIS_P[i][1][j][k])+np.pi
                        
                        if j==0:
                            temp_list.append(VIS_P[i][1][j][k])
                            
                VIS_P[i][1]=np.array(VIS_P[i][1],dtype='float')
                
                for j in np.arange(VIS_P[i][1].shape[0]):
                    for k in np.arange(VIS_P[i][1].shape[1]):
                        VIS_P[i][1][j][k]=VIS_P[i][1][j][k]/(max(temp_list))
                
                #Plots: Phase distribution and spectrum
                if(displayed_in_freq):
                    fig=plt.figure('Phase distribution @'+str(VIS_P[i][0])+'THz_Spectra')        
                else:
                    fig=plt.figure('Phase distribution @'+str(VIS_P[i][0])+'$\mu$m_Spectra')
                    
                grid=plt.GridSpec(20,13,wspace=0.8)
                
                if(molec_vertical):
                    format_adjustment=10
                else:
                    format_adjustment=9
                                
                ax4 = fig.add_subplot(grid[1:-2,0:format_adjustment])
                ax4.set_aspect('equal')
                
                if(displayed_in_freq):        
                    ax4.set_title('Phase distribution @'+str(VIS_P[i][0])+'THz',loc='left',fontdict={'fontsize':10})
                else:        
                    ax4.set_title('Phase distribution @'+str(VIS_P[i][0])+'$\mu$m',loc='left',fontdict={'fontsize':10})
                
                ax4.axis("off")
                
                if(color_level_coords!=['X']):
                    color_arrow(0.999999999,1,ax4, [color_level_coords[0][0],0.4*(color_level_coords[0][1]-0.3)], [color_level_coords[1][0],0.4*(color_level_coords[1][1]+0.3)], cmap=cm.get_cmap('gray'), n=500, lw=waveguide_graphic_width,phase_plotting=True, include_tags=False)
                    color_arrow(0, 1, ax4, [color_level_coords[0][0],0.4*(color_level_coords[0][1])], [color_level_coords[1][0],0.4*(color_level_coords[1][1])], cmap=colormap_style_phase, n=500, lw=waveguide_graphic_width,phase_plotting=True, include_tags=False)
                    ax4.text(color_level_coords[0][0]-0.8,0.4*(color_level_coords[0][1]+1),'0',rotation=90,fontsize=7)
                    ax4.text(color_level_coords[0][0]-1,0.4*(color_level_coords[0][1]+1.2),'0',rotation=90,fontsize=7)
                    ax4.text(color_level_coords[0][0]-0.9,0.4*(color_level_coords[0][1]+1.1),'0',rotation=90,fontsize=7,color='white')
                    ax4.text(color_level_coords[1][0]-0.8,0.4*(color_level_coords[1][1]-9),'2$\pi$',rotation=90,fontsize=7)
                    ax4.text(color_level_coords[1][0]-1,0.4*(color_level_coords[1][1]-8.8),'2$\pi$',rotation=90,fontsize=7)
                    ax4.text(color_level_coords[1][0]-0.9,0.4*(color_level_coords[1][1]-8.9),'2$\pi$',rotation=90,fontsize=7,color='white')
            
                for iii in np.arange(A_u.shape[0]):
                    for jjj in np.arange(A_u.shape[1]):
                        if(abs(A_u[iii][jjj])==1):
                            if(abs(A_u[jjj][iii])==1):
                                if((XY[iii][0]-XY[jjj][0])!=0):
                                    pendiente=-1/((XY[iii][1]-XY[jjj][1])/(XY[iii][0]-XY[jjj][0]))
                                    dx=2*np.cos(np.arctan(pendiente))
                                    dy=2*np.sin(np.arctan(pendiente))
                                else:
                                    dx=2.5
                                    dy=0
                                if(iii<jjj):
                                    triangle(VIS_P[i][1][0][iii],VIS_P[i][1][0][jjj],ax4, (XY[iii]+[dx,dy]).tolist(), (XY[jjj]+[dx,dy]).tolist(),opacity)
                                else:
                                    triangle(VIS_P[i][1][0][iii],VIS_P[i][1][0][jjj],ax4, (XY[iii]-[dx,dy]).tolist(), (XY[jjj]-[dx,dy]).tolist(),opacity)
                            else:
                                triangle(VIS_P[i][1][0][iii],VIS_P[i][1][0][jjj],ax4, XY[iii].tolist(), XY[jjj].tolist(),opacity)
                                      
                for iii in np.arange(A_u.shape[0]):
                    for jjj in np.arange(A_u.shape[1]):
                        if(A_u[iii][jjj]==1):
                            s_t='p'
                        elif(A_u[iii][jjj]==-1):
                            s_t='t'
                        if(abs(A_u[iii][jjj])==1):
                            if(abs(A_u[jjj][iii])==1):
                                if((XY[iii][0]-XY[jjj][0])!=0):
                                    pendiente=-1/((XY[iii][1]-XY[jjj][1])/(XY[iii][0]-XY[jjj][0]))
                                    dx=2*np.cos(np.arctan(pendiente))
                                    dy=2*np.sin(np.arctan(pendiente))
                                else:
                                    dx=2.5
                                    dy=0
                                if(iii<jjj):
                                    color_arrow(VIS_P[i][1][0][iii],VIS_P[i][1][0][jjj],ax4, (XY[iii]+[dx,dy]).tolist(), (XY[jjj]+[dx,dy]).tolist(), cmap=colormap_style_phase, n=500,lw=waveguide_graphic_width,phase_plotting=True, segment_type=s_t, include_tags=INCLUDE_LOCAL_TAGS, eval_tags_fontsize=EVAL_TAGS_FONTSIZE)
                                else:
                                    color_arrow(VIS_P[i][1][0][iii],VIS_P[i][1][0][jjj],ax4, (XY[iii]-[dx,dy]).tolist(), (XY[jjj]-[dx,dy]).tolist(), cmap=colormap_style_phase, n=500,lw=waveguide_graphic_width,phase_plotting=True, segment_type=s_t, include_tags=INCLUDE_LOCAL_TAGS, eval_tags_fontsize=EVAL_TAGS_FONTSIZE)
                            else:
                                color_arrow(VIS_P[i][1][0][iii],VIS_P[i][1][0][jjj],ax4, XY[iii].tolist(), XY[jjj].tolist(), cmap=colormap_style_phase, n=500,lw=waveguide_graphic_width,phase_plotting=True, segment_type=s_t, include_tags=INCLUDE_LOCAL_TAGS, eval_tags_fontsize=EVAL_TAGS_FONTSIZE)
                              
                plt.tick_params(left=False,
                                bottom=False,
                                labelleft=False,
                                labelbottom=False)
                
                fig.add_subplot(grid[1:7,10:])
                
                if (displayed_in_freq):
                    if(PLOT_LOG_T):
                        plt.plot(FREQ,10*np.log10(TRANSMITTANCE))
                        plt.vlines(VIS_P[i][0],min(10*np.log10(TRANSMITTANCE)),max(10*np.log10(TRANSMITTANCE)),colors='crimson',zorder=5,linestyle='dashed')
                        plt.ylabel('Transmission (dB)',fontdict={'fontsize':7})
                    else:
                        plt.plot(FREQ,TRANSMITTANCE)
                        plt.vlines(VIS_P[i][0],min(TRANSMITTANCE),max(TRANSMITTANCE),colors='crimson',zorder=5,linestyle='dashed')
                        plt.ylabel('Transmission (adim.)',fontdict={'fontsize':7})
                    plt.xticks(fontsize=7)
                    plt.xlabel('Frequency (THz)',fontdict={'fontsize':7})
                    
                else:
                    if(PLOT_LOG_T):
                        plt.plot(LDA0,10*np.log10(TRANSMITTANCE),zorder=1)
                        plt.vlines(VIS_P[i][0],min(10*np.log10(TRANSMITTANCE)),max(10*np.log10(TRANSMITTANCE)),colors='crimson',zorder=5,linestyle='dashed')
                        plt.ylabel('Transmission (dB)',fontdict={'fontsize':7})
                    else: 
                        plt.plot(LDA0,TRANSMITTANCE,zorder=1)
                        plt.vlines(VIS_P[i][0],min(TRANSMITTANCE),max(TRANSMITTANCE),colors='crimson',zorder=5,linestyle='dashed')
                        plt.ylabel('Transmission (adim.)',fontdict={'fontsize':7})
                    plt.xticks(fontsize=7)
                    plt.xlabel('Wavelength in vacuum ($\mu$m)',fontdict={'fontsize':7})
                    
                plt.ticklabel_format(useOffset=False)  
                plt.yticks(fontsize=7)
                if(len(T_PLOT_YLIMS[0])!=0 and len(T_PLOT_YLIMS[1])!=0):
                    plt.ylim(bottom=T_PLOT_YLIMS[0][0],top=T_PLOT_YLIMS[1][0])
                
                if(save_figs):
                    plt.savefig(carpeta+str('\PhaseDistribution@')+str(VIS_P[i][0])+str('_Spectra.pdf'),dpi=image_dpi, transparent=True, bbox_inches='tight')
                plt.show()
                print('⏳ [ '+str(100*(tot_figs_per_value*(cont-1)+cont2)/(tot_figs_per_value*len(VIS)))[0:5]+' % ]')
                cont2+=1
                
                #Plots: Just phase distribution
                
                if(displayed_in_freq):
                    fig=plt.figure('Phase distribution @'+str(VIS_P[i][0])+'THz')        
                else:
                    fig=plt.figure('Phase distribution @'+str(VIS_P[i][0])+'um')
                
                ax6 = fig.add_subplot(1,1,1)
                ax6.set_aspect('equal')
                
                if(displayed_in_freq):        
                    ax6.set_title('Phase distribution @'+str(VIS_P[i][0])+'THz',loc='left',fontdict={'fontsize':10})
                else:        
                    ax6.set_title('Phase distribution @'+str(VIS_P[i][0])+'um',loc='left',fontdict={'fontsize':10})
                
                ax6.axis("off")
                
                if(color_level_coords!=['X']):
                    color_arrow(0.999999999,1,ax6, [color_level_coords[0][0],0.4*(color_level_coords[0][1]-0.3)], [color_level_coords[1][0],0.4*(color_level_coords[1][1]+0.3)], cmap=cm.get_cmap('gray'), n=500, lw=waveguide_graphic_width,phase_plotting=True, include_tags=False)
                    color_arrow(0, 1, ax6, [color_level_coords[0][0],0.4*(color_level_coords[0][1])], [color_level_coords[1][0],0.4*(color_level_coords[1][1])], cmap=colormap_style_phase, n=500, lw=waveguide_graphic_width,phase_plotting=True, include_tags=False)
                    ax6.text(color_level_coords[0][0]-0.8,0.4*(color_level_coords[0][1]+1),'0',rotation=90,fontsize=7)
                    ax6.text(color_level_coords[0][0]-1,0.4*(color_level_coords[0][1]+1.2),'0',rotation=90,fontsize=7)
                    ax6.text(color_level_coords[0][0]-0.9,0.4*(color_level_coords[0][1]+1.1),'0',rotation=90,fontsize=7,color='white')
                    ax6.text(color_level_coords[1][0]-0.8,0.4*(color_level_coords[1][1]-9),'2$\pi$',rotation=90,fontsize=7)
                    ax6.text(color_level_coords[1][0]-1,0.4*(color_level_coords[1][1]-8.8),'2$\pi$',rotation=90,fontsize=7)
                    ax6.text(color_level_coords[1][0]-0.9,0.4*(color_level_coords[1][1]-8.9),'2$\pi$',rotation=90,fontsize=7,color='white')
            
                for iii in np.arange(A_u.shape[0]):
                    for jjj in np.arange(A_u.shape[1]):
                        if(abs(A_u[iii][jjj])==1):
                            if(abs(A_u[jjj][iii])==1):
                                if((XY[iii][0]-XY[jjj][0])!=0):
                                    pendiente=-1/((XY[iii][1]-XY[jjj][1])/(XY[iii][0]-XY[jjj][0]))
                                    dx=2*np.cos(np.arctan(pendiente))
                                    dy=2*np.sin(np.arctan(pendiente))
                                else:
                                    dx=2.5
                                    dy=0
                                if(iii<jjj):
                                    triangle(VIS_P[i][1][0][iii],VIS_P[i][1][0][jjj],ax6, (XY[iii]+[dx,dy]).tolist(), (XY[jjj]+[dx,dy]).tolist(),opacity)
                                else:
                                    triangle(VIS_P[i][1][0][iii],VIS_P[i][1][0][jjj],ax6, (XY[iii]-[dx,dy]).tolist(), (XY[jjj]-[dx,dy]).tolist(),opacity)
                            else:
                                triangle(VIS_P[i][1][0][iii],VIS_P[i][1][0][jjj],ax6, XY[iii].tolist(), XY[jjj].tolist(),opacity)
                                          
                for iii in np.arange(A_u.shape[0]):
                    for jjj in np.arange(A_u.shape[1]):
                        if(A_u[iii][jjj]==1):
                            s_t='p'
                        elif(A_u[iii][jjj]==-1):
                            s_t='t'
                        if(abs(A_u[iii][jjj])==1):
                            if(abs(A_u[jjj][iii])==1):
                                if((XY[iii][0]-XY[jjj][0])!=0):
                                    pendiente=-1/((XY[iii][1]-XY[jjj][1])/(XY[iii][0]-XY[jjj][0]))
                                    dx=2*np.cos(np.arctan(pendiente))
                                    dy=2*np.sin(np.arctan(pendiente))
                                else:
                                    dx=2.5
                                    dy=0
                                if(iii<jjj):
                                    color_arrow(VIS_P[i][1][0][iii],VIS_P[i][1][0][jjj],ax6, (XY[iii]+[dx,dy]).tolist(), (XY[jjj]+[dx,dy]).tolist(), cmap=colormap_style_phase, n=1000,lw=waveguide_graphic_width,phase_plotting=True, segment_type=s_t, include_tags=INCLUDE_LOCAL_TAGS, eval_tags_fontsize=EVAL_TAGS_FONTSIZE)
                                else:
                                    color_arrow(VIS_P[i][1][0][iii],VIS_P[i][1][0][jjj],ax6, (XY[iii]-[dx,dy]).tolist(), (XY[jjj]-[dx,dy]).tolist(), cmap=colormap_style_phase, n=1000,lw=waveguide_graphic_width,phase_plotting=True, segment_type=s_t, include_tags=INCLUDE_LOCAL_TAGS, eval_tags_fontsize=EVAL_TAGS_FONTSIZE)
                            else:
                                color_arrow(VIS_P[i][1][0][iii],VIS_P[i][1][0][jjj],ax6, XY[iii].tolist(), XY[jjj].tolist(), cmap=colormap_style_phase, n=1000,lw=waveguide_graphic_width,phase_plotting=True, segment_type=s_t, include_tags=INCLUDE_LOCAL_TAGS, eval_tags_fontsize=EVAL_TAGS_FONTSIZE)  
                      
                plt.tick_params(left=False,
                                bottom=False,
                                labelleft=False,
                                labelbottom=False)
                
                if(save_figs):
                    plt.savefig(carpeta+str('\PhaseDistribution@')+str(VIS_P[i][0])+str('.pdf'),dpi=image_dpi, transparent=True, bbox_inches='tight')
                plt.show()
                print('⏳ [ '+str(100*(tot_figs_per_value*(cont-1)+cont2)/(tot_figs_per_value*len(VIS)))[0:5]+' % ]')
                cont2+=1
            
            cont+=1
            
            
        print('✓  Visualization graphics completed.')
        print()
    
          
    #Transmittance plot output file generation
    if(GENERATE_T_FILE):
        T_file=open(carpeta+'/Transmission.txt','a')
        T_file.write(str(datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S"))+'\n\n')
        if (displayed_in_freq):
            T_file.write('Frequency (THz), Transmission (a.u.)\n')
            for i in np.arange(len(TRANSMITTANCE)):
                T_file.write(str(FREQ[i])+', '+str(np.real(TRANSMITTANCE[i]))+'\n')
        else:
            T_file.write('Lambda0 (um), Transmission (a.u.)\n')
            for i in np.arange(len(TRANSMITTANCE)):
                T_file.write(str(LDA0[i])+', '+str(np.real(TRANSMITTANCE[i]))+'\n')
        T_file.write('\n\n')
        T_file.close()
        print('✓  Transmission data saved in text file.')
        print()
    
    #Phase plot output file generation
    if(GENERATE_P_FILE):
        P_file=open(carpeta+'/Phase.txt','a')
        P_file.write(str(datetime.datetime.now().strftime("%d-%m-%Y %H:%M:%S"))+'\n\n')
        if (displayed_in_freq):
            P_file.write('Frequency (THz), Phase (rad)\n')
            for i in np.arange(len(PHASE)):
                P_file.write(str(FREQ[i])+', '+str(np.real(PHASE[i]))+'\n')
        else:
            P_file.write('Lambda0 (um), Phase (rad)\n')
            for i in np.arange(len(PHASE)):
                P_file.write(str(LDA0[i])+', '+str(np.real(PHASE[i]))+'\n')
        P_file.write('\n\n')
        P_file.close()
        print('✓  Phase-shift data saved in text file.')
        print()
    
    #Transmission plot
    plt.figure("Evaluated transmission plot")
    AX=plt.gca()
    AX.set_prop_cycle(great_coloring)
    if(PLOT_LOG_T):
        TRANSMITTANCE=10*np.log10(TRANSMITTANCE)
        plt.ylabel('Transmission (dB)')
    else:
        plt.ylabel('Transmission (a.u.)')
    
    if (displayed_in_freq):
        plt.plot(FREQ,TRANSMITTANCE)
        plt.xlabel('Frequency (THz)')
    else:
        plt.plot(LDA0,TRANSMITTANCE)
        plt.xlabel('Wavelength in vacuum ($\mu$m)')
    
    plt.ticklabel_format(useOffset=False)
    if(len(T_PLOT_YLIMS[0])!=0 and len(T_PLOT_YLIMS[1])!=0):
        plt.ylim(bottom=T_PLOT_YLIMS[0][0],top=T_PLOT_YLIMS[1][0])
    
    if(save_figs):
        plt.savefig(carpeta+str('\Transmission')+str('.png'),dpi=image_dpi, transparent=True, bbox_inches='tight')
    plt.show()
    
    #Phase plot
    plt.figure("Evaluated phase plot")
    AX=plt.gca()
    AX.set_prop_cycle(great_coloring)
    plt.ylabel('Phase state (rad)')
    if (displayed_in_freq):
        plt.plot(FREQ,PHASE)
        plt.xlabel('Frequency (THz)')
    else:
        plt.plot(LDA0,PHASE)
        plt.xlabel('Wavelength in vacuum ($\mu$m)')
    plt.ticklabel_format(useOffset=False)
    if(save_figs):
        plt.savefig(carpeta+str('\Phase')+str('.png'),dpi=image_dpi, transparent=True, bbox_inches='tight')
    plt.show()
    
    #Device parameters report generation
    report_file=open(carpeta+'/device_parameters.txt','a')
    report_file.writelines(used_parameters_report)
    report_file.close()
    
    chronos=time.perf_counter()-chronos
    
    #Closing information printing
    print('╭──────────────────────╮')
    print('│                      │')
    print('╰──────────────────────╯')
    print('╭────────────────────────────────────────────╮')
    print('│ ✅ ⁙Simulation successfully completed⁙ 💯 │')
    print('╰────────────────────────────────────────────╯')
    print('                    ╭────────────────────────╮')
    print('                    │       in '+str(chronos)[0:6]+' s      │')
    print('                    ╰────────────────────────╯')
    
    if(not OPEN_INTERACTIVE_W):
        plt.close('all')