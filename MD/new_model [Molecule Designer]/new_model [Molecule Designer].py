#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#          Molecule Designer - User panel          |
#                                                  |
#                           © 2021 Juan José Arango|
#                  Universidad Nacional de Colombia|
#                                                  |
#    Licensed under the Apache License, Version 2.0|
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Spectrum and field distribution simulation for photonic circuits#


import sys
sys.path.insert(1, './../MoleculeDesigner')
sys.path.insert(2, './../MoleculeDesigner/software')
from molecule_designer import run_engine
from obj_arrow import emerge_cmap as emerge
from obj_arrow import birds_cmap as birds


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# (( Model Title ))
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 


# Data input settings 📂
directory_dev_config=     './Device_config'          #Path of the directory in which modelled circuit's graph information is stored.
adj_matrix_file=          '/adj_M.txt'               #TXT file within directory_dev_config. [*] It is required that adjacency-matrix's first and last rows/columns correspond to the reference and signal evaluation nodes. Matrix file format: rows separated in different text lines, column elements separated by commas.

coord_matrix_from_image=  True                       #If True: only 1 circuit's graph image should be stored in directory_dev_config. PNG and JPG formats supported. If False: provide an external TXT file with graph nodes coordinates in external_coord_matrix.
external_coord_matrix=    '/external_coords.txt'     #(optional) TXT file within directory_dev_config. [*] It is required that external_coord_matrix's rows preserve the row/column order set by adjacency matrix (both refer to a common graph nodes numbering). External coordinates format: coordinates ordered pairs separated in different text lines, each coordinate comprises two numbers separated by a comma, following the conventional pixel coordinate system.


#Study sweep settings 📐                             #Wavelength/Frequency mode selection is made during code execution, just include the sweep information of the mode to be selected.
lambda0_start=             0  #um (micrometers)        #Wavelength domain study
lambda0_stop=              0  #um (micrometers)
lambda0_samples=           1000

nu_start=                  0 #THz (e12 Hz)             #Frequency domain study
nu_stop=                   0 #THz (e12 Hz)
nu_samples=                1000


# Evaluation points for visualization 📌
val_vis=                  []  #um || THz             #Include evaluation points within brackets, separated by commas. If blank, no visualization is performed.


# Results settings 📋
directory_results=        './Results'     #If directory does not exist, a new one is created. [*] If the same model is executed several times without updating directory_results, results files might get overwritten.
open_interactive_window=  True            #If True, plots and graphics are presented in interactive visualization interfaces after completing execution.


#     ∙Transmission plot 📊
transmission_plot_ylims=  [[],[]]         #Transmittance spectrum y-axis limits: [[<y_min>],[<y_max>]]. [*]Leave it blank for auto-set. [**]It is recommended to plot first without restrictions, in order to check data range.
plot_log_transmission=    False           #If True, transmittance spectrum is plotted in dB.
generate_T_file=          True            #If True, transmittance plot data are stored in a text file.

#     ∙Phase plot 📈
generate_P_file=          True            #If True, phase-shift plot data are stored in a text file.

#     ∙Graphics 👁
generate_energy_dist=     False           #If True, mode's complex amplitude norm graphic is generated.
colormap_style_energy=    birds           #Energy graphics cmap [*]default: 'YlGnBu' [**]Compatible with python's cmaps strings and (customizable) cmap objects (additional cmaps can be included in and imported from MoleculeDesigner\software\obj_arrow).
saturation_level_energy=  1               #Field norm distribution saturation level; [*] necessarily within (0,1).

generate_phase_dist=      False           #If True, mode's phase-shift graphic is generated. [*] In order to get a correct phase distribution, it is necessary to provide a sufficiently narrow discretization for propagation sections (at graph definition stage). Required to avoid sampling error of the periodic phase distribution.
colormap_style_phase=     emerge          #Phase graphics cmap [*]default: 'twilight' [**]Compatible with python's cmaps strings and (customizable) cmap objects (additional cmaps can be included in and imported from MoleculeDesigner\software\obj_arrow).

waveguide_graphic_width=  5               #Waveguide schematic line width.
include_color_bar=        True            #If True, a color bar legend is included.
include_direction_arrows= True            #If True, signal-flow direction is indicated with arrows.
include_local_eval_tags=  False           #If True, field norm or phase-shift local values are explicitly plotted as text labels.
eval_tags_fontsize=       5               #Local evaluation tags fontsize.
image_dpi=                500             #Field distribution graphic detail degree.

save_figs=                True            #If True, all generated plots and graphics are saved as PNG or PDF files, in directory_results.

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

#  🚀
run_engine(lambda0_start, lambda0_stop, lambda0_samples, nu_start, nu_stop, nu_samples, val_vis, directory_dev_config, directory_results, adj_matrix_file, coord_matrix_from_image, external_coord_matrix, generate_energy_dist, generate_phase_dist, colormap_style_energy, colormap_style_phase, include_color_bar, include_direction_arrows, saturation_level_energy, waveguide_graphic_width, save_figs, image_dpi, include_local_eval_tags, eval_tags_fontsize, open_interactive_window, transmission_plot_ylims, plot_log_transmission, generate_T_file, generate_P_file)