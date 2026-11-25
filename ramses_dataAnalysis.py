import numpy as np
import matplotlib.pyplot as plt
import yt
# from mpi4py import MPI 
import argparse
import math
import matplotlib as mpl
from dataclasses import dataclass
from typing import Any
from scipy.fftpack import fft, fftfreq
from friendsOfFriends.matplotlibrc import changePlotSettings
from raghav_core import *

changePlotSettings();

plt.rcParams['text.usetex'] = False

'''
The name of the new particle type will be 'stars'
'''
@yt.particle_filter(requires=["particle_mass"], filtered_type="io")
def stars_user(pfilter, data):
    filter= data[pfilter.filtered_type, "particle_mass"]<0.01;
    return filter;

@yt.particle_filter(requires=["particle_mass"], filtered_type="io")
def darkMatter_user(pfilter, data):
    filter= data[pfilter.filtered_type, "particle_mass"]>0.01;
    return filter;

def open_disc_file(dirs, fileNumber):
    ds = yt.load('{}info_{}.txt'.format(dirs, fileNumber));
    print('reading ramses file:', '{}info_{}.txt'.format(dirs, fileNumber));
    return ds, fileNumber;

def read_ramses_file(file_read_directory, fileNumber):
    #Prepare the file number with leading zeroes
    zeroes = (4 - math.floor(np.log10(fileNumber)))*'0';
    fileNumber = '{}{}'.format(zeroes, fileNumber)   

    file_read_directory = '{}output_{}/'.format(file_read_directory_prefix, fileNumber); #yt reads the directory and gathers the input form individual files inside it

    ds, fileNumber = open_disc_file(file_read_directory, fileNumber);
    return ds, fileNumber;

def plot_particle_mass_projection(nmin, nmax, saveName, file_read_directory_prefix,
                                 plot_dir, axis='z', width=(100, 'kpc')):
    """Save a yt projection of all particles' mass for each RAMSES output.

    The projected field is particle mass, displayed in solar masses. ``width``
    is passed to yt and defaults to a 100 kpc-wide field of view.
    """
    field = ('all', 'particle_mass')
    for output_number in range(nmin, nmax + 1):
        file_number = '{:04d}'.format(output_number)
        output_directory = '{}output_{}/'.format(file_read_directory_prefix, file_number)
        ds, _ = open_disc_file(output_directory, file_number)

        projection = yt.ParticleProjectionPlot(ds, axis, [field], width=width)
        projection.set_unit(field, 'Msun')
        output_name = '{}{}_{}'.format(plot_dir, saveName, file_number)
        print('Saving particle mass projection:', output_name)
        projection.save(output_name)

    return

# ====================================
# Some constants
# =====================================

avogadro_no = 6.0221e23;
mean_mass = 1.4;
one_kpc = 3.085677581E+21 #cms
one_pc = one_kpc*0.001
featherLengthScale = one_kpc;
solar_mass = 1.989e33 #grams
one_myr = 3.15576e+13; #seconds

# ====================================
# Runtime parameters
# =====================================

'''
The file numbers that you wish to plot. 

# file_read_directory_prefix = '/home/ra6154ar/ramses_oscar/simulations/isolatedGalaxy/'
# plot_dir = '/home/ra6154ar/ramses_oscar/simulations/isolatedGalaxy/'
# file_read_directory_prefix = '/home/ra6154ar/ramses_updated/simulations/sedov_mhd/';
# plot_dir = '/home/ra6154ar/ramses_updated/simulations/sedov_mhd/';
'''
nmin = 1;
nmax = 11;
# file_read_directory_prefix ='/home/ra6154ar/ramses_nexus/simulations/hbd_10_hydro/'
# file_read_directory_prefix ='/home/ra6154ar/ramses_oscar/simulations/sedov3d_hydro_test/'

fileNumber = 1;
plot_dir='/home/ra6154ar/ramses_oscar/simulations/sedov3d_mhd_test/'
file_read_directory_prefix ='/home/ra6154ar/ramses_oscar/simulations/sedov3d_mhd_test/'
saveName = 'sedov3d_hyrdo'

ds, fileNumber = read_ramses_file(file_read_directory_prefix, fileNumber);

ad = ds.all_data();
density_ppcm3 = ad[('gas', 'density')]*avogadro_no/mean_mass;
radius_pc = ad[('index', 'spherical_radius')]/one_pc;
pres = ad[('gas', 'pressure')]
temp = (pres.v*mean_mass)/(ad[('gas', 'density')].v*const.boltzmann_cgs*const.avogadro_no)
vely = ad[('gas', 'velocity_y')].v/1e5; #km/sec
bz = ad[('gas', 'magnetic_field_z')].v*1e6; #micro-Gauss

#Rho vs R
fig, ax = plt.subplots();
ax.scatter(radius_pc, density_ppcm3, s = 1);
ax.set_xlabel(r'$R \,(\rm pc)$');
ax.set_ylabel(r'$n \,(\rm cm^{-3})$');
print('Saving figure: {}{}_density_R'.format(plot_dir, saveName));
fig.savefig('{}{}_density_R'.format(plot_dir, saveName));

fig, ax = plt.subplots();
ax.scatter(radius_pc, ad[('gas', 'density')], s = 1);
ax.set_xlabel(r'$R \,(\rm pc)$');
ax.set_ylabel(r'$\rho \,(\rm gm\, cm^{-3})$');
print('Saving figure: {}{}_density_cgs_R'.format(plot_dir, saveName));
fig.savefig('{}{}_density_cgs_R'.format(plot_dir, saveName));


#Pressure vs R
fig, ax = plt.subplots();
ax.scatter(radius_pc, temp, s = 1);
ax.set_xlabel(r'$R \,(\rm pc)$');
ax.set_ylabel(r'$T \,(\rm K)$');
print('Saving figure: {}{}_temp_R'.format(plot_dir, saveName));
fig.savefig('{}{}_temp_R'.format(plot_dir, saveName));


#vely vs R
fig, ax = plt.subplots();
ax.scatter(radius_pc, vely, s = 1);
ax.set_xlabel(r'$R \,(\rm pc)$');
ax.set_ylabel(r'$v_y \,(\rm km\, s^{-1})$');
print('Saving figure: {}{}_vely_R'.format(plot_dir, saveName));
fig.savefig('{}{}_vely_R'.format(plot_dir, saveName));


#Bz vs R
fig, ax = plt.subplots();
ax.scatter(radius_pc, bz, s = 1);
ax.set_xlabel(r'$R \,(\rm pc)$');
ax.set_ylabel(r'$B_z \,(\rm \mu G)$');
print('Saving figure: {}{}_Bz_R'.format(plot_dir, saveName));
fig.savefig('{}{}_Bz_R'.format(plot_dir, saveName));

# plot_midPlane_density_profile(fileNumber, width, "density")











# x_axis_label = r'$x$ (\rm code units)'
# y_axis_label = r'$y$ (\rm code units)'
# colobar_label = r'\rm Projected  Density $\left ( code\, units \right )$'
# cmap = 'viridis';
# resolution=512 #taken as [resolution, resolution]
# dpi = 200;
# width = 1.0; #boxlen is 1
# res = [resolution, resolution]




# # ======================================
# # yt projection plot
# # ======================================
# normi =  mpl.colors.LogNorm(vmin = 0.5, vmax = 1.5);
# for i in range(nmin, nmax+1):
#     #Decide the file number
#     fileNumber = i
#     if(fileNumber==0):
#         zeroes = 4*'0'
#     else:
#         zeroes = (4 - math.floor(np.log10(fileNumber)))*'0';
#     fileNumber = '{}{}'.format(zeroes, fileNumber)
#     saveName = 'sedov3d_{}'.format(fileNumber);
#     file_read_directory = '/scratch/ek9/ra8040/ramses/simulations/sedov3d/output_{}/'.format(fileNumber);


#     #Read the file and plot
#     ds, fileNumber = open_disc_file(file_read_directory, fileNumber);


#     proj = ds.proj(('gas', 'density'), 'z');
#     frb = proj.to_frb(width, res);


#     fig, ax = plt.subplots(dpi = dpi);
#     ax.tick_params(color = 'grey', which = 'both');
#     im = ax.imshow(np.array(frb[('gas', 'density')]), cmap = cmap, extent = (-width/2, width/2, -width/2, width/2), norm = normi, origin = 'lower');
#     ax.set_xlabel(x_axis_label);
#     ax.set_ylabel(y_axis_label);
#     ax.set_xlim(-width/2, width/2);

#     clb = fig.colorbar(im);
#     clb.set_label(colobar_label, fontsize = 10);

#     ax.set_ylim(-width/2, width/2);        
#     fig.savefig('{}{}'.format(plot_dir, saveName), dpi = dpi);


#     # ======================================
#     #  Line profiles [scatter plots]
#     # ======================================

