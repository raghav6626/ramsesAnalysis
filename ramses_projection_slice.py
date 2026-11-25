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
from rapack import *
import rapack.constants as const
import cmasher as cmr
from ramses_helper import units_density, units_length, units_time

changePlotSettings();

plt.rcParams['text.usetex'] = False
plt.rcParams.update({'axes.labelsize': 18, 'xtick.labelsize': 14, 'ytick.labelsize': 14,
                     'legend.fontsize': 14});

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
    print('reading ramses file'.format(dirs, fileNumber))
    return ds, fileNumber;

'''
Trying to do standard slice.
Will make profile plots for the density now.
Slicing syntax: 
yt plays you with its conventions...
'''
def plot_midPlane_density_profile(fileNumber, width, field):
    boxLen = 600*const.one_kpc;
    leftEdge = (boxLen-width)/2;
    rightEdge = (boxLen+width)/2;

    zeroes = (4 - math.floor(np.log10(fileNumber)))*'0';
    fileNumber = '{}{}'.format(zeroes, fileNumber);
    saveName = "profilePlot"
    saveName = saveName+'_{}'.format(fileNumber)      

    file_read_directory = '{}output_{}/'.format(file_read_directory_prefix, fileNumber);

    ds, fileNumber = open_disc_file(file_read_directory, fileNumber);  
    # slc = ds.slice('z', coord= 0.5);

    # slc = ds.r[leftEdge:rightEdge, leftEdge:rightEdge, 0]; #this is flattened array
    slc = ds.r[:,:, 0.5]; #note that the 0.5 selects the z = z_mid
    density = slc[("gas", field)];
    x = slc[("gas", "x")];
    y = slc[("gas", "y")];
    z = slc[("gas", "z")];

    fig, ax = plt.subplots();
    ax.hist(x/const.one_kpc);
    ax.set_xlabel(r'$x (kpc)$');
    ax.set_ylabel('PDF');
    fig.savefig('xHist')

    fig, ax = plt.subplots();
    ax.hist(y/const.one_kpc);
    ax.set_xlabel(r'$x (kpc)$');
    ax.set_ylabel('PDF');
    fig.savefig('yHist')

    fig, ax = plt.subplots();
    ax.hist(z/const.one_kpc);
    ax.set_xlabel(r'$z (kpc)$');
    ax.set_ylabel('PDF');
    fig.savefig('zHist')

    frb = yt.FixedResolutionBuffer(slc,(0, 0.5, 0, 0.5) , (512, 512)) #this is the imag.  (xmin, xmax, ymin, ymax) in code_units
    image = frb[("gas", field)];
    print(image);
    x = frb[("gas", 'x')].v.flatten();
    y = frb[("gas", 'y')].v.flatten();

    fig, ax = plt.subplots();
    ax.hist(x/const.one_kpc);
    ax.set_xlabel(r'$x (kpc)$');
    ax.set_ylabel('PDF');
    fig.savefig('xHist_frb')

    fig, ax = plt.subplots();
    ax.hist(y/const.one_kpc);
    ax.set_xlabel(r'$x (kpc)$');
    ax.set_ylabel('PDF');
    fig.savefig('yHist_frb')

    # exit();



    
    trot = 100*const.one_myr;
    cmap = 'viridis';
    width_kpc = width/const.one_kpc;
    normi = mpl.colors.LogNorm(vmin = np.amin(density), vmax = np.amax(density));
    textColor = 'white';
    projectionDirection = 'z';

    time = (ds.current_time.v)/const.one_myr;
    time = time/trot;	
    fontsize = 12;
    fig, ax = plt.subplots();
    im = ax.imshow(image, cmap = cmap, norm = normi, origin = 'lower'); # extent = (-width/2, width/2, -width/2, width/2),
    ax.set_ylabel(r'$y~(\rm kpc)$');
    ax.text(0.07, 0.10, r"$t = %.2f  \, t_{{\rm rot}}$"%time, transform=ax.transAxes, size = fontsize, color = textColor);

    # fig, im, ax = plot_slice(fig, ax,ds,trot, image, cmap, width_kpc, normi, textColor, projectionDirection); 

    clb = fig.colorbar(im);  
    fig.savefig('sliceTest')

    exit();
    density_part = density*const.avogadro_no;
    x = slc[("gas", "x")];
    y = slc[("gas", "y")];
    r_cyl = np.sqrt(x**2 + y**2);

    print(density.units)

    fig, ax = plt.subplots();
    ax.hist(r_cyl/const.one_kpc);
    ax.set_xlabel(r'$\log_{10}\rho (code_units)$');
    ax.set_ylabel('PDF');
    fig.savefig('densityHist')

    return;

def sedovSliceMovie(nmin, nmax, saveName, file_read_directory_prefix, plot_dir, pbeta = False):

    for i in range(nmin, nmax+1):
        print('plotting file number: ', i);
        fileNumber = i;
        zeroes = (4 - math.floor(np.log10(fileNumber)))*'0';
        fileNumber = '{}{}'.format(zeroes, fileNumber)
        saveName = saveName+'_{}'.format(fileNumber)      
    
        file_read_directory = '{}output_{}/'.format(file_read_directory_prefix, fileNumber); #yt reads the directory and gathers the input form individual files inside it

        ds, fileNumber = open_disc_file(file_read_directory, fileNumber);
        slc = yt.SlicePlot(ds, 'z', ("gas", "density"));
        slc.annotate_timestamp(text_args = {'color':'white'});
        slc.annotate_grids();
        slc.save('{}slice_{}'.format(plot_dir, fileNumber));
    
    return;
'''
        ad = ds.all_data();
        x = ad[('x')].v;
        y = ad[('y')].v;
        z = ad[('z')].v;
        density = ad[('density')].
        density_particle_cc = density*avogadro_no; #assuming mean molecular mass = 1.4
        sound_speed = ad[('sound_speed')].v
        magnetic_field_magnitude = ad[('magnetic_field_magnitude')].v

        thermal_pressure = density*(sound_speed)**2;
        magnetic_pressure = magnetic_field_magnitude**2/(8*np.pi);
        ##Re-centering the coordinates to the center of the domain. 
        #Note that the maximum radial coordinate will be sqrt(3 x 0.5^2);
        center = [0.5*const.one_pc, 0.5*const.one_pc, 0.5*const.one_pc];

        x = x - center[0];
        y = y - center[0];
        z = z - center[0];

        radius_spherical = np.sqrt(x**2 + y**2 + z**2);

        fig, ax = plt.subplots();
        ax.scatter(radius_spherical/const.one_pc, density_particle_cc, s = 1);
        ax.set_xlabel(r'$R \,(\rm pc)$');
        ax.set_ylabel(r'$\rho \,(\rm gm \, cm^{-3})$');
        fig.savefig('{}densityProfile'.format(plot_dir));

        fig, ax = plt.subplots();
        ax.hist(sound_speed/1e5);
        ax.set_xlabel(r'$cs (km/sec)$');
        fig.savefig('{}cs_hist'.format(plot_dir));

        # fig, ax = plt.subplots();
        # ax.scatter(radius_spherical/const.one_pc, pressure, s = 1);
        # ax.set_xlabel(r'$R \,(\rm kpc)$');
        # ax.set_ylabel(r'$\log_{10} P \,(\rm code \, units)$');
        # fig.savefig('{}soundSpeed'.format(plot_dir));

        if(pbeta):
            physicalField = 'plasma_beta'
            slc = yt.SlicePlot(ds, 'y', ("gas", '{}'.format(physicalField))); 
            
            slc.annotate_streamlines( ("gas", "magnetic_field_z"), ("gas", "magnetic_field_x"), density = 3, color = 'white');
            slc.annotate_timestamp(text_args = {'color':'white'});
            slc.set_cmap(("gas", physicalField), "RdBu_r");
            slc.save('{}slice_pbeta_{}'.format(plot_dir, fileNumber))

            plasma_beta = magnetic_field_magnitude/(6.459534296116588e-7)
            # plamsa_beta = thermal_pressure/magnetic_pressure
            fig, ax = plt.subplots();
            ax.hist(plasma_beta);
            ax.set_xlabel(r'\log_{10} P beta');
            fig.savefig('{}pbeta_hist'.format(plot_dir));
'''
def stellar_projection_plot(nmin, nmax, saveName,file_read_directory_prefix, plot_dir):
    for i in range(nmin, nmax+1):
        print('plotting file number: ', i);
        fileNumber = i;
        zeroes = (4 - math.floor(np.log10(fileNumber)))*'0';
        fileNumber = '{}{}'.format(zeroes, fileNumber)
        saveName = saveName+'_{}'.format(fileNumber)      

        file_read_directory = '{}output_{}/'.format(file_read_directory_prefix, fileNumber);

        ds, fileNumber = open_disc_file(file_read_directory, fileNumber); 
        ds.add_particle_filter("stars_user");
        ds.add_particle_filter("darkMatter_user");
  
        print('Printing out the field list in the RAMSES dataset', ds.field_list);
        print('Particle types in the data: ', ds.particle_types_raw);

        ad = ds.all_data();

        particle_mass = ad[("stars_user", "particle_mass")];
        print('Particle count dict', ds.particle_type_counts);
        # print(particle_mass)

        fig, ax = plt.subplots();
        ax.hist(particle_mass);
        ax.set_ylabel('PDF');
        ax.set_xlabel('M');
        fig.savefig('particleMass_distribution');


        p = yt.ParticleProjectionPlot(ds,'z', [("stars_user", "particle_mass")], width = (100, "kpc"));
        p.set_cmap([("stars_user", "particle_mass")], 'copper')
        p.save('{}projected_stellar_map_{}'.format(plot_dir, fileNumber));
    return;

'''
Gas density projections of a single RAMSES output: face-on (along z, the disc axis)
and edge-on (along x, perpendicular to the disc axis). Two files are saved; yt appends
'_Projection_z_density' / '_Projection_x_density' to the name.
rmax is the maximum cylindrical radius, in kpc. Only gas inside the cylinder
R < rmax (spanning the full box height) enters the projection. In the edge-on view
the line of sight crosses the cylinder, so the column is limited to the chord length.
'''
def plot_gas_projection(fileNumber, file_read_directory_prefix, plot_dir, rmax):
    fileNumber = '{:05d}'.format(fileNumber);
    file_read_directory = '{}output_{}/'.format(file_read_directory_prefix, fileNumber);
    ds, fileNumber = open_disc_file(file_read_directory, fileNumber);

    center = ds.domain_center;
    rmax_code = ds.quan(rmax, 'kpc');
    half_height = ds.domain_width[2]/2;

    # Cylinder with axis along z, centred on the box centre.
    cylinder = ds.disk(center, [0, 0, 1], rmax_code, half_height);

    field = ("gas", "density");
    for axis in ['z', 'x']: # 'z' = face-on, 'x' = edge-on
        p = yt.ProjectionPlot(ds, axis, field, center = center, width = (2*rmax, 'kpc'),
                              data_source = cylinder);
        p.set_cmap(field, 'viridis');
        p.set_font_size(18);
        p.annotate_timestamp(text_args = {'color':'white'});
        p.save('{}projection_gas_{}'.format(plot_dir, fileNumber));
    return;

'''
Slices of the gas number density in particles per cm^3: the mid-plane (z = box centre,
face-on) and the plane x = box centre (edge-on, perpendicular to the disc). Two files
are saved; yt appends '_Slice_z_...' / '_Slice_x_...' to the name.
n = rho/(mu m_H), with mu = mean_mass (module-level default) unless given.
rmax (kpc) sets the half-width of the plotted square.
'''
def plot_gas_slice_ppcc(fileNumber, file_read_directory_prefix, plot_dir, rmax, mu = None):
    if mu is None:
        mu = mean_mass;
    fileNumber = '{:05d}'.format(fileNumber);

    file_read_directory = '{}output_{}/'.format(file_read_directory_prefix, fileNumber);
    ds, fileNumber = open_disc_file(file_read_directory, fileNumber);

    def _number_density_ppcc(field, data):
        return ds.arr(data["gas", "density"].v*const.avogadro_no/const.mean_molecular_mass, "cm**-3");

    ds.add_field(("gas", "number_density_ppcc"), function = _number_density_ppcc,
                 sampling_type = "local", units = "cm**-3", force_override = True);

    field = ("gas", "number_density_ppcc");
    for axis in ['z', 'x']: # 'z' = face-on, 'x' = edge-on
        s = yt.SlicePlot(ds, axis, field, center = ds.domain_center, width = (2*rmax, 'kpc'));
        # s.set_cmap(field, 'viridis'); 
        s.set_cmap(field, cmr.eclipse); 
        s.set_colorbar_label(field, r'$n\;(\rm cm^{-3})$');
        s.set_font_size(18);
        s.annotate_timestamp(text_args = {'color':'white'});
        s.save('{}slice_gas_ppcc_{}'.format(plot_dir, fileNumber));
    return;

'''
Radial profiles in the mid-plane (z = box centre, i.e. galactic z = 0) out to the
cylindrical radius rmax (kpc): temperature (K), number density (cm^-3), |B| (Gauss).
Grey points are individual AMR cells; the line is the mean in nbins radial bins.

Caveat: the cells lie on a Cartesian grid, so an annulus is only approximated by the
cells whose centres fall in it (the cell size sets a minimum useful bin width, and
inner bins contain few cells). Radial and azimuthal field components (e.g. v_R, v_phi,
B_R, B_phi) obtained by projecting the Cartesian components onto (R, phi) therefore
carry discretisation errors, largest near the centre and wherever the cells are coarse.

To do: Include the analytical radial profiles that I use in the initialisation procedure. 


'''

def plot_analytical_density_profile(ax, r):
    mg_disk = (8.593*(10)**9)*const.msun
    rd= 3.43218*units_length;
    zd = 0.343218*units_length;
    rho0 = mg_disk/(4*np.pi*(rd**2)*zd);
    rcgs = r*const.one_kpc;

    rho_r = rho0*np.exp(-rcgs/rd) #Safely assuming mid-plane.
    n_r = rho_r*const.avogadro_no/const.mean_molecular_mass

    ax.plot(r, n_r, linestyle = ':', color = 'green', alpha = 0.5);


    return ax;

def plot_midplane_radial_profiles(fileNumber, file_read_directory_prefix, plot_dir, rmax, mu = None, nbins = 50):
    if mu is None:
        mu = mean_mass;
    fileNumber = '{:05d}'.format(fileNumber);
    file_read_directory = '{}output_{}/'.format(file_read_directory_prefix, fileNumber);
    ds, fileNumber = open_disc_file(file_read_directory, fileNumber);


    center = ds.domain_center;
    slc = ds.slice('z', center[2]); # cells intersecting the plane z = z_centre

    x = (slc["gas", "x"] - center[0]).to('kpc').v;
    y = (slc["gas", "y"] - center[1]).to('kpc').v;
    r = np.sqrt(x**2 + y**2);

    temperature = (slc["gas", "pressure"].v*const.mean_molecular_mass)/(slc["gas", "density"].v*const.boltzmann_cgs*const.avogadro_no);
    n_ppcc = (slc["gas", "density"].v * const.avogadro_no/const.mean_molecular_mass)
    B_gauss = slc["gas", "magnetic_field_strength"].to('gauss').v;

    # Cylindrical azimuthal velocity v_phi = (x v_y - y v_x)/R, in km/s.
    vx = slc["gas", "velocity_x"].to('km/s').v;
    vy = slc["gas", "velocity_y"].to('km/s').v;
    v_phi = (x*vy - y*vx)/np.where(r > 0, r, np.nan);

    # ---- Radial binning ----------------------------------------------------
    # Goal: replace the cloud of per-cell values by one mean value per radial annulus.

    # Boolean mask: True for cells with R <= rmax. Cells outside are ignored everywhere below.
    inside = r <= rmax;

    # nbins+1 bin edges spanning [0, rmax], equally spaced.
    # Bin b is the annulus edges[b] <= R < edges[b+1].
    edges = np.linspace(0, rmax, nbins+1);

    # Bin centres (the average of each pair of neighbouring edges); used as the x-values of the line.
    r_mid = 0.5*(edges[1:] + edges[:-1]);

    # For every selected cell, find which bin it falls in.
    # np.digitize returns i such that edges[i-1] <= R < edges[i], with i = 1..nbins
    # (i = nbins+1 only for R == rmax exactly). Subtracting 1 turns it into a 0-based
    # bin index b, so that cell k lies in the annulus edges[b] <= R_k < edges[b+1].
    # 'which' has the same length as r[inside], i.e. it is aligned with q[inside] below.
    which = np.digitize(r[inside], edges) - 1;

    def binned_mean(q):
        # Keep only the cells inside rmax, so that q is aligned with 'which'.
        q = q[inside];
        # One output value per bin; NaN marks bins with no cells (matplotlib skips NaNs).
        out = np.full(nbins, np.nan);
        for b in range(nbins):
            # Boolean mask selecting the cells that belong to bin b.
            sel = which == b;
            if sel.any():
                # Unweighted arithmetic mean of the cell values in this annulus
                # (every cell counts equally, regardless of its size or volume).
                out[b] = q[sel].mean();
        return out;

    quantities = [(temperature, r'$T\;(\rm K)$', 'temperature'),
                  (n_ppcc, r'$n\;(\rm cm^{-3})$', 'density_ppcc'),
                  (B_gauss, r'$|B|\;(\rm G)$', 'B_gauss'),
                  (v_phi, r'$v_\phi\;(\rm km\,s^{-1})$', 'v_phi')];

    fig, axes = plt.subplots(4, 1, sharex = True, figsize = (6, 12));
    for ax, (q, label, name) in zip(axes, quantities):
        ax.scatter(r[inside], q[inside], s = 1, color = 'grey', alpha = 0.3);
        ax.plot(r_mid, binned_mean(q), color = 'k');
        if name != 'v_phi': # v_phi can be negative, keep it linear
            ax.set_yscale('log');
        ax.set_ylabel(label);


    plot_analytical_density_profile(axes[1], r[inside])


    axes[-1].set_xlabel(r'$R\;(\rm kpc)$');
    axes[-1].set_xlim(0, rmax);
    fig.tight_layout();


    fig.savefig('{}radial_profiles_midplane_{}'.format(plot_dir, fileNumber));
    plt.close(fig);
    return;


'''
Demonstration of Cartesian discretisation errors in cylindrical components.
Analytical field: A = |A| e_phi, |A| = 1, on an n x n grid covering [-1, 1]^2 (the origin
sits on a cell corner, as in a grid with an even number of cells).
As in a finite-volume code, the Cartesian components (A_x, A_y) stored in a cell are the
cell averages of the analytical field, while the unit vectors e_R, e_phi used for the
decomposition are evaluated at the cell centre:
    A_R   =  A_x cos(phi_c) + A_y sin(phi_c)
    A_phi = -A_x sin(phi_c) + A_y cos(phi_c)
Exactly, A_R = 0 and A_phi = 1. (If A_x, A_y were also point values at the cell centre,
the decomposition would be exact; the error comes from averaging over the cell.)
nsub x nsub sub-samples per cell are used to evaluate the cell average.
'''
def plot_cartesian_discretization_error(plot_dir = '', n = 1024, nsub = 8):
    L = 1.0; # half-width of the box
    dx = 2*L/n;
    xc = -L + (np.arange(n) + 0.5)*dx; # cell centres
    X, Y = np.meshgrid(xc, xc, indexing = 'ij');

    # Cell-averaged Cartesian components of e_phi = (-y, x)/R
    Ax = np.zeros((n, n));
    Ay = np.zeros((n, n));
    offsets = (np.arange(nsub) + 0.5)/nsub - 0.5; # sub-sample positions within a cell, in units of dx
    for ox in offsets:
        for oy in offsets:
            xs = X + ox*dx;
            ys = Y + oy*dx;
            Rs = np.hypot(xs, ys);
            Ax += -ys/Rs;
            Ay += xs/Rs;
    Ax /= nsub**2;
    Ay /= nsub**2;

    # Unit vectors at the cell centres
    phi = np.arctan2(Y, X);
    A_R = Ax*np.cos(phi) + Ay*np.sin(phi);
    A_phi = -Ax*np.sin(phi) + Ay*np.cos(phi);

    print('max |A_R| = {:.3e}, min A_phi = {:.5f}, max |A_phi - 1| = {:.3e}'.format(
          np.abs(A_R).max(), A_phi.min(), np.abs(A_phi - 1).max()));

    # The errors are confined to the few cells nearest the origin, so show only the central
    # ncrop x ncrop cells (cell edges visible), using the full range of the data as colour limits.
    ncrop = 40;
    sl = slice(n//2 - ncrop//2, n//2 + ncrop//2);
    A_R = A_R[sl, sl];
    A_phi = A_phi[sl, sl];
    h = ncrop//2*dx;
    extent = (-h, h, -h, h);
    # Symmetric log colour scale: keeps the sign of the errors A_R - 0 and A_phi - 1.
    # Values with |error| < 1e-6 of the maximum are shown on a linear scale around zero.
    err_R = A_R;
    err_phi = A_phi - 1;

    fig, axes = plt.subplots(1, 2, figsize = (12, 5));

    for ax, err, title in [(axes[0], err_R, r'$A_R$ (exact: 0)'),
                           (axes[1], err_phi, r'$A_\phi - 1$ (exact: 0)')]:
        vmax = np.abs(err).max();
        norm = mpl.colors.SymLogNorm(linthresh = 1e-6*vmax, vmin = -vmax, vmax = vmax, base = 10);
        im = ax.imshow(err.T, origin = 'lower', extent = extent, cmap = 'RdBu_r',
                       norm = norm, interpolation = 'nearest');
        ax.set_title(title);
        fig.colorbar(im, ax = ax);

    for ax in axes:
        ax.set_xlabel(r'$x$');
        ax.set_ylabel(r'$y$');
    fig.tight_layout();
    fig.savefig('{}cartesian_discretization_error'.format(plot_dir));
    plt.close(fig);
    return;

def plot_analytical_density_vertical_profile(ax, z, r0):
    mg_disk = (8.593*(10)**9)*const.msun
    rd= 3.43218*units_length;
    zd = 0.343218*units_length;
    rho0 = mg_disk/(4*np.pi*(rd**2)*zd);
    zcgs = z*const.one_kpc;
    r0cgs = r0*const.one_kpc;


    rho_z = rho0*np.exp(-np.abs(zcgs)/zd)*np.exp(-r0cgs/rd) #Safely assuming mid-plane.
    n_z = rho_z*const.avogadro_no/const.mean_molecular_mass

    ax.plot(z, n_z, linestyle = ':', color = 'green', alpha = 0.5);


    return ax;
'''
Vertical (z) profiles of temperature (K), number density (cm^-3) and |B| (Gauss).
Uses the cells in the annulus |R - R0| <= dR/2 (kpc), |z| <= zmax (kpc), centred on the
box centre, with the disc axis along z. A single exact radius contains no cells, so the
profile is taken over a thin annulus of width dR around R0 (default 3 kpc, 2.75-3.25 kpc).
Grey points are individual AMR cells; the line is the mean in nbins equal-width bins in
signed z.
'''
def plot_vertical_profiles(fileNumber, file_read_directory_prefix, plot_dir, zmax, R0 = 3.0, dR = 0.5, mu = None, nbins = 50):
    if mu is None:
        mu = mean_mass;
    fileNumber = '{:05d}'.format(fileNumber);
    file_read_directory = '{}output_{}/'.format(file_read_directory_prefix, fileNumber);
    ds, fileNumber = open_disc_file(file_read_directory, fileNumber);

    center = ds.domain_center;
    # Cylinder with axis along z, just large enough to contain the annulus;
    # the last argument is the half-height. The annulus itself is cut out below.
    cyl = ds.disk(center, [0, 0, 1], ds.quan(R0 + dR, 'kpc'), ds.quan(zmax, 'kpc'));

    x = (cyl["gas", "x"] - center[0]).to('kpc').v;
    y = (cyl["gas", "y"] - center[1]).to('kpc').v;
    z = (cyl["gas", "z"] - center[2]).to('kpc').v;
    R = np.sqrt(x**2 + y**2); # cylindrical radius of each cell centre

    temperature = (cyl["gas", "pressure"].v*const.mean_molecular_mass)/(cyl["gas", "density"].v*const.boltzmann_cgs*const.avogadro_no);
    n_ppcc = (cyl["gas", "density"]*const.avogadro_no/const.mean_molecular_mass)
    B_gauss = cyl["gas", "magnetic_field_strength"].to('gauss').v;

    # ---- Vertical binning --------------------------------------------------
    # Same algorithm as in plot_midplane_radial_profiles, but in signed z over [-zmax, zmax].
    # Keep only cells whose centres lie in the annulus |R - R0| <= dR/2 and in |z| <= zmax
    # (the disk selects by cell overlap, so it contains cells outside both limits).
    inside = (np.abs(R - R0) <= dR/2) & (np.abs(z) <= zmax);
    edges = np.linspace(-zmax, zmax, nbins+1);
    z_mid = 0.5*(edges[1:] + edges[:-1]);
    which = np.digitize(z[inside], edges) - 1; # 0-based bin index of each selected cell

    def binned_mean(q):
        q = q[inside];
        out = np.full(nbins, np.nan); # NaN for empty bins
        for b in range(nbins):
            sel = which == b;
            if sel.any():
                out[b] = q[sel].mean(); # unweighted mean over the cells in the bin
        return out;

    quantities = [(temperature, r'$T\;(\rm K)$'),
                  (n_ppcc, r'$n\;(\rm cm^{-3})$'),
                  (B_gauss, r'$|B|\;(\rm G)$')];

    fig, axes = plt.subplots(3, 1, sharex = True, figsize = (6, 9));
    for ax, (q, label) in zip(axes, quantities):
        ax.scatter(z[inside], q[inside], s = 1, color = 'grey', alpha = 0.3);
        ax.plot(z_mid, binned_mean(q), color = 'k');
        ax.set_yscale('log');
        ax.set_ylabel(label);

    ax = plot_analytical_density_vertical_profile(axes[1], z[inside], R0)

    axes[-1].set_xlabel(r'$z\;(\rm kpc)$');
    axes[-1].set_xlim(-zmax, zmax);
    axes[0].set_title(r'$R = %.2f \pm %.2f\;\rm kpc$'%(R0, dR/2));
    fig.tight_layout();
    fig.savefig('{}vertical_profiles_R{}kpc_{}.png'.format(plot_dir, R0, fileNumber));
    plt.close(fig);
    return;

'''
Mid-plane (z = box centre) slices of the magnetic field in cylindrical components
(B_R, B_phi, B_z), one figure per component, in microGauss. The slice is sampled onto an
npix x npix grid covering the square |x|, |y| <= rmax (kpc) about the box centre. The
components come from yt's derived fields magnetic_field_cylindrical_{radius,theta,z}, which
use the 'center' and 'normal' field parameters (box centre, z axis). Saved as
'{plot_dir}midplane_{BR,Bphi,Bz}_{fileNumber}.png'.
'''
def plot_midplane_B_cylindrical(fileNumber, file_read_directory_prefix, plot_dir, rmax, npix = 512):
    fileNumber = '{:05d}'.format(fileNumber);
    file_read_directory = '{}output_{}/'.format(file_read_directory_prefix, fileNumber);
    ds, fileNumber = open_disc_file(file_read_directory, fileNumber);

    center = ds.domain_center;
    slc = ds.slice('z', center[2]);
    slc.set_field_parameter('center', center);
    slc.set_field_parameter('normal', [0, 0, 1]);
    half = ds.quan(rmax, 'kpc').to('code_length');
    bounds = (center[0] - half, center[0] + half, center[1] - half, center[1] + half);
    frb = yt.FixedResolutionBuffer(slc, bounds, (npix, npix));

    components = [("magnetic_field_cylindrical_radius", r'B_R', 'BR'),
                  ("magnetic_field_cylindrical_theta", r'B_\phi', 'Bphi'),
                  ("magnetic_field_cylindrical_z", r'B_z', 'Bz')];
    for field, label, name in components:
        B = frb["gas", field].to('G').v;
        vmax = np.nanpercentile(np.abs(B), 99);
        fig, ax = plt.subplots();
        im = ax.imshow(B, origin = 'lower', extent = (-rmax, rmax, -rmax, rmax),
                       cmap = 'RdBu_r', vmin = -vmax, vmax = vmax);
        ax.set_xlabel(r'$x\;(\rm kpc)$');
        ax.set_ylabel(r'$y\;(\rm kpc)$');
        clb = fig.colorbar(im);
        clb.set_label(r'$%s\;(\mu\rm G)$'%label);
        fig.tight_layout();
        fig.savefig('{}midplane_{}_{}.png'.format(plot_dir, name, fileNumber));
        plt.close(fig);
    return;

'''
Mid-plane (z = box centre) slice of yt's magnetic field divergence.
For RAMSES in yt, this field is built from face-centred components as
divB = sum_axis(B_axis_right - B_axis_left)/dx, so it is a discrete finite-volume
divergence evaluated at cell centres.
Saved as '{plot_dir}midplane_divB_{fileNumber}.png'.
'''
def plot_midplane_divB_slice(fileNumber, file_read_directory_prefix, plot_dir, rmax, npix = 512):
    fileNumber = '{:05d}'.format(fileNumber);
    file_read_directory = '{}output_{}/'.format(file_read_directory_prefix, fileNumber);
    ds, fileNumber = open_disc_file(file_read_directory, fileNumber);

    center = ds.domain_center;
    slc = ds.slice('z', center[2]);
    half = ds.quan(rmax, 'kpc').to('code_length');
    bounds = (center[0] - half, center[0] + half, center[1] - half, center[1] + half);
    frb = yt.FixedResolutionBuffer(slc, bounds, (npix, npix));

    divB = frb["gas", "magnetic_field_divergence"].to('gauss/cm').v;
    vmax = np.nanpercentile(np.abs(divB), 99);
    if (not np.isfinite(vmax)) or (vmax <= 0):
        vmax = np.nanmax(np.abs(divB));
    if (not np.isfinite(vmax)) or (vmax <= 0):
        vmax = 1e-30;

    norm = mpl.colors.SymLogNorm(linthresh = 1e-3*vmax, vmin = -vmax, vmax = vmax, base = 10);

    fig, ax = plt.subplots();
    im = ax.imshow(divB, origin = 'lower', extent = (-rmax, rmax, -rmax, rmax),
                   cmap = 'RdBu_r', norm = norm);
    ax.set_xlabel(r'$x\;(\rm kpc)$');
    ax.set_ylabel(r'$y\;(\rm kpc)$');
    clb = fig.colorbar(im);
    clb.set_label(r'$\nabla\cdot\mathbf{B}\;(\rm G\,cm^{-1})$');
    fig.tight_layout();
    fig.savefig('{}midplane_divB_{}.png'.format(plot_dir, fileNumber));
    plt.close(fig);
    return;

def plot_particle_mass_projection(nmin, nmax, saveName, file_read_directory_prefix,
                                 plot_dir, axis='z', width=(100, 'kpc')):
    """Save a yt projection of all particles' mass for each RAMSES output.

    The projected field is particle mass, displayed in solar masses. ``width``
    is passed to yt and defaults to a 100 kpc-wide field of view.
    """
    field = ('all', 'particle_mass')
    for output_number in range(nmin, nmax + 1):
        file_number = '{:05d}'.format(output_number)
        output_directory = '{}output_{}/'.format(file_read_directory_prefix, file_number)
        ds, _ = open_disc_file(output_directory, file_number)

        for projection_axis in dict.fromkeys((axis, 'x')):
            projection = yt.ParticleProjectionPlot(ds, projection_axis, [field], width=width)
            projection.set_unit(field, 'Msun')
            output_name = '{}{}_{}_{}'.format(plot_dir, saveName, file_number, projection_axis)
            print('Saving particle mass projection:', output_name)
            projection.save(output_name)

    return

# ====================================
# Some constants
# =====================================

avogadro_no = 6.0221e23;
mean_mass = const.mean_molecular_mass;
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
'''
# nmin = 1;
# nmax = 11;
file_read_directory_prefix = '/home/ra6154ar/ramses_updated/simulations/isolated_galaxy/'
plot_dir = '/home/ra6154ar/ramses_updated/simulations/isolated_galaxy/';

fileNumber = 1
rmax = 25 #kpc


# plot_gas_slice_ppcc(fileNumber, file_read_directory_prefix, plot_dir, rmax, mu = mean_mass);

# plot_midplane_radial_profiles(fileNumber, file_read_directory_prefix, plot_dir, rmax, mu = mean_mass, nbins = 50);

# plot_vertical_profiles(fileNumber, file_read_directory_prefix, plot_dir, zmax=10, mu = mean_mass, nbins = 50)

# plot_midplane_B_cylindrical(fileNumber, file_read_directory_prefix, plot_dir, rmax)

# plot_midplane_divB_slice(fileNumber, file_read_directory_prefix, plot_dir, rmax)

plot_particle_mass_projection(1, 1, 'particlePlot',file_read_directory_prefix, 
                                 plot_dir, axis='z', width=(rmax, 'kpc'))

# plot_cartesian_discretization_error(plot_dir = plot_dir, n = 1024, nsub = 8)
