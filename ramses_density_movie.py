'''
Movies of the gas number-density evolution in two slices through the box centre:
the mid-plane (z slice, face-on) and the plane x = box centre (x slice, edge-on).

Reads every output_XXXXX directory found in the simulation directory, makes one slice
per output and per orientation (same style as plot_gas_slice_ppcc in
ramses_projection_slice.py), and stitches the frames into one mp4 per orientation with
ffmpeg. The frames are kept in <frame_dir> as frame_{face_on,edge_on}_XXXXX.png.

Usage:
    python3 ramses_density_movie.py
    python3 ramses_density_movie.py --rmax 25 --fps 10 --vmin 1e-6 --vmax 1e2
If --vmin/--vmax are not given, the colour limits are taken from the first snapshot (for
each orientation separately) so that all frames of a movie share one fixed colour scale.
'''
import argparse
import glob
import os
import re
import shutil
import subprocess

import numpy as np
import matplotlib as mpl
import yt
import cmasher as cmr
import rapack.constants as const

# Importing rapack (rapack/matplotlib_settings.py) sets text.usetex = True, which needs a
# LaTeX installation. Turn it off again (must come after the rapack import).
mpl.rcParams['text.usetex'] = False;

default_dir = '/home/ra6154ar/ramses_updated/simulations/isolated_galaxy/';


def find_outputs(sim_dir):
    # Returns a list of (output number, info file path), sorted by output number.
    outputs = [];
    for d in glob.glob(os.path.join(sim_dir, 'output_*')):
        m = re.fullmatch(r'output_(\d+)', os.path.basename(d));
        if m is None:
            continue;
        num = m.group(1);
        info = os.path.join(d, 'info_{}.txt'.format(num));
        if os.path.isfile(info):
            outputs.append((int(num), info));
    outputs.sort();
    return outputs;


def add_number_density(ds):
    def _number_density_ppcc(field, data):
        return ds.arr(data["gas", "density"].v*const.avogadro_no/const.mean_molecular_mass, "cm**-3");
    ds.add_field(("gas", "number_density_ppcc"), function = _number_density_ppcc,
                 sampling_type = "local", units = "cm**-3", force_override = True);
    return ("gas", "number_density_ppcc");


# slice normal axis -> label used in the file names
orientations = {'z': 'face_on', 'x': 'edge_on'};


def make_movie(sim_dir, frame_dir, movie_base, rmax, fps, vmin = None, vmax = None):
    outputs = find_outputs(sim_dir);
    if len(outputs) == 0:
        raise SystemExit('No outputs found in {}'.format(sim_dir));
    print('Found {} outputs: {} ... {}'.format(len(outputs), outputs[0][0], outputs[-1][0]));

    os.makedirs(frame_dir, exist_ok = True);
    for old in glob.glob(os.path.join(frame_dir, 'frame_*.png')):
        os.remove(old);

    limits = {axis: (vmin, vmax) for axis in orientations}; # per-orientation colour limits

    for iframe, (num, info) in enumerate(outputs):
        print('plotting output {} ({}/{})'.format(num, iframe + 1, len(outputs)));
        ds = yt.load(info);
        field = add_number_density(ds);

        for axis, label in orientations.items():
            lo, hi = limits[axis];
            if lo is None or hi is None:
                # Fixed colour limits from the first snapshot's slice in this orientation.
                vals = ds.slice(axis, ds.domain_center["xyz".index(axis)])[field].v;
                vals = vals[vals > 0];
                lo = vals.min() if lo is None else lo;
                hi = vals.max() if hi is None else hi;
                limits[axis] = (lo, hi);
                print('{} colour limits: {:.3e} - {:.3e} cm^-3'.format(label, lo, hi));

            s = yt.SlicePlot(ds, axis, field, center = ds.domain_center, width = (2*rmax, 'kpc'));
            s.set_cmap(field, cmr.eclipse);
            s.set_zlim(field, lo, hi);
            s.set_colorbar_label(field, r'$n\;(\rm cm^{-3})$');
            s.set_font_size(18);
            s.annotate_timestamp(text_args = {'color':'white'});

            saved = s.save(os.path.join(frame_dir, 'tmp_frame'));
            shutil.move(saved[0], os.path.join(frame_dir, 'frame_{}_{:05d}.png'.format(label, iframe)));

    ffmpeg = shutil.which('ffmpeg');
    if ffmpeg is None:
        print('ffmpeg not found: frames are in {}'.format(frame_dir));
        return;
    for label in orientations.values():
        movie_name = '{}_{}.mp4'.format(movie_base, label);
        # Even dimensions are required by yuv420p; the pad filter guarantees them.
        cmd = [ffmpeg, '-y', '-framerate', str(fps), '-i', os.path.join(frame_dir, 'frame_{}_%05d.png'.format(label)),
               '-vf', 'pad=ceil(iw/2)*2:ceil(ih/2)*2', '-c:v', 'libx264', '-pix_fmt', 'yuv420p', movie_name];
        subprocess.run(cmd, check = True);
        print('wrote {}'.format(movie_name));
    return;


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description = 'Mid-plane density movie from RAMSES outputs');
    parser.add_argument('--sim_dir', default = default_dir);
    parser.add_argument('--frame_dir', default = None, help = 'default: <sim_dir>/movie_frames');
    parser.add_argument('--movie', default = None,
                        help = 'movie base name; _face_on.mp4 / _edge_on.mp4 are appended. default: <sim_dir>/density');
    parser.add_argument('--rmax', type = float, default = 25.0, help = 'half-width of the plot, kpc');
    parser.add_argument('--fps', type = int, default = 5);
    parser.add_argument('--vmin', type = float, default = None, help = 'lower colour limit, cm^-3');
    parser.add_argument('--vmax', type = float, default = None, help = 'upper colour limit, cm^-3');
    args = parser.parse_args();

    frame_dir = args.frame_dir or os.path.join(args.sim_dir, 'movie_frames');
    movie = args.movie or os.path.join(args.sim_dir, 'density');
    if movie.endswith('.mp4'):
        movie = movie[:-4];
    make_movie(args.sim_dir, frame_dir, movie, args.rmax, args.fps, args.vmin, args.vmax);
