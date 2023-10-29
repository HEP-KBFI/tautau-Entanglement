#!/usr/bin/env python3

# Example usage:
#
# plot_entanglement_vars.py -i ~/Entanglement/analysis/SuperKEKB/2023Oct22_1M -o ~/Entanglement/analysis/SuperKEKB/2023Oct22_1M/json_plots -E png
#
# The output directory is created automatically with plots in case it doesn't exist. Both the entanglement variable and its significance are
# plotted. You can change that by supplying "-S only" if you want only the plots for significance or by supplying "-S excluded" if you want
# plots only for the entanglement variable itself. By default only Rchsh is plotted, but you can also change that via -e option.
# Option -E can be ignored if you want plots with pdf and png extensions. Finally, plots can be created for the following cases:
# 1) "individual" -- if you want a plot for specific initialization mode, tau decay mode and spin analyzer method;
# 2) "init_method" -- if you want a plot that groups together all initialization modes for a particular tau decay mode and spin analyzer method;
# 3) "decay_mode" -- if you want a plot that groups together all tau decay modes for a particular initialization mode and spin analyzer method;
# 4) "spin_analyzer" -- if you want a plot that groups together all spin analyzer methods for a particular initialization and tau decay mode.
# All four options are exercised by default. Options 2-4 allow to compare initialization modes, tau decay modes and spin analyzer methods.

from TauAnalysis.Entanglement.tools.aux import ENTANGLEMENT_VARS

import argparse
import json
import matplotlib.pyplot as plt
import re
import os

import matplotlib as mpl
CMS = {
  "mathtext.default": "regular",
  "figure.figsize": (10.0, 10.0),
  "font.size": 26,
  "axes.labelsize": "medium",
  "axes.unicode_minus": False,
  "xtick.labelsize": "small",
  "ytick.labelsize": "small",
  "legend.fontsize": "small",
  "legend.handlelength": 1.5,
  "legend.borderpad": 0.5,
  "legend.frameon": True,
  "xtick.direction": "in",
  "xtick.major.size": 12,
  "xtick.minor.size": 6,
  "xtick.major.pad": 6,
  "xtick.top": True,
  "xtick.major.top": True,
  "xtick.major.bottom": True,
  "xtick.minor.top": True,
  "xtick.minor.bottom": True,
  "xtick.minor.visible": True,
  "ytick.direction": "in",
  "ytick.major.size": 12,
  "ytick.minor.size": 6.0,
  "ytick.right": True,
  "ytick.major.left": True,
  "ytick.major.right": True,
  "ytick.minor.left": True,
  "ytick.minor.right": True,
  "ytick.minor.visible": True,
  "grid.alpha": 0.8,
  "grid.linestyle": ":",
  "axes.linewidth": 2,
  "savefig.transparent": False,
  "xaxis.labellocation": "right",
  "yaxis.labellocation": "top",
  "legend.title_fontsize": "medium",
}

for k, v in CMS.items():
  if k in mpl.rcParams:
    mpl.rcParams[k] = v

fn_rgx = re.compile(
  r"analyzeEntanglementNtuple_(?P<sample_name>\w+)_(?P<mode>\w+)Mode_(?P<axis>\w+)Axis_(?P<dm>\w+)DecayMode_" + \
  r"by_(?P<spin_analyzer>\w+)_absCosTheta(?P<cut>.+).json"
)

# key = entanglement variable
# values = lower threshold for observing the QM effects; symbol in latex
CENTRAL_CHOICES = [ 'nominal', 'median' ]
AXIS_CHOICES = [ 'beam', 'higgs' ]

SPIN_ANALYZERS = {
  'asymmetry'          : 'FB asymm.',
  'differentialXsec1d' : '1d distr.',
  'differentialXsec2d' : '2d distr.',
  'mlfit'              : 'ML-fit',
  'summation'          : 'Exp. value',
}
DECAY_MODES = {
  'pi_pi'   : { 'label' : r'$\pi^+\pi^-$',       'marker' : 'o' },
  'pi_rho'  : { 'label' : r'$\pi^\pm\rho^\mp$',  'marker' : 'v' },
  'pi_a1'   : { 'label' : r'$\pi^\pm a_1^\mp$',  'marker' : '^' },
  'rho_rho' : { 'label' : r'$\rho^+\rho^-$',     'marker' : 's' },
  'rho_a1'  : { 'label' : r'$\rho^\pm a_1^\mp$', 'marker' : 'X' },
  'a1_a1'   : { 'label' : r'$a_1^+a_1^-$',       'marker' : 'd' },
  'comb'    : { 'label' : 'Combination',         'marker' : '*' },
#  'had_had' : { 'label' : 'All hadronic',        'marker' : 'r$\club$' },
}
INIT_METHODS = {
  'gen'      : 'MC-truth level',
#  'startPos' : 'Analytic',
  'kinFit'   : 'KF after smearing',
}
PLOT_CHOICES = [ 'individual', 'init_method', 'decay_mode', 'spin_analyzer' ]
EXTENSIONS = [ 'pdf', 'png' ]

DEFAULT_CENTRAL = 'median'
DEFAULT_AXIS = 'beam'
PLOT_CHOICES_DEFAULT = [ "decay_mode" ]

def str_to_float(s):
  return float(s.replace('p', '.'))

def comp_significance(src, entanglement_var, central):
  assert(entanglement_var in src)
  central_value = src[entanglement_var][central]
  error = src[entanglement_var]['error']
  threshold = ENTANGLEMENT_VARS[entanglement_var]['threshold']
  return (central_value - threshold) / error if (error > 0 and central_value > threshold) else 0

def read_data(input_dirs):
  data = {}
  for input_dir in input_dirs:
    for fn in os.listdir(input_dir):
      rgx_match = fn_rgx.match(fn)
      if not rgx_match:
        continue
      sample_name = rgx_match.group('sample_name')
      if sample_name not in data:
        data[sample_name] = {}
      axis = rgx_match.group('axis')
      if axis not in data[sample_name]:
        data[sample_name][axis] = {}
      mode = rgx_match.group('mode')
      if mode not in data[sample_name][axis]:
        data[sample_name][axis][mode] = {}
      dm = rgx_match.group('dm')
      if dm not in data[sample_name][axis][mode]:
        data[sample_name][axis][mode][dm] = {}
      spin_analyzer = rgx_match.group('spin_analyzer')
      if spin_analyzer not in data[sample_name][axis][mode][dm]:
        data[sample_name][axis][mode][dm][spin_analyzer] = {}
      cut = rgx_match.group('cut')
      if cut in data[sample_name][axis][mode][dm][spin_analyzer]:
        raise ValueError(f"Duplicate found: {fn}")
      with open(os.path.join(input_dir, fn), 'r') as f:
        json_data = json.load(f)
      for entanglement_var in ENTANGLEMENT_VARS:
        for central in CENTRAL_CHOICES:
          json_data[entanglement_var][f'{central}_significance'] = comp_significance(json_data, entanglement_var, central)
      data[sample_name][axis][mode][dm][spin_analyzer][cut] = json_data
  for sample_name in data:
    for axis in data[sample_name]:
      for mode in data[sample_name][axis]:
        dms, spin_analyzers, cuts = [], [], []
        for dm in data[sample_name][axis][mode]:
          if dm not in dms:
            dms.append(dm)
          for spin_analyzer in data[sample_name][axis][mode][dm]:
            if spin_analyzer not in spin_analyzers:
              spin_analyzers.append(spin_analyzer)
            for cut in data[sample_name][axis][mode][dm][spin_analyzer]:
              if cut not in cuts:
                cuts.append(cut)
        for dm in dms:
          for spin_analyzer in spin_analyzers:
            if spin_analyzer not in data[sample_name][axis][mode][dm]:
              print(f"Warning: no data found for sample_name={sample_name}, axis={axis}, mode={mode}, dm={dm} spin_analyzer={spin_analyzer}")
              continue
            for cut in cuts:
              if cut not in data[sample_name][axis][mode][dm][spin_analyzer]:
                print(f"Warning: no data found for sample_name={sample_name}, axis={axis}, mode={mode}, dm={dm} spin_analyzer={spin_analyzer} cut={cut}")
        data[sample_name][axis][mode]['comb'] = {}
        for spin_analyzer in spin_analyzers:
          data[sample_name][axis][mode]['comb'][spin_analyzer] = {}
          for cut in cuts:
            data[sample_name][axis][mode]['comb'][spin_analyzer][cut] = {}
            for entanglement_var in ENTANGLEMENT_VARS:
              data[sample_name][axis][mode]['comb'][spin_analyzer][cut][entanglement_var] = {}
              for central in CENTRAL_CHOICES:
                label = f'{central}_significance'
                data[sample_name][axis][mode]['comb'][spin_analyzer][cut][entanglement_var][label] = sum(
                  data[sample_name][axis][mode][dm][spin_analyzer][cut][entanglement_var][label]**2 for dm in dms if \
                  dm != 'had_had' and \
                  spin_analyzer in data[sample_name][axis][mode][dm] and \
                  cut in data[sample_name][axis][mode][dm][spin_analyzer]
                )**0.5
  return data

def plot(data, entanglement_var, is_significance, central = DEFAULT_CENTRAL, title = '', legend_title = '', output_fns = None):
  assert(central in CENTRAL_CHOICES)

  markersize = 10
  linewidth = 2
  plt.figure(figsize = (10, 8), dpi = 100)
  plot_threshold = False
  threshold = ENTANGLEMENT_VARS[entanglement_var]['threshold']
  max_point = (-1, -999)
  for label, src in data.items():
    cosTheta_bins = src.keys()
    plot_data = list(sorted(
      [ (str_to_float(cosTheta_bin), src[cosTheta_bin]) for cosTheta_bin in cosTheta_bins ],
      key = lambda pair: pair[0],
    ))
    xcoords = [ pair[0] for pair in plot_data ]
    ycoords = [ pair[1][entanglement_var][f'{central}_significance' if is_significance else central] for pair in plot_data ]
    marker = 'o'
    for dm in DECAY_MODES:
      if DECAY_MODES[dm]['label'] == label:
        marker = DECAY_MODES[dm]['marker']
    if is_significance:
      plt.plot(xcoords, ycoords, marker = marker, markersize = markersize, lw = linewidth, label = label)
    else:
      if any(yval < threshold for yval in ycoords):
        plot_threshold = True
      yerr = [ pair[1][entanglement_var]['error'] for pair in plot_data ]
      p = plt.scatter(xcoords, ycoords, marker = marker, s = markersize * 5, label = label)
      plt.errorbar(xcoords, ycoords, yerr, fmt = marker, capsize = 4, color = p.get_facecolor())
    max_point_it = sorted(zip(xcoords, ycoords), key = lambda xy: xy[1], reverse = True)[0]
    if max_point_it[1] > max_point[1]:
      max_point = max_point_it
  # if max_point[1] > 0:
  #   plt.plot(max_point[0], max_point[1], 'o', ms = 30, mec = 'lime', mfc = 'none', mew = 2)
  if plot_threshold:
    plt.axhline(y = threshold, ls = '--', color = 'black', lw = 2)
  plt.grid(True)
  plt.xlim(0, 1)
  plt.xlabel(r'Upper limit on $\left|\cos\theta^*\right|$')
  ylabel = f"${ENTANGLEMENT_VARS[entanglement_var]['symbol']}$"
  if is_significance:
    ylabel = f'Significance of {ylabel}'
  plt.ylabel(ylabel)
  if title:
    plt.title(title, fontsize = 20)
  if len(data) > 1:
    legend = plt.legend(
      title = legend_title,
      bbox_to_anchor = (0, 1, 1, 0),
      loc = "lower left",
      mode = "expand",
      fontsize = 18,
      ncol = max(len(data) // 2, 1),
    )
    plt.setp(legend.get_title(), fontsize = 'x-small')
  if output_fns:
    for output_fn in output_fns:
      plt.savefig(output_fn, bbox_inches = 'tight')
      print(f'Saved figure to {output_fn}')
  #plt.show()
  plt.close()

def has_data(src, init_method, decay_mode, spin_analyzer):
  return init_method in src and decay_mode in src[init_method] and spin_analyzer in src[init_method][decay_mode]

def construct_fns(entanglement_var, is_significance, output_dir, init_method = None, decay_mode = None, spin_analyzer = None, exts = None):
  assert(exts)
  entanglement_varName = entanglement_var
  if is_significance:
    entanglement_varName += 'Significance'
  baseName = f"plot_{entanglement_varName}_"
  mode_str = 'Mode'
  dm_str = 'DecayMode'
  analyzer_str = 'SpinAnalyzer'
  if init_method and decay_mode and spin_analyzer:
    baseName += f"{'_'.join([ init_method + mode_str, decay_mode + dm_str, spin_analyzer + analyzer_str ])}"
  elif init_method and decay_mode:
    baseName += f"{'_'.join([ init_method + mode_str, decay_mode + dm_str, 'by' + analyzer_str ])}"
  elif init_method and spin_analyzer:
    baseName += f"{'_'.join([ init_method + mode_str, spin_analyzer + analyzer_str, 'by' + dm_str ])}"
  elif decay_mode and spin_analyzer:
    baseName += f"{'_'.join([ decay_mode + dm_str, spin_analyzer + analyzer_str, 'by' + mode_str ])}"
  else:
    raise ValueError("Need to specify at least two of the three following arguments: method, decay_mode, spin_analyzer")
  return [ os.path.join(output_dir, f'{baseName}.{ext}') for ext in exts ]

def plot_by(data, entanglement_var, is_significance, init_method = None, decay_mode = None, spin_analyzer = None, central = DEFAULT_CENTRAL,
            output_dir = '', exts = None, show_title = False):
  if init_method:
    assert(init_method in INIT_METHODS)
  if decay_mode:
    assert(decay_mode in DECAY_MODES)
  if spin_analyzer:
    assert(spin_analyzer in SPIN_ANALYZERS)

    assert(os.path.isdir(output_dir))
    assert(exts)

  if decay_mode == "comb" and not is_significance:
    return False
  if init_method and decay_mode and spin_analyzer:
    if not has_data(data, init_method, decay_mode, spin_analyzer):
      print(f'Warning: could not find data for method = {init_method}, decay mode = {decay_mode} and spin analyzer = {spin_analyzer}')
      return False
    title = f"{INIT_METHODS[init_method]}, {DECAY_MODES[decay_mode]['label']} decay mode, {SPIN_ANALYZERS[spin_analyzer]} method" if show_title else ''
    src = { title : data[init_method][decay_mode][spin_analyzer] }
    fns = construct_fns(entanglement_var, is_significance, output_dir, init_method, decay_mode, spin_analyzer, exts)
    plot(src, entanglement_var, is_significance, central, title = title, output_fns = fns)
  elif init_method and decay_mode:
    title = f"{INIT_METHODS[init_method]}, {DECAY_MODES[decay_mode]['label']} decay mode" if show_title else ''
    src = {}
    for spin_analyzer in SPIN_ANALYZERS:
      if has_data(data, init_method, decay_mode, spin_analyzer):
        src[SPIN_ANALYZERS[spin_analyzer]] = data[init_method][decay_mode][spin_analyzer]
    if not src:
      print(f'Warning: could not find data for method = {init_method} and decay mode = {decay_mode}')
      return False
    fns = construct_fns(entanglement_var, is_significance, output_dir, init_method, decay_mode, None, exts)
    plot(src, entanglement_var, is_significance, central, title = title, legend_title = 'Spin analyzer' if show_title else '', output_fns = fns)
  elif init_method and spin_analyzer:
    title = f'{INIT_METHODS[init_method]}, {SPIN_ANALYZERS[spin_analyzer]} method' if show_title else ''
    src = {}
    for decay_mode in DECAY_MODES:
      if decay_mode == 'comb' and not is_significance:
        continue
      if has_data(data, init_method, decay_mode, spin_analyzer):
        src[DECAY_MODES[decay_mode]['label']] = data[init_method][decay_mode][spin_analyzer]
    if not src:
      print(f'Warning: could not find data for method = {init_method} and spin analyzer = {spin_analyzer}')
      return False
    fns = construct_fns(entanglement_var, is_significance, output_dir, init_method, None, spin_analyzer, exts)
    plot(src, entanglement_var, is_significance, central, title = title, legend_title = 'Decay mode' if show_title else '', output_fns = fns)
  elif decay_mode and spin_analyzer:
    title = f"{DECAY_MODES[decay_mode]['label']} decay mode, {SPIN_ANALYZERS[spin_analyzer]} method" if show_title else ''
    src = {}
    for init_method in INIT_METHODS:
      if has_data(data, init_method, decay_mode, spin_analyzer):
        src[INIT_METHODS[init_method]] = data[init_method][decay_mode][spin_analyzer]
    if not src:
      print(f'Warning: could not find data for decay mode = {decay_mode} and spin analyzer = {spin_analyzer}')
      return False
    fns = construct_fns(entanglement_var, is_significance, output_dir, None, decay_mode, spin_analyzer, exts)
    plot(src, entanglement_var, is_significance, central, title = title, legend_title = 'Starting point' if show_title else '', output_fns = fns)
  else:
    raise ValueError("Need to specify at least two of the three following arguments: method, decay_mode, spin_analyzer")
  return True

def plot_all(data, entanglement_vars, significance, by = None, exts = None, central = DEFAULT_CENTRAL, output_dir = '', show_title = False):
  if not by:
    by = PLOT_CHOICES
  else:
    assert(all(plot_choice in PLOT_CHOICES for plot_choice in by))
  if not exts:
    exts = EXTENSIONS
  args = []
  assert(all(entanglement_var in ENTANGLEMENT_VARS for entanglement_var in entanglement_vars))

  assert(output_dir)
  output_dir_abspath = os.path.abspath(output_dir)
  if not os.path.isdir(output_dir_abspath):
    os.makedirs(output_dir_abspath)
    print(f'Created directory {output_dir_abspath}')

  if significance == 'only':
    significance_opts = [ True ]
  elif significance == 'excluded':
    significance_opts = [ False ]
  elif significance == 'both':
    significance_opts = [ False, True ]
  else:
    raise ValueError(f"Invalid option: {significance}")
  for entanglement_var in entanglement_vars:
    for is_significance in significance_opts:
      def get_args(init_method = None, decay_mode = None, spin_analyzer = None):
        arg = {
          'data'             : data,
          'entanglement_var' : entanglement_var,
          'is_significance'  : is_significance,
          'central'          : central,
          'output_dir'       : output_dir_abspath,
          'exts'             : exts,
          'show_title'       : show_title,
        }
        if init_method:
          arg['init_method'] = init_method
        if decay_mode:
          arg['decay_mode'] = decay_mode
        if spin_analyzer:
          arg['spin_analyzer'] = spin_analyzer
        return arg
      if 'individual' in by:
        for init_method in INIT_METHODS:
          for decay_mode in DECAY_MODES:
            for spin_analyzer in SPIN_ANALYZERS:
              args.append(get_args(init_method, decay_mode, spin_analyzer))
      if 'init_method' in by:
        for init_method in INIT_METHODS:
          for decay_mode in DECAY_MODES:
            args.append(get_args(init_method, decay_mode, None))
      if 'decay_mode' in by:
        for init_method in INIT_METHODS:
          for spin_analyzer in SPIN_ANALYZERS:
            args.append(get_args(init_method, None, spin_analyzer))
      if 'spin_analyzer' in by:
        for decay_mode in DECAY_MODES:
          for spin_analyzer in SPIN_ANALYZERS:
            args.append(get_args(None, decay_mode, spin_analyzer))

  num_plots = 0
  for arg in args:
    num_plots += plot_by(**arg)
  print(f'Created {num_plots} figure(s) in total')

if __name__ == '__main__':
  parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('-i', '--input', nargs = '*', required = True, help = 'List of directories containing the JSON files')
  parser.add_argument('-s', '--sample-name', type = str, default = 'dy_lo_pythia8_ext', help = 'Sample name')
  parser.add_argument('-c', '--central', type = str, choices = CENTRAL_CHOICES, default = DEFAULT_CENTRAL, help = 'What to use as the central value')
  parser.add_argument('-e', '--entanglement-vars', nargs = '*', choices = ENTANGLEMENT_VARS.keys(), default = [ "Rchsh", "concurrence" ], help = 'Entanglement variables')
  parser.add_argument('-p', '--plot-targets', nargs = '*', choices = PLOT_CHOICES, default = PLOT_CHOICES_DEFAULT, help = 'What to plot')
  parser.add_argument('-S', '--significance', type = str, choices = [ 'only', 'excluded', 'both' ], default = 'only', help = 'Entanglement variable or its significance')
  parser.add_argument('-o', '--output', type = str, required = True, help = 'Output directory')
  parser.add_argument('-a', '--axis', type = str, choices = AXIS_CHOICES, default = DEFAULT_AXIS, help = 'Coordinate system')
  parser.add_argument('-E', '--extensions', nargs = '*', choices = EXTENSIONS, default = EXTENSIONS, help = 'File extensions')
  parser.add_argument('-t', '--show-title', action = 'store_true', default = False, help = 'Show title on every plot')
  args = parser.parse_args()

  # need to define:
  input_dirs = args.input
  sample_name = args.sample_name
  central = args.central
  entanglement_vars = args.entanglement_vars
  plot_targets = args.plot_targets
  significance = args.significance
  output_dir = args.output
  axis = args.axis
  exts = args.extensions
  show_title = args.show_title

  data = read_data(input_dirs)
  assert(sample_name in data)
  assert(axis in data[sample_name])
  sample_data = data[sample_name][axis]
  plot_all(sample_data, entanglement_vars, significance, by = plot_targets, exts = exts, central = central, output_dir = output_dir, show_title = show_title)
