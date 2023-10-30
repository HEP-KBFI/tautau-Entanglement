#!/usr/bin/env python3

from TauAnalysis.Entanglement.tools.aux import ENTANGLEMENT_VARS, CENTRAL_CHOICES, AXIS_CHOICES, SPIN_ANALYZERS, DECAY_MODES, INIT_METHODS, EXTENSIONS
from TauAnalysis.Entanglement.tools.plot_settings import CMS
from TauAnalysis.Entanglement.tools.jobTools import mkdir

import argparse
import ROOT
import matplotlib.pyplot as plt
import os
import re

import matplotlib as mpl
for k, v in CMS.items():
  if k in mpl.rcParams:
    mpl.rcParams[k] = v

XPLOT_SETTINGS = {
  'cosThetaStar' : { 'symbol' : r'$\cos\theta^*$', 'xmin' : -1, 'xmax' : +1 },
  'zPlus'        : { 'symbol' : r'$z^+$',          'xmin' :  0, 'xmax' : +1 },
  'zMinus'       : { 'symbol' : r'$z^-$',          'xmin' :  0, 'xmax' : +1 },
}
YPLOT_SETTINGS = {
  'Rchsh' : { 'symbol' : r'$\mathfrak{m}_{12}$', 'ref' : +1 },
  'C_nn'  : { 'symbol' : r'$C_{nn}$',            'ref' : None },
  'C_rr'  : { 'symbol' : r'$C_{rr}$',            'ref' : None },
  'C_kk'  : { 'symbol' : r'$C_{kk}$',            'ref' : None },
}

FN_RGX = re.compile(
  r"analyzeEntanglementNtuple_(?P<sample_name>\w+)_(?P<init_mode>\w+)Mode_(?P<axis>\w+)Axis_(?P<dm>\w+)DecayMode_by_(?P<spin_analyzer>\w+).root"
)

MARKERS = [ 'o', 'v', '^', 's', 'X', 'd', '*' ]

def read_1d_histogram(fn, hname):
  assert(os.path.isfile(fn))
  f = ROOT.TFile.Open(fn, 'read')
  h = f.Get(hname)
  assert(h)
  assert('TH1' in str(type(h)))
  xaxis = h.GetXaxis()
  nbins = xaxis.GetNbins()
  xcoords = [ xaxis.GetBinCenter(i) for i in range(1, nbins + 1) ]
  ycoords = [ h.GetBinContent(i) for i in range(1, nbins + 1) ]
  yerr =  [ h.GetBinError(i) for i in range(1, nbins + 1) ]
  f.Close()
  return { 'x' : xcoords, 'y' : ycoords, 'yerr' : yerr }

def plot(data, outfn, title, data_labels = None, xlabel = '', ylabel = '', xmin = None, xmax = None, yref = None, exts = None, verbose = True):
  if exts is None:
    exts = [ 'pdf', 'png' ]

  if data_labels:
    assert(len(data) == len(data_labels))

  plt.figure(figsize = (10, 8), dpi = 100)

  for idx, entry in enumerate(data):
    p = plt.scatter(entry['x'], entry['y'], marker = MARKERS[idx], s = 40, label = data_labels[idx] if data_labels else '')
    plt.errorbar(entry['x'], entry['y'], yerr = entry['yerr'], ls = 'none', capsize = 5, c = p.get_facecolor())

  if yref is not None:
    plt.axhline(y = yref, ls = '--', c = 'black')
  plt.grid(True)
  if xmin is not None:
    plt.xlim(xmin = xmin)
  if xmax is not None:
    plt.xlim(xmax = xmax)
  if xlabel:
    plt.xlabel(xlabel)
  if ylabel:
    plt.ylabel(ylabel)
  plt.title(title, fontsize = 14)
  if data_labels:
    legend = plt.legend(bbox_to_anchor = (1.3, 0.6))
    plt.setp(legend.get_title(), fontsize = 'xx-small')
  for ext in exts:
    fn = f'{outfn}.{ext}'
    plt.savefig(fn, bbox_inches = 'tight')
    if verbose:
      print(f'Saved file {fn}')
  plt.close()

NOT_AVAILABLE = 'N/A'

if __name__ == "__main__":

  parser = argparse.ArgumentParser()
  group = parser.add_mutually_exclusive_group()

  parser.add_argument('-i', '--input-dirs', type = str, nargs = '+', default = [], help = 'Directories where ROOT files reside')
  parser.add_argument('-o', '--output-dir', type = str, required = True, help = 'Output directory where plots are saved')
  parser.add_argument('-a', '--axis', type = str, choices = AXIS_CHOICES, default = 'beam', help = 'Coordinate system')
  parser.add_argument('-x', '--xvar', type = str, choices = XPLOT_SETTINGS.keys(), default = 'cosThetaStar', help = 'Variable on the x axis')
  parser.add_argument('-y', '--yvar', type = str, choices = YPLOT_SETTINGS.keys(), default = 'Rchsh', help = 'Variable on the y axis')
  parser.add_argument('-e', '--extensions', type = str, nargs = '+', choices = EXTENSIONS, default = [ 'png' ], help = 'File extension of output files')
  group.add_argument("-I", "--init-modes", type = str, choices = INIT_METHODS.keys(), nargs = '*', default = NOT_AVAILABLE, help = "Initialization modes")
  group.add_argument("-D", "--decay-modes", type = str, choices = DECAY_MODES.keys(), nargs = '*', default = NOT_AVAILABLE, help = "Decay modes")
  group.add_argument("-S", "--spin-analyzers", type = str, choices = SPIN_ANALYZERS.keys(), nargs = '*', default = NOT_AVAILABLE, help = "Spin analyzer mode")
  group.add_argument("-s", "--samples", type = str, nargs = '+', help = "MC samples")
  args = parser.parse_args()

  input_dirs = args.input_dirs
  output_dir = os.path.abspath(args.output_dir)
  axis = args.axis
  init_modes = args.init_modes
  decay_modes = args.decay_modes
  spin_analyzers = args.spin_analyzers
  samples = args.samples
  xvar = args.xvar
  yvar = args.yvar
  extensions = args.extensions

  if init_modes == []:
    init_modes = INIT_METHODS.keys()
  elif init_modes == NOT_AVAILABLE:
    init_modes = []
  if decay_modes == []:
    decay_modes = DECAY_MODES.keys()
  elif decay_modes == NOT_AVAILABLE:
    decay_modes = []
  if spin_analyzers == []:
    spin_analyzers = SPIN_ANALYZERS.keys()
  elif spin_analyzers == NOT_AVAILABLE:
    spin_analyzers = []

  hname = f'{xvar}_{yvar}'
  mkdir(output_dir, True)
  xlabel = XPLOT_SETTINGS[xvar]['symbol']
  ylabel = YPLOT_SETTINGS[yvar]['symbol']
  xmin, xmax = XPLOT_SETTINGS[xvar]['xmin'], XPLOT_SETTINGS[xvar]['xmax']
  yref = YPLOT_SETTINGS[yvar]['ref']

  data = {}
  available_sample_names = []
  available_init_modes = []
  available_decay_modes = []
  available_spin_analyzers = []

  for input_dir in input_dirs:
    for input_filename in os.listdir(input_dir):
      match = FN_RGX.match(input_filename)
      if not match:
        continue
      if match.group('axis') != axis:
        continue
      # sample name
      sample_name = match.group('sample_name')
      if samples and sample_name not in samples:
        continue
      if sample_name not in data:
        data[sample_name] = {}
      if sample_name not in available_sample_names:
        available_sample_names.append(sample_name)
      # init mode
      init_mode = match.group('init_mode')
      if init_modes and init_mode not in init_modes:
        continue
      if init_mode not in data[sample_name]:
        data[sample_name][init_mode] = {}
      if init_mode not in available_init_modes:
        available_init_modes.append(init_mode)
      # decay mode
      dm = match.group('dm')
      if decay_modes and dm not in decay_modes:
        continue
      if dm not in data[sample_name][init_mode]:
        data[sample_name][init_mode][dm] = {}
      if dm not in available_decay_modes:
        available_decay_modes.append(dm)
      # spin analyzer
      spin_analyzer = match.group('spin_analyzer')
      if spin_analyzers and spin_analyzer not in spin_analyzers:
        continue
      data[sample_name][init_mode][dm][spin_analyzer] = read_1d_histogram(os.path.join(input_dir, input_filename), hname)
      if spin_analyzer not in available_spin_analyzers:
        available_spin_analyzers.append(spin_analyzer)

  common_args = {
    'xlabel' : xlabel,
    'ylabel' : ylabel,
    'xmin'   : xmin,
    'xmax'   : xmax,
    'yref'   : yref,
    'exts'   : extensions,
  }
  if init_modes:
    for sample_name in available_sample_names:
      for dm in available_decay_modes:
        for spin_analyzer in available_spin_analyzers:
          plot_data = []
          plot_labels = []
          for init_mode in init_modes:
            if sample_name in data and \
               init_mode in data[sample_name] and \
               dm in data[sample_name][init_mode] and \
               spin_analyzer in data[sample_name][init_mode][dm]:
              plot_data.append(data[sample_name][init_mode][dm][spin_analyzer])
              plot_labels.append(INIT_METHODS[init_mode])
          if not plot_data:
            continue
          outfn = os.path.join(output_dir, f'plot_{hname}_{sample_name}_{dm}DecayMode_{spin_analyzer}SpinAnalyzer_byMode')
          plot(
            data = plot_data,
            outfn = outfn,
            title = f'Sample {sample_name}, decay mode = {dm}, analyzer = {spin_analyzer}',
            data_labels = plot_labels,
            **common_args
          )
  elif decay_modes:
    for sample_name in available_sample_names:
      for init_mode in available_init_modes:
        for spin_analyzer in available_spin_analyzers:
          plot_data = []
          plot_labels = []
          for dm in decay_modes:
            if sample_name in data and \
               init_mode in data[sample_name] and \
               dm in data[sample_name][init_mode] and \
               spin_analyzer in data[sample_name][init_mode][dm]:
              plot_data.append(data[sample_name][init_mode][dm][spin_analyzer])
              plot_labels.append(DECAY_MODES[dm]['label'])
          if not plot_data:
            continue
          outfn = os.path.join(output_dir, f'plot_{hname}_{sample_name}_{init_mode}Mode_{spin_analyzer}SpinAnalyzer_byDecayMode')
          plot(
            data = plot_data,
            outfn = outfn,
            title = f'Sample {sample_name}, mode = {init_mode}, analyzer = {spin_analyzer}',
            data_labels = plot_labels,
            **common_args
          )
  elif spin_analyzers:
    for sample_name in available_sample_names:
      for init_mode in available_init_modes:
        for dm in available_decay_modes:
          plot_data = []
          plot_labels = []
          for spin_analyzer in spin_analyzers:
            if sample_name in data and \
               init_mode in data[sample_name] and \
               dm in data[sample_name][init_mode] and \
               spin_analyzer in data[sample_name][init_mode][dm]:
              plot_data.append(data[sample_name][init_mode][dm][spin_analyzer])
              plot_labels.append(SPIN_ANALYZERS[spin_analyzer])
          if not plot_data:
            continue
          outfn = os.path.join(output_dir, f'plot_{hname}_{sample_name}_{init_mode}Mode_{dm}DecayMode_bySpinAnalyzer')
          plot(
            data = plot_data,
            outfn = outfn,
            title = f'Sample {sample_name}, mode = {init_mode}, decay mode = {dm}',
            data_labels = plot_labels,
            **common_args
          )
  elif samples:
    for init_mode in available_init_modes:
      for dm in available_decay_modes:
        for spin_analyzer in available_spin_analyzers:
          plot_data = []
          plot_labels = []
          for sample_name in samples:
            if init_mode in data[sample_name] and \
               dm in data[sample_name][init_mode] and \
               spin_analyzer in data[sample_name][init_mode][dm]:
              plot_data.append(data[sample_name][init_mode][dm][spin_analyzer])
              plot_labels.append(sample_name)
          if not plot_data:
            continue
          outfn = os.path.join(output_dir, f'plot_{hname}_{init_mode}Mode_{dm}DecayMode_{spin_analyzer}SpinAnalyzer_bySample')
          plot(
            data = plot_data,
            outfn = outfn,
            title = f'Mode = {init_mode}, decay mode = {dm}, analyzer = {spin_analyzer}',
            data_labels = plot_labels,
            **common_args
          )
  else:
    for sample_name in data:
      for init_mode in data[sample_name]:
        for dm in data[sample_name][init_mode]:
          for spin_analyzer in data[sample_name][init_mode][dm]:
            outfn = os.path.join(output_dir, f'plot_{hname}_{sample_name}_{init_mode}Mode_{dm}DecayMode_{spin_analyzer}SpinAnalyzer')
            plot(
              data = [ data[sample_name][init_mode][dm][spin_analyzer] ],
              outfn = outfn,
              title = f'Sample {sample_name}, mode = {init_mode}, decay mode = {dm}, analyzer = {spin_analyzer}',
              data_labels = None,
              **common_args
            )
