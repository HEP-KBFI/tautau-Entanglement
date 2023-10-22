#!/usr/bin/env python3

# Example usage:
#
# get_latex_tables.py -i ~/Entanglement/analysis/SuperKEKB/2023Oct17
#
# In case the JSON files are stored in multiple directories, you can list those directories as argument to
# -i/--input as long as there are no files with the same name in those directories.
# If you want to produce Latex tables for a particular cut in |cos(theta)|, you can do that by adding
#
# -c <cut>
#
# to the command, where <cut> can be "opt" if you want tables for the optimal cut, or something like
# "0p40" if you want tables for the cut |cos(theta)| <= 0.4.
#
# To colorize the terminal output, pipe it to: pygmentize -l latex

import argparse
import jinja2
import json
import os
import collections

CENTRAL_CHOICES = [ 'nominal', 'median' ]
ENTANGLEMENT_VARS = [ 'Rchsh', 'concurrence', 'steerability' ]
MODE_CHOICES = [ 'gen', 'startPos', 'kinFit' ]
DM_CHOICES = collections.OrderedDict([
  ('pi_pi',   r'$\Pgpp\Pgpm$'),
  ('pi_rho',  r'$\Pgppm\Pgrmp$'),
  ('pi_a1',   r'$\Pgppm\Pamp$'),
  ('rho_rho', r'$\Pgrp\Pgrm$'),
  ('rho_a1',  r'$\Pgrpm\Pamp$'),
  ('a1_a1',   r'$\Pap\Pam$'),
  ('had_had', 'All channels'),
])
SPIN_ANALYZER_CHOICES = collections.OrderedDict([
  ('summation',          'Exp. value'),
  ('differentialXsec2d', '2d distr.'),
  ('differentialXsec1d', '1d distr.'),
  ('asymmetry',          'FB asymm.'),
  ('mlfit',              'ML fit'),
])
AXIS_CHOICES = [ 'beam', 'higgs' ]
COORDINATE_CHOICES = [ 'n', 'r', 'k' ]

DEFAULT_CENTRAL = 'nominal'
DEFAULT_AXIS = 'beam'
DEFAULT_SPIN_ANALYZER = 'mlfit'

TABLE2_TEMPLATE = r"""
<b>if comment</b>% <v>comment</v><b>endif</b>
\begin{tabular}{c|r@{$ \,\,\pm\,\, $}rr@{$ \,\,\pm\,\, $}rr@{$ \,\,\pm\,\, $}r}
Element <b>for i in range(3)</b> & \multicolumn{2}{c}{<v>A[i]</v>}<b>endfor</b> \\
\hline
<b>for i in range(3)</b><v>A[i]</v><b>for j in range(3)</b> & $<v>C[A[i]+A[j]]|fmt(1)</v>$ & $<v>E[A[i]+A[j]]|fmt(1)</v>$<b>endfor</b>\\
<b>endfor</b>\end{tabular}
"""

TABLE3_TEMPLATE = r"""
<b>if comment</b>% <v>comment</v><b>endif</b>
\begin{tabular}{l|r@{$ \,\,\pm\,\, $}rr@{$ \,\,\pm\,\, $}rr@{$ \,\,\pm\,\, $}rr@{$ \,\,\pm\,\, $}r}
Method<b>for e in E</b> & \multicolumn{2}{c}{$\<v>e</v>$}<b>endfor</b> \\
\hline
<b>for s in S</b><v>S[s]|ljust(10)</v><b>for e in E</b> & $<v>D[s][e][C]|fmt</v>$ & $<v>D[s][e]["error"]|fmt</v>$<b>endfor</b> \\
<b>endfor</b>\end{tabular}
"""

TABLE5_TEMPLATE = r"""
<b>if comment</b>% <v>comment</v><b>endif</b>
\begin{tabular}{c|r@{$ \,\,\pm\,\, $}rr@{$ \,\,\pm\,\, $}rr@{$ \,\,\pm\,\, $}rr@{$ \,\,\pm\,\, $}r}
Decay channel<b>for e in E</b> & \multicolumn{2}{c}{$\<v>e</v>$}<b>endfor</b> \\
\hline
<b>for m in M</b><v>M[m]|ljust(14)</v><b>for e in E</b> & $<v>D[m][S][e][C]|fmt</v>$ & $<v>D[m][S][e]["error"]|fmt</v>$<b>endfor</b> \\
<b>endfor</b>\end{tabular}
"""

jinja2_env = jinja2.Environment(
  block_start_string = '<b>',
  block_end_string = '</b>',
  variable_start_string = '<v>',
  variable_end_string = '</v>',
)

def fmt(s, repl_plus = False):
  return f"{s:+.2f}".replace('+', ' ') if repl_plus else f"{s:.2f}"

jinja2_env.filters['fmt'] = fmt
jinja2_env.filters['ljust'] = lambda s, pad: s.ljust(pad)

def read_data(input_dirs, sample_name, modes = None, dms = None, spin_analyzers = None, axis = DEFAULT_AXIS, cut = ''):
  assert(sample_name)
  assert(all(os.path.isdir(input_dir) for input_dir in input_dirs))
  assert(axis in AXIS_CHOICES)
  if modes is None:
    modes = MODE_CHOICES
  else:
    assert(all(mode in MODE_CHOICES for mode in modes) and len(modes) == len(set(modes)))
  if dms is None:
    dms = DM_CHOICES
  else:
    assert(all(dm in DM_CHOICES for dm in dms) and len(dms) == len(set(dms)))
  if spin_analyzers is None:
    spin_analyzers = SPIN_ANALYZER_CHOICES
  else:
    assert(all(spin_analyzer in SPIN_ANALYZER_CHOICES for spin_analyzer in spin_analyzers) and \
              len(spin_analyzers) == len(set(spin_analyzers)))
  data = {}
  for mode in modes:
    data[mode] = {}
    for dm in dms:
      data[mode][dm] = {}
      for spin_analyzer in spin_analyzers:
        fn = f'analyzeEntanglementNtuple_{sample_name}_{mode}Mode_{axis}Axis_{dm}DecayMode_by_{spin_analyzer}'
        if cut:
          fn += f'_{cut}'
        fn += '.json'
        input_dir_selected = []
        for input_dir in input_dirs:
          fp = os.path.join(input_dir, fn)
          if os.path.isfile(fp):
            input_dir_selected.append(input_dir)
        if not input_dir_selected:
          raise RuntimeError("Could not find {} in none of the following dirs: {}".format(fn, ', '.join(input_dirs)))
        if len(input_dir_selected) > 1:
          raise RuntimeError("Found multiple instances of {} in these dirs: {}".format(fn, ', '.join(input_dir_selected)))
        fp = os.path.join(input_dir_selected[0], fn)
        with open(fp, 'r') as f:
          data[mode][dm][spin_analyzer] = json.load(f)
  return data

def get_Cmatrix(data, spin_analyzer = DEFAULT_SPIN_ANALYZER, central = DEFAULT_CENTRAL, coordinates = None, comment = ''):
  assert(central in CENTRAL_CHOICES)
  assert(spin_analyzer in SPIN_ANALYZER_CHOICES)
  if coordinates is None:
    coordinates = COORDINATE_CHOICES
  else:
    assert(all(coordinate in COORDINATE_CHOICES for coordinate in coordinates) and len(coordinates) == len(set(coordinates)))
  Cmatrix_data = data[spin_analyzer]['C']
  return jinja2_env.from_string(TABLE2_TEMPLATE).render(
    C = Cmatrix_data[central],
    E = Cmatrix_data['error'],
    A = coordinates,
    comment = comment,
  )

def get_entanglementVariables(data, central = DEFAULT_CENTRAL, spin_analyzers = None, comment = ''):
  assert(central in CENTRAL_CHOICES)
  if spin_analyzers is None:
    spin_analyzers = SPIN_ANALYZER_CHOICES
  else:
    assert(all(spin_analyzer in SPIN_ANALYZER_CHOICES for spin_analyzer in spin_analyzers) and \
           len(spin_analyzers) == len(set(spin_analyzers)))
  return jinja2_env.from_string(TABLE3_TEMPLATE).render(
    D = data,
    C = central,
    E = ENTANGLEMENT_VARS,
    S = spin_analyzers,
    comment = comment,
  )

def get_decayModes(data, central = DEFAULT_CENTRAL, entanglement_vars = None, spin_analyzer = None, dms = None, comment = ''):
  assert(central in CENTRAL_CHOICES)
  if entanglement_vars is None:
    entanglement_vars = ENTANGLEMENT_VARS
  else:
    assert(all(entanglement_var in ENTANGLEMENT_VARS for entanglement_var in entanglement_vars) and len(entanglement_vars) == len(set(entanglement_vars)))
  assert(spin_analyzer in SPIN_ANALYZER_CHOICES)
  if dms is None:
    dms = DM_CHOICES
  else:
    assert(all(dm in DM_CHOICES for dm in dms) and len(dms) == len(set(dms)))
  return jinja2_env.from_string(TABLE5_TEMPLATE).render(
    D = data,
    C = central,
    E = entanglement_vars,
    S = spin_analyzer,
    M = dms,
    comment = comment,
  )

if __name__ == '__main__':
  parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('-i', '--input', nargs = '*', required = True, help = 'List of directories containing the JSON files')
  parser.add_argument('-s', '--sample-name', type = str, default = 'dy_lo_pythia8_ext', help = 'Sample name')
  parser.add_argument('-c', '--central', type = str, choices = CENTRAL_CHOICES, default = DEFAULT_CENTRAL, help = 'What to quote as the central value')
  parser.add_argument('-S', '--default-spin-analyzer', type = str, choices = SPIN_ANALYZER_CHOICES, default = DEFAULT_SPIN_ANALYZER, help = 'Default spin analyzer')
  parser.add_argument('-a', '--spin-analyzers', nargs = '*', choices = SPIN_ANALYZER_CHOICES, default = SPIN_ANALYZER_CHOICES, help = 'Spin analyzers')
  parser.add_argument('-d', '--decay-modes', nargs = '*', type = str, choices = DM_CHOICES.keys(), default = DM_CHOICES.keys(), help = 'Decay modes')
  parser.add_argument('-e', '--entanglement-vars', nargs = '*', choices = ENTANGLEMENT_VARS, default = ENTANGLEMENT_VARS, help = 'Entanglement variables')
  parser.add_argument('-A', '--axis', type = str, choices = AXIS_CHOICES, default = DEFAULT_AXIS, help = 'Coordinate system')
  parser.add_argument('-C', '--cut', type = str, default = '', help = 'Print tables after applying a cut')
  args = parser.parse_args()

  input_dirs = args.input
  sample_name = args.sample_name
  central = args.central
  default_spin_analyzer = args.default_spin_analyzer
  spin_analyzers = args.spin_analyzers
  decay_modes = collections.OrderedDict([ (dm, DM_CHOICES[dm]) for dm in args.decay_modes ])
  entanglement_vars = args.entanglement_vars
  axis = args.axis
  cut = args.cut

  data = read_data(input_dirs, sample_name, modes = [ 'gen', 'kinFit' ], dms = decay_modes, spin_analyzers = spin_analyzers, axis = axis)
  if cut:
    data_cut = read_data(input_dirs, sample_name, modes = [ 'gen', 'kinFit' ], dms = decay_modes, spin_analyzers = spin_analyzers, axis = axis, cut = cut)
  print(get_Cmatrix(data['gen']['pi_pi'], default_spin_analyzer, central, comment = 'Table 2'))
  print(get_entanglementVariables(data['gen']['pi_pi'], central, spin_analyzers = spin_analyzers, comment = 'Table 3'))
  if cut:
    print(get_entanglementVariables(data_cut['gen']['pi_pi'], central, spin_analyzers = spin_analyzers, comment = 'Table 4'))
    print(get_decayModes(data_cut['gen'], central, entanglement_vars, default_spin_analyzer, decay_modes, comment = 'Table 5'))
    print(get_decayModes(data_cut['kinFit'], central, entanglement_vars, default_spin_analyzer, decay_modes, comment = 'Table 6'))
