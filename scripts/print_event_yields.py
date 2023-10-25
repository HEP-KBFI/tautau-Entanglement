#!/usr/bin/env python3

import argparse
import os
import re
import json
import prettytable

DMs = [ "had_had", "pi_pi", "pi_rho", "pi_a1", "rho_rho", "rho_a1", "a1_a1" ]
INIT_MODES = [ 'gen', 'startPos', 'kinFit' ]
AXIS_CHOICES = [ 'beam', 'higgs' ]
DEFAULT_AXIS = 'beam'

def fmt_cut(cut):
  return 'Inclusive' if cut == 'incl' else f"|cosTheta| < {cut.replace('p', '.')}"

if __name__ == '__main__':
  parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
  parser.add_argument('-i', '--input', nargs = '*', required = True, help = 'List of directories containing the JSON files')
  parser.add_argument('-s', '--sample-name', type = str, default = 'dy_lo_pythia8_ext', help = 'Sample name')
  parser.add_argument('-m', '--mode', type = str, required = True, choices = INIT_MODES, help = 'Init mode')
  parser.add_argument('-a', '--axis', type = str, choices = AXIS_CHOICES, default = DEFAULT_AXIS, help = 'Coordinate system')
  args = parser.parse_args()

  input_dirs = args.input
  sample_name = args.sample_name
  init_mode = args.mode
  axis = args.axis

  rgx_incl = re.compile(
    rf'analyzeEntanglementNtuple_{sample_name}_{init_mode}Mode_{axis}Axis_(?P<dm>.+)DecayMode_by_(?P<spin_analyzer>\w+).json'
  )
  rgx_bin = re.compile(
    rf'analyzeEntanglementNtuple_{sample_name}_{init_mode}Mode_{axis}Axis_(?P<dm>.+)DecayMode_by_(?P<spin_analyzer>\w+)_absCosTheta(?P<cut>.+).json'
  )

  data = {}
  dms = []
  for input_dir in input_dirs:
    if not os.path.isdir(input_dir):
      raise ValueError(f'No such directory: {input_dir}')

    for fn in os.listdir(input_dir):
      match_incl = rgx_incl.match(fn)
      match_bin = rgx_bin.match(fn)
      fp = os.path.join(input_dir, fn)

      match = None
      if match_bin:
        match = match_bin
      elif match_incl:
        match = match_incl
      else:
        continue
      assert(match)

      cut = match.group('cut') if 'cut' in match.groupdict() else 'incl'
      if cut not in data:
        data[cut] = {}

      with open(fp, 'r') as f:
        fd = json.load(f)
        ev_yield = int(fd['metadata']['sample_size'])

      dm = match.group('dm')
      if dm not in dms:
        dms.append(dm)

      if dm not in data[cut]:
        data[cut][dm] = ev_yield
      else:
        assert(ev_yield == data[cut][dm])

  assert(all(dm in DMs for dm in dms))
  cuts_sorted = list(sorted(data.keys(), key = lambda cut: float(cut.replace('p', '.')) if 'p' in cut else 1, reverse = True))
  dms_sorted = [ dm for dm in DMs if dm in dms ]
  cut_yields = { cut : data[cut]['had_had'] for cut in cuts_sorted if 'had_had' in data[cut] }
  dm_yields = { dm : data['incl'][dm] for dm in dms_sorted }

  table = prettytable.PrettyTable([ r'Cut \/ DM ->' ] + dms_sorted)
  for cut in cuts_sorted:
    row = [ fmt_cut(cut) ]
    for dm in dms_sorted:
      if dm not in data[cut]:
        row.append('N/A')
      else:
        row.append(f'{data[cut][dm]:.0f} ({data[cut][dm] / cut_yields[cut] * 100.:.2f}%, {data[cut][dm] / dm_yields[dm] * 100.:.2f}%)')
    table.add_row(row)
  print(table)
