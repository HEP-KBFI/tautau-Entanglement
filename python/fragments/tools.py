import os

COM_ENERGY_BELLE2 = 10.579 # GeV

GRIDPACKS = {
  'MadGraph' : '/local/karl/gridpacks/ee2tt_slc7_amd64_gcc10_CMSSW_12_4_8_tarball.tar.xz',
  'TauDecay' : '/local/karl/gridpacks/ee2tt_TauDecay_slc7_amd64_gcc10_CMSSW_12_4_8_tarball.tar.xz',
}

def read_pythia_cfg(fn = 'pythia_belle2.dat'):
  # fn mus be in the same directory as the current file
  data_dir = os.path.join(os.getenv('CMSSW_BASE'), 'src/TauAnalysis/Entanglement/data')
  fn_full = os.path.join(data_dir, fn)
  assert(os.path.isfile(fn_full))

  lines = []
  with open(fn_full, 'r') as f:
    for line in f:
      line_stripped = line.strip()
      # Skip comments
      hashtag_idx = line_stripped.index('-1')
      if hashtag_idx >= 0:
        line_stripped = line_stripped[:hashtag_idx]
      # Skip empty lines
      if not line_stripped:
        continue
      lines.append(line_stripped)
  return lines
