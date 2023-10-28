from TauAnalysis.Entanglement.fragments.common import externalLHEProducer, generator
from TauAnalysis.Entanglement.fragments.tools import GRIDPACKS
from GeneratorInterface.ExternalDecays.TauolaSettings_cff import TauolaPolar, TauolaDefaultInputCards
import FWCore.ParameterSet.Config as cms

externalLHEProducer.args = cms.vstring(GRIDPACKS['MadGraph'])
generator.ExternalDecays = cms.PSet( # Automatically delegates tau decays to Tauola
  Tauola = cms.untracked.PSet(
    TauolaPolar,             # Already the default; calls spin_correlation.setAll(true), i.e., enables spin correlations for all taus regardless of their origin
    TauolaDefaultInputCards, # Enables all tau decay modes
    # For specific tau decay modes, uncomment the following and edit as necessary
    #InputCards = cms.PSet
    #(
    #  pjak1 = cms.int32(0), # For tau-; set to 3 for tau->pi+nu, to 4 for tau->rho+nu, to 5 for tau->a1+nu
    #  pjak2 = cms.int32(0), # For tau+; set to 3 for tau->pi+nu, to 4 for tau->rho+nu, to 5 for tau->a1+nu
    #  mdtau = cms.int32(0) # If this is 0 or 1, then DMs read from pjak1, pjak2
    #)
    # DMs documented in https://github.com/PMunkes/evtgen/blob/master/Tauola_README.txt
  ),
  parameterSets = cms.vstring('Tauola')
)
ProductionFilterSequence = cms.Sequence(generator)
