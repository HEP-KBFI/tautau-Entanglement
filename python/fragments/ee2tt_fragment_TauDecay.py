from TauAnalysis.Entanglement.fragments.common import externalLHEProducer, generator
from TauAnalysis.Entanglement.fragments.tools import GRIDPACKS
import FWCore.ParameterSet.Config as cms

externalLHEProducer.args = cms.vstring(GRIDPACKS['TauDecay'])
ProductionFilterSequence = cms.Sequence(generator)
