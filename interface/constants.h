#ifndef TauAnalysis_Entanglement_constants_h
#define TauAnalysis_Entanglement_constants_h

// define charged and neutral pion, proton and tau lepton mass.
// The values are taken from Prog. Theor. Exp. Phys. 2022 (2022) 083C01
// and are referred to as "PDG values"
const double mElectron                   =   0.511e-3; // 511 keV
const double mChargedPion                =   0.139571; // [GeV]
// CV: use Higgs boson mass used in Monte Carlo production instead of "PDG value"
//const double mHiggs                      = 125.25;      // [GeV]
const double mHiggs                      = 125.0;      // [GeV]
const double mNeutralPion                =   0.134977; // [GeV]
const double mProton                     =   0.938272; // [GeV]
const double mTau                        =   1.77686;  // [GeV]

// define tau lepton lifetime 
// The value is defined as the expected tau decay distance in the restframe of the tau lepton
// and is taken from Prog. Theor. Exp. Phys. 2022 (2022) 083C01
const double ct                          =   8.7e-3;   // 87 micrometer

// define electroweak coupling constant
// cf. Eq. (15) in Comp. Phys. Commun. 64 (1991) 275
const double gamma_va                    =   1.;

// define B field (assumed to be in z direction)
const double Bfield                      =   4.;       // [T]

// define reference point for track parametrization to coincide with nominal interaction point (0,0,0)
const double xr                          =   0.;       // [cm]
const double yr                          =   0.;       // [cm]
const double zr                          =   0.;       // [cm]

// define beam energies for LHC and SuperKEKB (Belle)
enum { kLHC, kSuperKEKB };
const double beamEnergy_LHC              = 7.e+3;      // [GeV]
const double beamEnergy_SuperKEKB_ePlus  = 3.844;      // [GeV]
const double beamEnergy_SuperKEKB_eMinus = 6.734;      // [GeV]

#endif // TauAnalysis_Entanglement_constants_h
