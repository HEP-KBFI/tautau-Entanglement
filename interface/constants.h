#ifndef TauAnalysis_Entanglement_constants_h
#define TauAnalysis_Entanglement_constants_h

// define charged and neutral pion, proton and tau lepton mass.
// The values are taken from Prog. Theor. Exp. Phys. 2022 (2022) 083C01
const double mChargedPion = 0.139571;
const double mNeutralPion = 0.134977;
const double mProton      = 0.938272;
const double mTau         = 1.77686;

// define electroweak coupling constant
// cf. Eq. (15) in Comp. Phys. Commun. 64 (1991) 275
const double gamma_va = 1.;

// define B field (assumed to be in z direction)
const double Bfield = 4.; // [T]

// define reference point for track parametrization to coincide with nominal interaction point (0,0,0)
const double xr     = 0.; // [cm]
const double yr     = 0.; // [cm]
const double zr     = 0.; // [cm]

#endif // TauAnalysis_Entanglement_constants_h
