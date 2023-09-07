#ifndef TauAnalysis_Entanglement_fillWithOverFlow_h
#define TauAnalysis_Entanglement_fillWithOverFlow_h

// forward declarations
class TH1;
class TH2;

void
fill1D(TH1 * histogram,
       double x,
       double evtWeight = 1.);

void
fillWithOverFlow1D(TH1* histogram,
                   double x,
                   double evtWeight = 1.);

void
fill2D(TH2 * histogram,
       double x,
       double y,
       double evtWeight = 1.);

void
fillWithOverFlow2D(TH2* histogram,
                   double x,
                   double y,
                   double evtWeight = 1.);

#endif // TauAnalysis_Entanglement_fillWithOverFlow_h
