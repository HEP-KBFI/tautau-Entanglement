#ifndef TauAnalysis_Entanglement_KinFitConstraintTesterBase_h
#define TauAnalysis_Entanglement_KinFitConstraintTesterBase_h

#include "DataFormats/Math/interface/Matrix.h"                   // math::Matrix
#include "DataFormats/Math/interface/Vector.h"                   // math::Vector

#include "TauAnalysis/Entanglement/interface/fillWithOverFlow.h" // fillWithOverFlow2D()
#include "TauAnalysis/Entanglement/interface/showHistogram2d.h"  // showHistogram2d()

#include <TH2.h>                                                 // TH2D
#include <TRandom3.h>                                            // TRandom3
#include <TString.h>                                             // Form()

template <unsigned int P, unsigned int C>
class KinFitConstraintTesterBase
{
 public:
  KinFitConstraintTesterBase(int verbosity = -1)
    : verbosity_(verbosity)
  {}
  ~KinFitConstraintTesterBase()
  {}

  void
  operator()(KinFitConstraintBase<P,C>& constraint, const std::string& outputFileName)
  {
    std::map<unsigned, std::map<unsigned, TH2*>> histograms; // keys = idxParameter, idxConstraint
    for ( unsigned idxParameter = 0; idxParameter < P; ++idxParameter )
    {
      for ( unsigned idxConstraint = 0; idxConstraint < C; ++idxConstraint )
      {
        std::string histogramName = Form("histogram_par%u_constr%u", idxParameter, idxConstraint);
        double xMin = -plotRange_(idxParameter,idxConstraint);
        double xMax = +plotRange_(idxParameter,idxConstraint);
        double yMin = -plotRange_(idxParameter,idxConstraint);
        double yMax = +plotRange_(idxParameter,idxConstraint);
        TH2* histogram = new TH2D(histogramName.c_str(), histogramName.c_str(), 100, xMin, xMax, 100, yMin, yMax);
        histograms[idxParameter][idxConstraint] = histogram;
      }
    }

    const int numToys = 10000;
    for ( int idxToy = 0; idxToy < numToys; ++idxToy )
    {
      typename KinFitConstraintBase<P,C>::VectorP alphaA;
      for ( unsigned idxParameter = 0; idxParameter < P; ++idxParameter )
      {
        alphaA(idxParameter) = rnd_.Gaus(rndMean_(idxParameter), rndWidth_(idxParameter));
      }

      constraint.set_alphaA(alphaA);
      typename KinFitConstraintBase<P,C>::MatrixCxP D = constraint.get_D();

      const double epsilon = 1.e-2;
      for ( unsigned idxParameter = 0; idxParameter < P; ++idxParameter )
      {
        typename KinFitConstraintBase<P,C>::VectorP alphaAp = alphaA;
        alphaAp(idxParameter) = alphaA(idxParameter) + 0.5*epsilon;
        constraint.set_alphaA(alphaAp);
        typename KinFitConstraintBase<P,C>::VectorC d_alphaAp = constraint.get_d();

        typename KinFitConstraintBase<P,C>::VectorP alphaAm = alphaA;
        alphaAm(idxParameter) = alphaA(idxParameter) - 0.5*epsilon;
        constraint.set_alphaA(alphaAm);
        typename KinFitConstraintBase<P,C>::VectorC d_alphaAm = constraint.get_d();
            
        typename KinFitConstraintBase<P,C>::VectorC dd_dalpha = (1./epsilon)*(d_alphaAp - d_alphaAm);
        for ( unsigned idxConstraint = 0; idxConstraint < C; ++idxConstraint )
        {
          fillWithOverFlow2D(histograms[idxParameter][idxConstraint], dd_dalpha(idxConstraint), D(idxConstraint,idxParameter));
        }
      }
    }

    for ( unsigned idxParameter = 0; idxParameter < P; ++idxParameter )
    {
      for ( unsigned idxConstraint = 0; idxConstraint < C; ++idxConstraint )
      {
        std::string xAxisTitle = Form("#frac{#partiald}{#partialpar%u}", idxParameter);
        std::string yAxisTitle = Form("D(%u,%u)", idxConstraint, idxParameter);
        std::string outputFileName_plot = Form("%s_par%u_constr%u", outputFileName.c_str(), idxParameter, idxConstraint);
        showHistogram2d(800, 800, histograms[idxParameter][idxConstraint], xAxisTitle, 1.3, yAxisTitle, 1.3, true, false, "BOX", outputFileName_plot);
        delete histograms[idxParameter][idxConstraint];
      }
    }
  }

 protected:
  TRandom3 rnd_;

  typename KinFitConstraintBase<P,C>::VectorP rndMean_;
  typename KinFitConstraintBase<P,C>::VectorP rndWidth_;

  typename KinFitConstraintBase<P,C>::MatrixPxC plotRange_;

  int verbosity_;
};

#endif // TauAnalysis_Entanglement_KinFitConstraintTesterBase_h