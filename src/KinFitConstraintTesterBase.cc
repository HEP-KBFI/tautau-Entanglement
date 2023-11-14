#include "TauAnalysis/Entanglement/interface/KinFitConstraintTesterBase.h"

#include "TauAnalysis/Entanglement/interface/cmsException.h"     // cmsException
#include "TauAnalysis/Entanglement/interface/fillWithOverFlow.h" // fillWithOverFlow2D()
#include "TauAnalysis/Entanglement/interface/showHistogram2d.h"  // showHistogram2d()

#include <TH2.h>                                                 // TH2D
#include <TString.h>                                             // Form()

#include <iostream>                                              // std::cout
#include <map>                                                   // std::map<>

KinFitConstraintTesterBase::KinFitConstraintTesterBase(const TVectorD& alpha0, int verbosity)
  : alpha0_(alpha0)
  , verbosity_(verbosity)
{}

KinFitConstraintTesterBase::~KinFitConstraintTesterBase()
{}

void
KinFitConstraintTesterBase::operator()(KinFitConstraintBase& constraint, const std::string& outputFileName)
{
  if ( verbosity_ >= 1 )
  {
    std::cout << "<KinFitConstraintTesterBase::operator()>:\n";
  }

  const unsigned int Np = constraint.get_Np();
  const unsigned int Nc_eq = constraint.get_Nc_eq();
  if ( alpha0_.GetNrows() != (int)Np )
    throw cmsException("KinFitConstraintTesterBase::operator()", __LINE__) 
      << "Dimension of parameter vector alpha0 does not match constraint !!\n";

  std::map<unsigned int, std::map<unsigned, TH2*>> histograms; // keys = idxConstraint, idxParameter
  for ( unsigned int idxConstraint = 0; idxConstraint < Nc_eq; ++idxConstraint )
  {
    for ( unsigned int idxParameter = 0; idxParameter < Np; ++idxParameter )
    {
      std::string histogramName = Form("histogram_constr%u_par%u", idxConstraint, idxParameter);
      double xMin = -plotRange_(idxConstraint,idxParameter);
      double xMax = +plotRange_(idxConstraint,idxParameter);
      double yMin = -plotRange_(idxConstraint,idxParameter);
      double yMax = +plotRange_(idxConstraint,idxParameter);
      TH2* histogram = new TH2D(histogramName.c_str(), histogramName.c_str(), 100, xMin, xMax, 100, yMin, yMax);
      histograms[idxConstraint][idxParameter] = histogram;
    }
  }

  const unsigned int numToys = 10000;
  unsigned int numToys_selected = 0;
  for ( unsigned int idxToy = 0; idxToy < numToys; ++idxToy )
  {
    if ( ((idxToy % 1000) == 0 && verbosity_ >= 1) || verbosity_ >= 3 )
    {
      std::cout << " processing toy " << idxToy << "/" << numToys << " (selected = " << numToys_selected << ")\n";
    }

    TVectorD alphaA(Np);
    for ( unsigned int idxParameter = 0; idxParameter < Np; ++idxParameter )
    {
      alphaA(idxParameter) = rnd_.Gaus(rndMean_(idxParameter), rndWidth_(idxParameter));
    }

    constraint.set_alphaA(alphaA);
    if ( constraint.get_errorFlag() )
    {
      if ( verbosity_ >= 1 )
      {
	std::cout << "skipping toy, because constraint returned an error !!\n";
      }
      continue;
    }	  
      
    const TMatrixD& D = constraint.get_D_eq();

    const double epsilon = 1.e-2;
    for ( unsigned int idxParameter = 0; idxParameter < Np; ++idxParameter )
    {
      TVectorD alphaAp = alphaA;
      alphaAp(idxParameter) = alphaA(idxParameter) + 0.5*epsilon;
      constraint.set_alphaA(alphaAp);
      TVectorD d_alphaAp = constraint.get_d_eq();

      TVectorD alphaAm = alphaA;
      alphaAm(idxParameter) = alphaA(idxParameter) - 0.5*epsilon;
      constraint.set_alphaA(alphaAm);
      TVectorD d_alphaAm = constraint.get_d_eq();
            
      TVectorD dd_dalpha = (1./epsilon)*(d_alphaAp - d_alphaAm);
      for ( unsigned int idxConstraint = 0; idxConstraint < Nc_eq; ++idxConstraint )
      {
        fillWithOverFlow2D(histograms[idxConstraint][idxParameter], dd_dalpha(idxConstraint), D(idxConstraint,idxParameter));
      }
    }

    ++numToys_selected;
  }

  if ( verbosity_ >= 1 )
  {
    std::cout << "#toys: analyzed = " << numToys << ", selected = " << numToys_selected << "\n";
  }

  for ( unsigned int idxConstraint = 0; idxConstraint < Nc_eq; ++idxConstraint )
  {
    for ( unsigned int idxParameter = 0; idxParameter < Np; ++idxParameter )
    {
      std::string xAxisTitle = Form("#frac{#partiald}{#partialpar%u}", idxParameter);
      std::string yAxisTitle = Form("D(%u,%u)", idxConstraint, idxParameter);
      size_t idx = outputFileName.find_last_of('.');
      std::string outputFileName_plot = std::string(outputFileName, 0, idx);
      outputFileName_plot.append("_");
      outputFileName_plot.append(Form("constr%u_par%u", idxConstraint, idxParameter));
      outputFileName_plot.append(std::string(outputFileName, idx));
      showHistogram2d(800, 800, histograms[idxConstraint][idxParameter], xAxisTitle, 1.3, yAxisTitle, 1.3, true, false, "BOX", outputFileName_plot, false, verbosity_);
      delete histograms[idxConstraint][idxParameter];
    }
  }
}
