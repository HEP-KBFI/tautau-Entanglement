#include "TauAnalysis/Entanglement/interface/SpinAlgo_by_differentialXsec1d.h"

#include "TauAnalysis/Entanglement/interface/DifferentialXsec_B.h" // DifferentialXsec_B, fit_DifferentialXsec_B()

#include "RooArgSet.h"                                             // RooArgSet
#include "RooGenericPdf.h"                                         // RooGenericPdf
#include "RooRealVar.h"                                            // RooRealVar
#include "RooDataHist.h"                                           // RooDataHist

#include <TH1.h>                                                   // TH1D
#include <TString.h>                                               // Form()

#include <string>                                                  // std::string

using namespace spin;

SpinAlgo_by_differentialXsec1d::SpinAlgo_by_differentialXsec1d(const edm::ParameterSet& cfg)
  : SpinAlgoBase(cfg)
{}

SpinAlgo_by_differentialXsec1d::~SpinAlgo_by_differentialXsec1d()
{}

namespace
{
  // CV: The class DifferentialXsec1d_C represents the single-differential cross section
  //       dsigma/dpsi
  //     given by Eq. (4.16) in the paper arXiv:1508.05271,
  //     with:
  //       psi = cosTheta_p*cosTheta_m
  //     as defined in Eq. (4.14) of the same paper
  class DifferentialXsec1d_C
  {
    public:
    DifferentialXsec1d_C(const std::string& label, int numBins = 40)
      : label_(label)
      , histogram_(nullptr)
    {
      std::string histogramName = Form("histogram1d_%s", label.c_str());
      histogram_ = new TH1D(histogramName.c_str(), histogramName.c_str(), numBins, -1., +1.);
      histogram_->Sumw2();
    }
    ~DifferentialXsec1d_C()
    {
      delete histogram_;
    }

    void
    fill(double cosTheta_p, double cosTheta_m, double evtWeight)
    {
      double psi = cosTheta_p*cosTheta_m;
      histogram_->Fill(psi, evtWeight);
    }

    const std::string&
    get_label() const
    {
      return label_;
    }

    const TH1*
    get_histogram() const
    {
      return histogram_;
    }

   private:
    std::string label_;
    TH1* histogram_;
  };

  double
  fit_DifferentialXsec1d_C(const DifferentialXsec1d_C& Xsec, int verbosity)
  {
    if ( verbosity >= 2 )
    {
      std::cout << "<fit_DifferentialXsec1d_C>: fitting " << Xsec.get_label() << "\n";
    }

    RooRealVar psi("psi", "psi", -1., +1.);
    
    RooDataHist data("data", "data", psi, Xsec.get_histogram());

    // CV: need to restrict range of fit variable C_ij to the interval [-1,+1],
    //     to avoid that PDF becomes less than zero
    //    (which is "unphysical" behaviour for a differential cross section and causes warnings from RooFit)
    RooRealVar C_ij("C_ij", "C_ij", 0., -1., +1.);

    // CV: added protection against small psi values when computing log(1/|psi|)
    RooGenericPdf pdf("pdf", "pdf", "0.5*(1. - C_ij*psi)*log(1./max(1.e-6, abs(psi)))", RooArgSet(C_ij, psi));
    pdf.fitTo(data, RooFit::PrintLevel(-1));
    if ( verbosity >= 2 )
    {
      std::cout << " " << Xsec.get_label() << " = " << C_ij.getVal() << " +/- " << C_ij.getError() << "\n";
    }

    return C_ij.getVal();
  }
}

spin::Measurement
SpinAlgo_by_differentialXsec1d::operator()(const spin::Dataset& dataset)
{
  if ( verbosity_ >= 3 )
  {
    std::cout << "<SpinAlgo_by_differentialXsec1d::operator()>:\n";
    std::cout << " #entries = " << dataset.size() << "\n";
  }

  DifferentialXsec_B Bp_n("Bp_n");
  DifferentialXsec_B Bp_r("Bp_r");
  DifferentialXsec_B Bp_k("Bp_k");

  DifferentialXsec_B Bm_n("Bm_n");
  DifferentialXsec_B Bm_r("Bm_r");
  DifferentialXsec_B Bm_k("Bm_k");

  DifferentialXsec1d_C C_nn("C_nn");
  DifferentialXsec1d_C C_rn("C_rn");
  DifferentialXsec1d_C C_kn("C_kn");
  DifferentialXsec1d_C C_nr("C_nr");
  DifferentialXsec1d_C C_rr("C_rr");
  DifferentialXsec1d_C C_kr("C_kr");
  DifferentialXsec1d_C C_nk("C_nk");
  DifferentialXsec1d_C C_rk("C_rk");
  DifferentialXsec1d_C C_kk("C_kk");

  size_t numEntries = dataset.size();
  for ( size_t idxEntry = 0; idxEntry < numEntries; ++idxEntry )
  {
    const spin::Data& entry = dataset.at(idxEntry);

    double hPlus_n = entry.get_hPlus_n();
    double hPlus_r = entry.get_hPlus_r();
    double hPlus_k = entry.get_hPlus_k();

    double hMinus_n = entry.get_hMinus_n();
    double hMinus_r = entry.get_hMinus_r();
    double hMinus_k = entry.get_hMinus_k();

    double evtWeight = entry.get_evtWeight();

    if ( verbosity_ >= 3 )
    {
      std::cout << "entry #" << idxEntry << ":\n";
      std::cout << " hPlus: n = " << hPlus_n  << ", r = " << hPlus_r  << ", k = " << hPlus_k  << "\n";
      std::cout << " hMinus: n = " << hMinus_n << ", r = " << hMinus_r << ", k = " << hMinus_k << "\n";
      std::cout << " evtWeight = " << evtWeight << "\n";
    }

    // CV: compute polarization vectors B+ and B- for tau+ and tau- according to Eq. (4.18)
    //     in the paper arXiv:1508.05271
    Bp_n.fill(hPlus_n, evtWeight);
    Bp_r.fill(hPlus_r, evtWeight);
    Bp_k.fill(hPlus_k, evtWeight);

    Bm_n.fill(hMinus_n, evtWeight);
    Bm_r.fill(hMinus_r, evtWeight);
    Bm_k.fill(hMinus_k, evtWeight);

    // CV: compute spin correlation matrix C according to Eq. (4.16)
    //     in the paper arXiv:1508.05271
    C_nn.fill(hPlus_n, hMinus_n, evtWeight);
    C_rn.fill(hPlus_r, hMinus_n, evtWeight);
    C_kn.fill(hPlus_k, hMinus_n, evtWeight);
    C_nr.fill(hPlus_n, hMinus_r, evtWeight);
    C_rr.fill(hPlus_r, hMinus_r, evtWeight);
    C_kr.fill(hPlus_k, hMinus_r, evtWeight);
    C_nk.fill(hPlus_n, hMinus_k, evtWeight);
    C_rk.fill(hPlus_r, hMinus_k, evtWeight);
    C_kk.fill(hPlus_k, hMinus_k, evtWeight);
  }

  math::Vector3 Bp;
  Bp(0) = fit_DifferentialXsec_B(Bp_n, verbosity_);
  Bp(1) = fit_DifferentialXsec_B(Bp_r, verbosity_);
  Bp(2) = fit_DifferentialXsec_B(Bp_k, verbosity_);

  math::Vector3 Bm;
  Bm(0) = fit_DifferentialXsec_B(Bm_n, verbosity_);
  Bm(1) = fit_DifferentialXsec_B(Bm_r, verbosity_);
  Bm(2) = fit_DifferentialXsec_B(Bm_k, verbosity_);

  math::Matrix3x3 C;
  C(0,0) = fit_DifferentialXsec1d_C(C_nn, verbosity_);
  C(0,1) = fit_DifferentialXsec1d_C(C_rn, verbosity_);
  C(0,2) = fit_DifferentialXsec1d_C(C_kn, verbosity_);
  C(1,0) = fit_DifferentialXsec1d_C(C_nr, verbosity_);
  C(1,1) = fit_DifferentialXsec1d_C(C_rr, verbosity_);
  C(1,2) = fit_DifferentialXsec1d_C(C_kr, verbosity_);
  C(2,0) = fit_DifferentialXsec1d_C(C_nk, verbosity_);
  C(2,1) = fit_DifferentialXsec1d_C(C_rk, verbosity_);
  C(2,2) = fit_DifferentialXsec1d_C(C_kk, verbosity_);

  spin::Measurement measurement(Bp, Bm, C);
  addEntanglementVariables(measurement);
  return measurement;
}
