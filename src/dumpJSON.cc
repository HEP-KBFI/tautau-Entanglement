#include "TauAnalysis/Entanglement/interface/dumpJSON.h" // dumpJSON()

#include "nlohmann/json.hpp" // nlohmann::json

std::string
dumpJSON(const spin::Measurement & measurement)
{
  // Elements ordered as: n, r, k
  const std::size_t n = 0;
  const std::size_t r = 1;
  const std::size_t k = 2;

  const auto C = measurement.get_C();
  const auto Cmedian = measurement.get_C(true);
  const auto Cerr = measurement.get_CErr();
  const auto Bm = measurement.get_Bm();
  const auto Bmmedian = measurement.get_Bm(true);
  const auto Bmerr = measurement.get_BmErr();
  const auto Bp = measurement.get_Bp();
  const auto Bpmedian = measurement.get_Bp(true);
  const auto Bperr = measurement.get_BpErr();

  typedef std::map<std::string, double> sd_map;
  typedef std::map<std::string, sd_map> ssd_map;
  typedef std::map<std::string, ssd_map> sssd_map;

  nlohmann::json j = sssd_map{
    { "C", ssd_map{
      { "nominal", sd_map{
          { "nn", C.At(n, n) },
          { "nr", C.At(n, r) },
          { "nk", C.At(n, k) },
          { "rn", C.At(r, n) },
          { "rr", C.At(r, r) },
          { "rk", C.At(r, k) },
          { "kn", C.At(k, n) },
          { "kr", C.At(k, r) },
          { "kk", C.At(k, k) },
        },
      },
      { "median", sd_map{
          { "nn", Cmedian.At(n, n) },
          { "nr", Cmedian.At(n, r) },
          { "nk", Cmedian.At(n, k) },
          { "rn", Cmedian.At(r, n) },
          { "rr", Cmedian.At(r, r) },
          { "rk", Cmedian.At(r, k) },
          { "kn", Cmedian.At(k, n) },
          { "kr", Cmedian.At(k, r) },
          { "kk", Cmedian.At(k, k) },
        },
      },
      { "error", sd_map{
            { "nn", Cerr.At(n, n) },
            { "nr", Cerr.At(n, r) },
            { "nk", Cerr.At(n, k) },
            { "rn", Cerr.At(r, n) },
            { "rr", Cerr.At(r, r) },
            { "rk", Cerr.At(r, k) },
            { "kn", Cerr.At(k, n) },
            { "kr", Cerr.At(k, r) },
            { "kk", Cerr.At(k, k) },
          },
        },
      },
    },
    { "Bm", ssd_map{
        { "nominal", sd_map{
            { "n", Bm.At(n) },
            { "r", Bm.At(r) },
            { "k", Bm.At(k) },
          },
        },
        { "median", sd_map{
            { "n", Bmmedian.At(n) },
            { "r", Bmmedian.At(r) },
            { "k", Bmmedian.At(k) },
          },
        },
        { "error", sd_map{
            { "n", Bmerr.At(n) },
            { "r", Bmerr.At(r) },
            { "k", Bmerr.At(k) },
          },
        },
      },
    },
    { "Bp", ssd_map{
        { "nominal", sd_map{
            { "n", Bp.At(n) },
            { "r", Bp.At(r) },
            { "k", Bp.At(k) },
          },
        },
        { "median", sd_map{
            { "n", Bpmedian.At(n) },
            { "r", Bpmedian.At(r) },
            { "k", Bpmedian.At(k) },
          },
        },
        { "error", sd_map{
            { "n", Bperr.At(n) },
            { "r", Bperr.At(r) },
            { "k", Bperr.At(k) },
          },
        },
      },
    },
  };
  j["concurrence"] = sd_map{
    { "nominal", measurement.get_concurrence() },
    { "median",  measurement.get_concurrence(true) },
    { "error",   measurement.get_concurrenceErr() },
  };
  j["Ek"] = sd_map{
    { "nominal", measurement.get_Ek() },
    { "median",  measurement.get_Ek(true) },
    { "error",   measurement.get_EkErr() },
  };
  j["Rchsh"] = sd_map{
    { "nominal", measurement.get_Rchsh() },
    { "median",  measurement.get_Rchsh(true) },
    { "error",   measurement.get_RchshErr() },
  };
  j["steerability"] = sd_map{
    { "nominal", measurement.get_steerability() },
    { "median",  measurement.get_steerability(true) },
    { "error",   measurement.get_steerabilityErr() },
  };
  j["metadata"] = sd_map{
    { "sample_size",              measurement.get_sample_size() },
    { "num_bootstrap_samples",    measurement.get_num_bootstrap_samples() },
    { "boostrap_size",            measurement.get_boostrap_size() },
    { "num_concurrence_failures", measurement.get_num_concurrence_failures() },
  };
  return j.dump(2, ' ', true);
}
