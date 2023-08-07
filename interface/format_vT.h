#ifndef TauAnalysis_Entanglement_format_vT_h
#define TauAnalysis_Entanglement_format_vT_h

#include <string> // std::string
#include <vector> // std::vector<>

typedef std::vector<bool>        vbool;
typedef std::vector<double>      vdouble;
typedef std::vector<float>       vfloat;
typedef std::vector<int>         vint;
typedef std::vector<std::string> vstring;
typedef std::vector<unsigned>    vunsigned;

std::string
format_vbool(const vbool& v);

std::string
format_vdouble(const vdouble& v);

std::string
format_vfloat(const vfloat& v);

std::string
format_vint(const vint& v);

std::string
format_vstring(const vstring& v);

std::string
format_vunsigned(const vunsigned& v);

#endif // TauAnalysis_Entanglement_format_vT_h
