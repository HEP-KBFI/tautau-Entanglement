#include "TauAnalysis/Entanglement/interface/format_vT.h"

#include <sstream> // std::ostringstream

template <typename T>
std::string
format_vT(const std::vector<T>& v)
{
  std::ostringstream os;
  os << "{ ";
  size_t numEntries = v.size();
  for ( size_t idx = 0; idx < numEntries; ++idx )
  {
    os << v[idx];
    if ( idx < (numEntries - 1) )
    {
      os << ", ";
    }
  }
  os << " }";
  return os.str();
}

std::string
format_vbool(const vbool& v)
{
  return format_vT(v);
}

std::string
format_vdouble(const vdouble& v)
{
  return format_vT(v);
}

std::string
format_vfloat(const vfloat& v)
{
  return format_vT(v);
}

std::string
format_vint(const vint& v)
{
  return format_vT(v);
}

std::string
format_vstring(const vstring& v)
{
  return format_vT(v);
}

std::string
format_vunsigned(const vunsigned& v)
{
  return format_vT(v);
}
