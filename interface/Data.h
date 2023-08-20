#ifndef TauAnalysis_Entanglement_Data_h
#define TauAnalysis_Entanglement_Data_h

namespace spin
{

class Data
{
 public:
  Data(float hPlus_r, float hPlus_n, float hPlus_k, 
       float hMinus_r, float hMinus_n, float hMinus_k,
       float evtWeight);
  ~Data();

  inline
  float
  get_hPlus_r() const
  {
    return hPlus_r_;
  }
  inline
  float
  get_hPlus_n() const
  {
    return hPlus_n_;
  }
  inline
  float
  get_hPlus_k() const
  {
    return hPlus_k_;
  }
  inline
  float
  get_hMinus_r() const
  {
    return hMinus_r_;
  }
  inline
  float
  get_hMinus_n() const
  {
    return hMinus_n_;
  }
  inline
  float
  get_hMinus_k() const
  {
    return hMinus_k_;
  }

  inline
  float
  get_evtWeight() const
  {
    return evtWeight_;
  }

 private:
  float hPlus_r_;
  float hPlus_n_;
  float hPlus_k_;
  float hMinus_r_;
  float hMinus_n_;
  float hMinus_k_;

  float evtWeight_;
};

}

#endif // TauAnalysis_Entanglement_Data_h
