#include "TauAnalysis/Entanglement/interface/comp_steerability.h"

using FuncType = std::function<double(double, double)>;

namespace {

  double
  integrate(FuncType integrand,
            unsigned n_theta = 200,
            unsigned n_phi = 200)
  {
    assert(n_theta > 1);
    assert(n_phi > 1);

    // Integration volume
    static const double V_theta = M_PI;
    static const double V_phi = 2 * M_PI;

    // Integration steps
    const double d_theta = V_theta / (n_theta - 1);
    const double d_phi = V_phi / (n_phi - 1);

    // Integrate
    double total = 0.;
    for(unsigned theta_idx = 0; theta_idx < n_theta; ++theta_idx)
    {
      const double theta = theta_idx * d_theta;
      double subtotal = 0.;
      for(unsigned phi_idx = 0; phi_idx < n_phi; ++phi_idx)
      {
        const double phi = phi_idx * d_phi;
        subtotal += integrand(theta, phi);
      }
      total += subtotal * std::sin(theta); // multiply subtotal with the Jacobian
    }

    // The result
    return total * V_theta * V_phi / (n_theta * n_phi);
  }

  double
  steerability_integrand(double theta,
                         double phi, const math::Matrix3x3& C)
  {
    const math::Vector3 unit_vector(
      std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), std::cos(theta)
    );
    const math::Vector3 scaled_vector = C * unit_vector;
    return ROOT::Math::Mag(scaled_vector);
  }
}

double
comp_steerability(const math::Matrix3x3& C,
                  unsigned n_theta,
                  unsigned n_phi)
{
  const auto steerability_bound = std::bind(::steerability_integrand, std::placeholders::_1, std::placeholders::_2, C);
  return integrate(steerability_bound, n_theta, n_phi) / (2 * M_PI);
}
