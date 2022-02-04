#ifndef PYMICROMEGAS_SUGRA_HPP
#define PYMICROMEGAS_SUGRA_HPP

#include "micromegas.hpp"
#include <array>
#include <cmath>
#include <stdexcept>

class SugraParameters {

  static constexpr double p_DEFAULT_M0 = 120.0;
  static constexpr double p_DEFAULT_A0 = -350.0;
  static constexpr double p_DEFAULT_MHF = 500.0;
  static constexpr double p_DEFAULT_TB = 40.0;
  static constexpr double p_DEFAULT_SGN = 1.0;

  using value_type = double;
  value_type p_m0 = p_DEFAULT_M0;
  value_type p_mhf = p_DEFAULT_MHF;
  value_type p_a0 = p_DEFAULT_A0;
  value_type p_tb = p_DEFAULT_TB;
  value_type p_sgn = p_DEFAULT_SGN;

public:
  SugraParameters() = default;

  // NOLINTNEXTLINE
  SugraParameters(value_type m0, value_type mhf, value_type a0, value_type tb,
                  value_type sgn)
      : p_m0(m0), p_mhf(mhf), p_a0(a0), p_tb(tb), p_sgn(sgn) {}

  [[nodiscard]] auto get_m0() const -> value_type { return p_m0; }
  [[nodiscard]] auto get_mhf() const -> value_type { return p_mhf; }
  [[nodiscard]] auto get_a0() const -> value_type { return p_a0; }
  [[nodiscard]] auto get_tb() const -> value_type { return p_tb; }
  [[nodiscard]] auto get_sgn() const -> value_type { return p_sgn; }

  auto set_m0(double val) -> void { p_m0 = val; }
  auto set_mhf(double val) -> void { p_mhf = val; }
  auto set_a0(double val) -> void { p_a0 = val; }
  auto set_tb(double val) -> void { p_tb = val; }
  auto set_sgn(double val) -> void { p_sgn = val; }
};

#endif // PYMICROMEGAS_SUGRA_HPP
