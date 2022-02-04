#ifndef PYMICROMEGAS_SETTINGS_HPP
#define PYMICROMEGAS_SETTINGS_HPP

#include <bitset>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

// alignas(32)
class MicromegasSettings {
  bool p_relic_density = true;
  bool p_masses = true;
  bool p_gmuon = true;
  bool p_bsg = false;
  bool p_bsmumu = false;
  bool p_btaunu = false;
  bool p_deltarho = false;
  bool p_rl23 = false;
  bool p_d_taunu_and_munu = false;
  bool p_masslimits = false;
  bool p_nucleon_amplitudes = false;
  bool p_direct_detection_pvalues = false;
  bool p_z_invisible = false;
  bool p_lsp_nlsp_lep = false;
  bool p_z_prime_limits = false;
  bool p_monojet = false;
  bool p_fast = true;
  double p_beps = 1e-4;

public:
  explicit MicromegasSettings(
      bool t_relic_density = true, bool t_masses = true, bool t_gmuon = true,
      bool t_bsg = false, bool t_bsmumu = false, bool t_btaunu = false,
      bool t_deltarho = false, bool t_rl23 = false,
      bool t_d_taunu_and_munu = false, bool t_masslimits = false,
      bool t_nucleon_amplitudes = false,
      bool t_direct_detection_pvalues = false, bool t_z_invisible = false,
      bool t_lsp_nlsp_lep = false, bool t_z_prime_limits = false,
      bool t_monojet = false, bool t_fast = true, double t_beps = 1e-4)
      : p_relic_density(t_relic_density), p_masses(t_masses), p_gmuon(t_gmuon),
        p_bsg(t_bsg), p_bsmumu(t_bsmumu), p_btaunu(t_btaunu),
        p_deltarho(t_deltarho), p_rl23(t_rl23),
        p_d_taunu_and_munu(t_d_taunu_and_munu), p_masslimits(t_masslimits),
        p_nucleon_amplitudes(t_nucleon_amplitudes),
        p_direct_detection_pvalues(t_direct_detection_pvalues),
        p_z_invisible(t_z_invisible), p_lsp_nlsp_lep(t_lsp_nlsp_lep),
        p_z_prime_limits(t_z_prime_limits), p_monojet(t_monojet),
        p_fast(t_fast), p_beps(t_beps) {}

  // clang-format off
  [[nodiscard]] auto get_relic_density() const -> const bool& {return p_relic_density ;}
  [[nodiscard]] auto get_masses() const -> const bool& {return p_masses ;}
  [[nodiscard]] auto get_gmuon() const -> const bool& {return p_gmuon ;}
  [[nodiscard]] auto get_bsg() const -> const bool& {return p_bsg ;}
  [[nodiscard]] auto get_bsmumu() const -> const bool& {return p_bsmumu ;}
  [[nodiscard]] auto get_btaunu() const -> const bool&{return p_btaunu ;}
  [[nodiscard]] auto get_deltarho() const -> const bool& {return p_deltarho ;}
  [[nodiscard]] auto get_rl23() const -> const bool& {return p_rl23 ;}
  [[nodiscard]] auto get_d_taunu_and_munu() const -> const bool& {return p_d_taunu_and_munu ;}
  [[nodiscard]] auto get_masslimits() const -> const bool& {return p_masslimits ;}
  [[nodiscard]] auto get_nucleon_amplitudes() const -> const bool& {return p_nucleon_amplitudes ;}
  [[nodiscard]] auto get_direct_detection_pvalues() const -> const bool& {return p_direct_detection_pvalues ;}
  [[nodiscard]] auto get_z_invisible() const -> const bool& {return p_z_invisible ;}
  [[nodiscard]] auto get_lsp_nlsp_lep() const -> const bool& {return p_lsp_nlsp_lep ;}
  [[nodiscard]] auto get_z_prime_limits() const -> const bool& {return p_z_prime_limits ;}
  [[nodiscard]] auto get_monojet() const -> const bool& {return p_monojet ;}
  [[nodiscard]] auto get_fast() const -> const bool& {return p_fast ;}
  [[nodiscard]] auto get_beps() const -> const double& {return p_beps ;}

  auto set_relic_density(bool val) -> void {p_relic_density=val ;}
  auto set_masses(bool val) -> void {p_masses=val ;}
  auto set_gmuon(bool val) -> void {p_gmuon=val ;}
  auto set_bsg(bool val) -> void {p_bsg=val ;}
  auto set_bsmumu(bool val) -> void {p_bsmumu=val ;}
  auto set_btaunu(bool val) -> void {p_btaunu=val ;}
  auto set_deltarho(bool val) -> void {p_deltarho=val ;}
  auto set_rl23(bool val) -> void {p_rl23=val ;}
  auto set_d_taunu_and_munu(bool val) -> void {p_d_taunu_and_munu=val ;}
  auto set_masslimits(bool val) -> void {p_masslimits=val ;}
  auto set_nucleon_amplitudes(bool val) -> void {p_nucleon_amplitudes=val ;}
  auto set_direct_detection_pvalues(bool val) -> void {p_direct_detection_pvalues=val ;}
  auto set_z_invisible(bool val) -> void {p_z_invisible=val ;}
  auto set_lsp_nlsp_lep(bool val) -> void {p_lsp_nlsp_lep=val ;}
  auto set_z_prime_limits(bool val) -> void {p_z_prime_limits=val ;}
  auto set_monojet(bool val) -> void {p_monojet=val ;}
  auto set_fast(bool val) -> void {p_fast=val ;}
  auto set_beps(double val) -> void {p_beps=val ;}

  // clang-format on
};

#endif // PYMICROMEGAS_SETTINGS_HPP
