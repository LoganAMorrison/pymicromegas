#ifndef PYMICROMEGAS_RESULTS_HPP
#define PYMICROMEGAS_RESULTS_HPP

#include "ewsb.hpp"
#include "settings.hpp"
#include "sugra.hpp"
#include <string>
#include <vector>

constexpr size_t NUM_MASSES = 31; // Electron-snuetrino
constexpr size_t IDX_MSNE = 0;    // Electron-snuetrino
constexpr size_t IDX_MSNM = 1;    // Muon-Snuetrino
constexpr size_t IDX_MSEL = 2;    // Left-handed Selectron
constexpr size_t IDX_MSER = 3;    // Right-handed Selectron
constexpr size_t IDX_MSML = 4;    // Left-handed Smuon
constexpr size_t IDX_MSMR = 5;    // Right-handed Smuon
constexpr size_t IDX_MSDL = 6;    // Left-handed down squark
constexpr size_t IDX_MSDR = 7;    // Right-handed up squark
constexpr size_t IDX_MSUL = 8;    // Left-handed up squark
constexpr size_t IDX_MSUR = 9;    // Right-handed up squark
constexpr size_t IDX_MSSL = 10;   // Left-handed strange squark
constexpr size_t IDX_MSSR = 11;   // Right-handed strange squark
constexpr size_t IDX_MSCL = 12;   // Left-handed charm squark
constexpr size_t IDX_MSCR = 13;   // Right-handed charm squark
constexpr size_t IDX_MSNL = 14;   // Left-handed tau snuetrino
constexpr size_t IDX_MSL1 = 15;   // Stau 1
constexpr size_t IDX_MSL2 = 16;   // Stau 2
constexpr size_t IDX_MSB1 = 17;   // Bottom squark 1
constexpr size_t IDX_MSB2 = 18;   // Bottom squark 2
constexpr size_t IDX_MST1 = 19;   // Top squark 1
constexpr size_t IDX_MST2 = 20;   // Top squark 2
constexpr size_t IDX_MG = 21;     // Gluino mass
constexpr size_t IDX_MNE1 = 22;   // Neutralino 1
constexpr size_t IDX_MNE2 = 23;   // Neutralino 2
constexpr size_t IDX_MNE3 = 24;   // Neutralino 3
constexpr size_t IDX_MNE4 = 25;   // Neutralino 4
constexpr size_t IDX_MC1 = 26;    // Chargino 1
constexpr size_t IDX_MC2 = 27;    // Chargino 2
constexpr size_t IDX_MHSM = 28;   // Higgs mass
constexpr size_t IDX_MH = 29;     // Other higgs mass
constexpr size_t IDX_MHC = 30;    // Charged Higgs mass

class MicromegasResults {
  using value_type = std::vector<double>;
  value_type p_omega{};
  value_type p_xf{};
  std::array<value_type, NUM_MASSES> p_masses;

  value_type p_gmuon{};
  value_type p_bsgsm{};
  value_type p_bsgnlo{};
  value_type p_bsmumu{};
  value_type p_btaunu{};
  value_type p_deltarho{};
  value_type p_rl23{};
  value_type p_dtaunu{};
  value_type p_dmunu{};
  value_type p_masslimits{};
  value_type p_proton_si_amp{};
  value_type p_neutron_si_amp{};
  value_type p_proton_sd_amp{};
  value_type p_neutron_sd_amp{};
  value_type p_pval_xenon1T{};
  value_type p_pval_darkside{};
  value_type p_pval_pico{};
  value_type p_pval_cresst{};
  value_type p_z_invisible{};
  value_type p_lsp_nlsp_lep{};
  value_type p_z_prime_limits{};
  value_type p_monojet{};

public:
  explicit MicromegasResults(size_t n) {
    p_omega.reserve(n);
    p_xf.reserve(n);
    p_gmuon.reserve(n);
    p_bsgsm.reserve(n);
    p_bsgnlo.reserve(n);
    p_bsmumu.reserve(n);
    p_btaunu.reserve(n);
    p_deltarho.reserve(n);
    p_rl23.reserve(n);
    p_dtaunu.reserve(n);
    p_dmunu.reserve(n);
    p_masslimits.reserve(n);
    p_proton_si_amp.reserve(n);
    p_neutron_si_amp.reserve(n);
    p_proton_sd_amp.reserve(n);
    p_neutron_sd_amp.reserve(n);
    p_pval_xenon1T.reserve(n);
    p_pval_darkside.reserve(n);
    p_pval_pico.reserve(n);
    p_pval_cresst.reserve(n);
    p_z_invisible.reserve(n);
    p_lsp_nlsp_lep.reserve(n);
    p_z_prime_limits.reserve(n);
    p_monojet.reserve(n);

    p_masses[IDX_MSNE].reserve(n);
    p_masses[IDX_MSNM].reserve(n);
    p_masses[IDX_MSEL].reserve(n);
    p_masses[IDX_MSER].reserve(n);
    p_masses[IDX_MSML].reserve(n);
    p_masses[IDX_MSMR].reserve(n);
    p_masses[IDX_MSDL].reserve(n);
    p_masses[IDX_MSDR].reserve(n);
    p_masses[IDX_MSUL].reserve(n);
    p_masses[IDX_MSUR].reserve(n);
    p_masses[IDX_MSSL].reserve(n);
    p_masses[IDX_MSSR].reserve(n);
    p_masses[IDX_MSCL].reserve(n);
    p_masses[IDX_MSCR].reserve(n);
    p_masses[IDX_MSNL].reserve(n);
    p_masses[IDX_MSL1].reserve(n);
    p_masses[IDX_MSL2].reserve(n);
    p_masses[IDX_MSB1].reserve(n);
    p_masses[IDX_MSB2].reserve(n);
    p_masses[IDX_MST1].reserve(n);
    p_masses[IDX_MST2].reserve(n);
    p_masses[IDX_MG].reserve(n);
    p_masses[IDX_MNE1].reserve(n);
    p_masses[IDX_MNE2].reserve(n);
    p_masses[IDX_MNE3].reserve(n);
    p_masses[IDX_MNE4].reserve(n);
    p_masses[IDX_MC1].reserve(n);
    p_masses[IDX_MC2].reserve(n);
    p_masses[IDX_MHSM].reserve(n);
    p_masses[IDX_MH].reserve(n);
    p_masses[IDX_MHC].reserve(n);
  }

  void compute_relic_density(const MicromegasSettings &);
  void compute_masses(const MicromegasSettings &);
  void compute_gmuon(const MicromegasSettings &);
  void compute_bsg(const MicromegasSettings &);
  void compute_bsmumu(const MicromegasSettings &);
  void compute_btaunu(const MicromegasSettings &);
  void compute_deltarho(const MicromegasSettings &);
  void compute_rl23(const MicromegasSettings &);
  void compute_d_taunu_and_munu(const MicromegasSettings &);
  void compute_masslimits(const MicromegasSettings &);
  void compute_nucleon_amplitudes(const MicromegasSettings &);
  void compute_direct_detection_pvalues(const MicromegasSettings &);
  void compute_z_invisible(const MicromegasSettings &);
  void compute_lsp_nlsp_lep(const MicromegasSettings &);
  void compute_z_prime_limits(const MicromegasSettings &);
  void compute_monojet(const MicromegasSettings &);

  void set_nans();

  // clang-format off
 [[nodiscard]] auto omega() const -> const value_type & { return p_omega; }
 [[nodiscard]] auto xf() const -> const value_type & { return p_xf; }
 [[nodiscard]] auto delta_rho() const -> const value_type & { return p_deltarho; }
 [[nodiscard]] auto b_sg_sm() const -> const value_type & { return p_bsgsm; }
 [[nodiscard]] auto b_sg_nlo() const -> const value_type & { return p_bsgnlo; }
 [[nodiscard]] auto b_smumu() const -> const value_type & { return p_bsmumu; }
 [[nodiscard]] auto b_taunu() const -> const value_type & { return p_btaunu; }
 [[nodiscard]] auto g_muon() const -> const value_type & { return p_gmuon; }
 [[nodiscard]] auto rl23() const -> const value_type & { return p_rl23; }
 [[nodiscard]] auto dtaunu() const -> const value_type & { return p_dtaunu; }
 [[nodiscard]] auto dmunu() const -> const value_type & { return p_dmunu; }
 [[nodiscard]] auto masslimits() const -> const value_type & { return p_masslimits; }
 [[nodiscard]] auto proton_si_amp() const -> const value_type & { return p_proton_si_amp; }
 [[nodiscard]] auto proton_sd_amp() const -> const value_type & { return p_proton_sd_amp; }
 [[nodiscard]] auto neutron_si_amp() const -> const value_type & { return p_neutron_si_amp; }
 [[nodiscard]] auto neutron_sd_amp() const -> const value_type & { return p_neutron_sd_amp; }
 [[nodiscard]] auto pval_xenon1T() const -> const value_type & { return p_pval_xenon1T; }
 [[nodiscard]] auto pval_cresst() const -> const value_type & { return p_pval_cresst; }
 [[nodiscard]] auto pval_darkside() const -> const value_type & { return p_neutron_si_amp; }
 [[nodiscard]] auto pval_pico() const -> const value_type & { return p_neutron_sd_amp; }
 [[nodiscard]] auto msne() const -> const value_type & { return p_masses[IDX_MSNE]; }
 [[nodiscard]] auto msnm() const -> const value_type & { return p_masses[IDX_MSNM]; }
 [[nodiscard]] auto msel() const -> const value_type & { return p_masses[IDX_MSEL]; }
 [[nodiscard]] auto mser() const -> const value_type & { return p_masses[IDX_MSER]; }
 [[nodiscard]] auto msml() const -> const value_type & { return p_masses[IDX_MSML]; }
 [[nodiscard]] auto msmr() const -> const value_type & { return p_masses[IDX_MSMR]; }
 [[nodiscard]] auto msdl() const -> const value_type & { return p_masses[IDX_MSDL]; }
 [[nodiscard]] auto msdr() const -> const value_type & { return p_masses[IDX_MSDR]; }
 [[nodiscard]] auto msul() const -> const value_type & { return p_masses[IDX_MSUL]; }
 [[nodiscard]] auto msur() const -> const value_type & { return p_masses[IDX_MSUR]; }
 [[nodiscard]] auto mssl() const -> const value_type & { return p_masses[IDX_MSSL]; }
 [[nodiscard]] auto mssr() const -> const value_type & { return p_masses[IDX_MSSR]; }
 [[nodiscard]] auto mscl() const -> const value_type & { return p_masses[IDX_MSCL]; }
 [[nodiscard]] auto mscr() const -> const value_type & { return p_masses[IDX_MSCR]; }
 [[nodiscard]] auto msnl() const -> const value_type & { return p_masses[IDX_MSNL]; }
 [[nodiscard]] auto msl1() const -> const value_type & { return p_masses[IDX_MSL1]; }
 [[nodiscard]] auto msl2() const -> const value_type & { return p_masses[IDX_MSL2]; }
 [[nodiscard]] auto msb1() const -> const value_type & { return p_masses[IDX_MSB1]; }
 [[nodiscard]] auto msb2() const -> const value_type & { return p_masses[IDX_MSB2]; }
 [[nodiscard]] auto mst1() const -> const value_type & { return p_masses[IDX_MST1]; }
 [[nodiscard]] auto mst2() const -> const value_type & { return p_masses[IDX_MST2]; }
 [[nodiscard]] auto mg() const -> const value_type & { return p_masses[IDX_MG]; }
 [[nodiscard]] auto mneut1() const -> const value_type & { return p_masses[IDX_MNE1]; }
 [[nodiscard]] auto mneut2() const -> const value_type & { return p_masses[IDX_MNE2]; }
 [[nodiscard]] auto mneut3() const -> const value_type & { return p_masses[IDX_MNE3]; }
 [[nodiscard]] auto mneut4() const -> const value_type & { return p_masses[IDX_MNE4]; }
 [[nodiscard]] auto mchg1() const -> const value_type & { return p_masses[IDX_MC1]; }
 [[nodiscard]] auto mchg2() const -> const value_type & { return p_masses[IDX_MC2]; }
 [[nodiscard]] auto mhsm() const -> const value_type & { return p_masses[IDX_MHSM]; }
 [[nodiscard]] auto mh() const -> const value_type & { return p_masses[IDX_MH]; }
 [[nodiscard]] auto mhc() const -> const value_type & { return p_masses[IDX_MHC]; }
  // clang-format on

  auto omega() -> value_type & { return p_omega; }
  auto xf() -> value_type & { return p_xf; }
  auto bsgsm() -> value_type & { return p_bsgsm; }
  auto bsgnlo() -> value_type & { return p_bsgnlo; }
  auto deltarho() -> value_type & { return p_deltarho; }
  auto bsmumu() -> value_type & { return p_bsmumu; }
  auto btaunu() -> value_type & { return p_btaunu; }
  auto gmuon() -> value_type & { return p_gmuon; }
  auto g_muon() -> value_type & { return p_gmuon; }
  auto rl23() -> value_type & { return p_rl23; }
  auto dtaunu() -> value_type & { return p_dtaunu; }
  auto dmunu() -> value_type & { return p_dmunu; }
  auto masslimits() -> value_type & { return p_masslimits; }
  auto proton_si_amp() -> value_type & { return p_proton_si_amp; }
  auto proton_sd_amp() -> value_type & { return p_proton_sd_amp; }
  auto neutron_si_amp() -> value_type & { return p_neutron_si_amp; }
  auto neutron_sd_amp() -> value_type & { return p_neutron_sd_amp; }
  auto pval_xenon1T() -> value_type & { return p_pval_xenon1T; }
  auto pval_cresst() -> value_type & { return p_pval_cresst; }
  auto pval_darkside() -> value_type & { return p_neutron_si_amp; }
  auto pval_pico() -> value_type & { return p_neutron_sd_amp; }
  auto msne() -> value_type & { return p_masses[IDX_MSNE]; }
  auto msnm() -> value_type & { return p_masses[IDX_MSNM]; }
  auto msel() -> value_type & { return p_masses[IDX_MSEL]; }
  auto mser() -> value_type & { return p_masses[IDX_MSER]; }
  auto msml() -> value_type & { return p_masses[IDX_MSML]; }
  auto msmr() -> value_type & { return p_masses[IDX_MSMR]; }
  auto msdl() -> value_type & { return p_masses[IDX_MSDL]; }
  auto msdr() -> value_type & { return p_masses[IDX_MSDR]; }
  auto msul() -> value_type & { return p_masses[IDX_MSUL]; }
  auto msur() -> value_type & { return p_masses[IDX_MSUR]; }
  auto mssl() -> value_type & { return p_masses[IDX_MSSL]; }
  auto mssr() -> value_type & { return p_masses[IDX_MSSR]; }
  auto mscl() -> value_type & { return p_masses[IDX_MSCL]; }
  auto mscr() -> value_type & { return p_masses[IDX_MSCR]; }
  auto msnl() -> value_type & { return p_masses[IDX_MSNL]; }
  auto msl1() -> value_type & { return p_masses[IDX_MSL1]; }
  auto msl2() -> value_type & { return p_masses[IDX_MSL2]; }
  auto msb1() -> value_type & { return p_masses[IDX_MSB1]; }
  auto msb2() -> value_type & { return p_masses[IDX_MSB2]; }
  auto mst1() -> value_type & { return p_masses[IDX_MST1]; }
  auto mst2() -> value_type & { return p_masses[IDX_MST2]; }
  auto mg() -> value_type & { return p_masses[IDX_MG]; }
  auto mneut1() -> value_type & { return p_masses[IDX_MNE1]; }
  auto mneut2() -> value_type & { return p_masses[IDX_MNE2]; }
  auto mneut3() -> value_type & { return p_masses[IDX_MNE3]; }
  auto mneut4() -> value_type & { return p_masses[IDX_MNE4]; }
  auto mchg1() -> value_type & { return p_masses[IDX_MC1]; }
  auto mchg2() -> value_type & { return p_masses[IDX_MC2]; }
  auto mhsm() -> value_type & { return p_masses[IDX_MHSM]; }
  auto mh() -> value_type & { return p_masses[IDX_MH]; }
  auto mhc() -> value_type & { return p_masses[IDX_MHC]; }
};

#endif // PYMICROMEGAS_RESULTS_HPP
