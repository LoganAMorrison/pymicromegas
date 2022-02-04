#ifndef PYMICROMEGAS_EWSB_HPP
#define PYMICROMEGAS_EWSB_HPP

#include "micromegas.hpp"

static constexpr double DEFAULT_MU = 350.0;
static constexpr double DEFAULT_MG1 = 200.0;
static constexpr double DEFAULT_MG2 = 400.0;
static constexpr double DEFAULT_MG3 = 800.0;
static constexpr double DEFAULT_ML1 = 500.0;
static constexpr double DEFAULT_ML2 = 500.0;
static constexpr double DEFAULT_ML3 = 500.0;
static constexpr double DEFAULT_MR1 = 200.0;
static constexpr double DEFAULT_MR2 = 200.0;
static constexpr double DEFAULT_MR3 = 200.0;
static constexpr double DEFAULT_MQ1 = 1000.0;
static constexpr double DEFAULT_MQ2 = 1000.0;
static constexpr double DEFAULT_MQ3 = 1000.0;
static constexpr double DEFAULT_MU1 = 300.0;
static constexpr double DEFAULT_MU2 = 300.0;
static constexpr double DEFAULT_MU3 = 300.0;
static constexpr double DEFAULT_MD1 = 300.0;
static constexpr double DEFAULT_MD2 = 300.0;
static constexpr double DEFAULT_MD3 = 300.0;
static constexpr double DEFAULT_AT = -1000.0;
static constexpr double DEFAULT_AB = 0.0;
static constexpr double DEFAULT_AL = 00;
static constexpr double DEFAULT_MH3 = 1000.0;
static constexpr double DEFAULT_TB = 10.0;
static constexpr double DEFAULT_AM = 0.0;

static constexpr double DEFAULT_ALFSMZ = 0.1184;
static constexpr double DEFAULT_MZ = 91.187;
static constexpr double DEFAULT_MW = 80.2;
static constexpr double DEFAULT_MTP = 173.07;
static constexpr double DEFAULT_MBMB = 4.23;
static constexpr double DEFAULT_MCMC = 1.27;
static constexpr double DEFAULT_GG = 1.117;
static constexpr double DEFAULT_Q = 100.0;
static constexpr double DEFAULT_EE = 0.3123;
static constexpr double DEFAULT_ML = 1.777;
static constexpr double DEFAULT_MQ = 0.05;

class EwsbParameters {
  using value_type = double;

  // Higgs Superpotential term: ℒ ⊃ μ Hu Hd
  value_type p_mu = DEFAULT_MU;

  // U(1)_Y gaugino (Bino) mass: ℒ ⊃ -1/2 M₁ B B
  value_type p_mg1 = DEFAULT_MG1;
  // SU(2) gaugino (Wino) mass: ℒ ⊃ -1/2 M₂ Wᵃ Wₐ
  value_type p_mg2 = DEFAULT_MG2;
  // SU(3) gaugino (Gluino) mass: ℒ ⊃ -1/2 M₃ Gᵃ Gₐ
  value_type p_mg3 = DEFAULT_MG3;

  // Mass of the 1st gen slepton doublet: ℒ ⊃ -ML₁ * ⟨L₁, L₁⟩
  value_type p_ml1 = DEFAULT_ML1;
  // Mass of the 2nd gen slepton doublet: ℒ ⊃ -ML₂ * ⟨L₂, L₂⟩
  value_type p_ml2 = DEFAULT_ML2;
  // Mass of the 3rd gen slepton doublet: ℒ ⊃ -ML₃ * ⟨L₃, L₃⟩
  value_type p_ml3 = DEFAULT_ML3;

  // Mass of the 1st gen RH slepton : ℒ ⊃ -Ml₁ * |lR₁|²
  value_type p_mr1 = DEFAULT_MR1;
  // Mass of the 2nd gen RH slepton : ℒ ⊃ -Ml₂ * |lR₂|²
  value_type p_mr2 = DEFAULT_MR2;
  // Mass of the 3rd gen RH slepton : ℒ ⊃ -Ml₃ * |lR₃|²
  value_type p_mr3 = DEFAULT_MR3;

  // Mass of the 1st gen squark doublet: ℒ ⊃ -MQ₁ * ⟨Q₁, Q₁⟩
  value_type p_mq1 = DEFAULT_MQ1;
  // Mass of the 2nd gen squark doublet: ℒ ⊃ -MQ₂ * ⟨Q₂, Q₂⟩
  value_type p_mq2 = DEFAULT_MQ2;
  // Mass of the 3rd gen squark doublet: ℒ ⊃ -MQ₃ * ⟨Q₃, Q₃⟩
  value_type p_mq3 = DEFAULT_MQ3;

  // Mass of the 1st gen RH up-squark : ℒ ⊃ -Mu₁ * |uR₁|²
  value_type p_mu1 = DEFAULT_MU1;
  // Mass of the 2nd gen RH up-squark : ℒ ⊃ -Mu₂ * |uR₂|²
  value_type p_mu2 = DEFAULT_MU2;
  // Mass of the 3rd gen RH up-squark : ℒ ⊃ -Mu₃ * |uR₃|²
  value_type p_mu3 = DEFAULT_MU3;

  // Mass of the 1st gen RH down-squark : ℒ ⊃ -Md₁ * |dR₁|²
  value_type p_md1 = DEFAULT_MD1;
  // Mass of the 2nd gen RH down-squark : ℒ ⊃ -Md₂ * |dR₂|²
  value_type p_md2 = DEFAULT_MD2;
  // Mass of the 3rd gen RH down-squark : ℒ ⊃ -Md₃ * |dR₃|²
  value_type p_md3 = DEFAULT_MD3;

  // Mass of the psuedo-scalar Higgs boson MA
  value_type p_mh3 = DEFAULT_MH3;
  // Tangent beta: tan(β) = v₂ / v₁
  value_type p_tb = DEFAULT_TB;

  // Trilinear top coupling: ℒ ⊃ At Yt tR ⟨Hu, Q₃⟩
  value_type p_at = DEFAULT_AT;
  // Trilinear bottom coupling: ℒ ⊃ Ab Yb bR ⟨Hd, Q₃⟩
  value_type p_ab = DEFAULT_AB;
  // Trilinear tau coupling: ℒ ⊃ Al Yl lR ⟨Hd, Q₃⟩
  value_type p_al = DEFAULT_AL;

  // Trilinear 1st/2nd gen up-quark couplings:
  //  ℒ ⊃ Au [ Yu uR ⟨Hu, Q₁⟩ + Yc cR ⟨Hu, Q₂⟩ ]
  value_type p_au = 0.0;
  // Trilinear 1st/2nd gen down-quark couplings:
  //  ℒ ⊃ Ad [ Yd dR ⟨Hd, Q₁⟩ + Ys sR ⟨Hd, Q₂⟩ ]
  value_type p_ad = 0.0;

  // alpha_QCD(MZ)
  value_type p_alfsmz = DEFAULT_ALFSMZ;
  // Z-boson mass.  PDG value   91.1876 +/- 0.0021
  value_type p_mz = DEFAULT_MZ;
  // W-boson mass.  PDG valeu   80.385  +/- 0.015
  value_type p_mw = DEFAULT_MW;
  // top quark pole mass
  value_type p_mtp = DEFAULT_MTP;
  // Mb(Mb) running mass
  value_type p_mbmb = DEFAULT_MBMB;
  // Mc(Mc) running mass
  value_type p_mcmc = DEFAULT_MCMC;
  // Drived by CalcHEP
  value_type p_gg = DEFAULT_GG;
  // QCD scale for running quark masses
  value_type p_q = DEFAULT_Q;
  // PDG value 4*pi/EE^2=128 corresponds to EE=0.3134
  value_type p_ee = DEFAULT_EE;
  // mass of tau-lepton
  value_type p_ml = DEFAULT_ML;
  // mass for light quarks
  value_type p_mq = DEFAULT_MQ;
  // No idea what this represents
  value_type p_am = DEFAULT_AM;

public:
  EwsbParameters() = default;

  explicit EwsbParameters(
      value_type mu, value_type mg1, value_type mg2,  // NOLINT
      value_type mg3, value_type ml1, value_type ml2, // NOLINT
      value_type ml3, value_type mr1, value_type mr2, // NOLINT
      value_type mr3, value_type mq1, value_type mq2, // NOLINT
      value_type mq3, value_type mu1, value_type mu2, // NOLINT
      value_type mu3, value_type md1, value_type md2, // NOLINT
      value_type md3, value_type mh3, value_type tb,  // NOLINT
      value_type at, value_type ab, value_type al)
      : p_mu(mu), p_mg1(mg1), p_mg2(mg2), p_mg3(mg3), p_ml1(ml1), p_ml2(ml2),
        p_ml3(ml3), p_mr1(mr1), p_mr2(mr2), p_mr3(mr3), p_mq1(mq1), p_mq2(mq2),
        p_mq3(mq3), p_mu1(mu1), p_mu2(mu2), p_mu3(mu3), p_md1(md1), p_md2(md2),
        p_md3(md3), p_mh3(mh3), p_tb(tb), p_at(at), p_ab(ab), p_al(al) {

    // p_alfsmz = micromegas::Micromegas::find_val("alfSMZ");
    // p_mz = micromegas::Micromegas::find_val("MZ");
    // p_mw = micromegas::Micromegas::find_val("MW");
    // p_mtp = micromegas::Micromegas::find_val("Mtp");
    // p_mbmb = micromegas::Micromegas::find_val("MbMb");
    // p_mcmc = micromegas::Micromegas::find_val("McMc");
    // p_gg = micromegas::Micromegas::find_val("GG");
    // p_q = micromegas::Micromegas::find_val("Q");
    // p_ee = micromegas::Micromegas::find_val("EE");
    // p_ml = micromegas::Micromegas::find_val("Ml");
    // p_mq = micromegas::Micromegas::find_val("Mq");
    // p_am = micromegas::Micromegas::find_val("Am");
    // p_au = micromegas::Micromegas::find_val("Au");
    // p_ad = micromegas::Micromegas::find_val("Ad");
  }

  void assign_all() const {
    micromegas::Micromegas::assign_val("mu", p_mu);
    micromegas::Micromegas::assign_val("MG1", p_mg1);
    micromegas::Micromegas::assign_val("MG2", p_mg2);
    micromegas::Micromegas::assign_val("MG3", p_mg3);
    micromegas::Micromegas::assign_val("Ml1", p_ml1);
    micromegas::Micromegas::assign_val("Ml2", p_ml2);
    micromegas::Micromegas::assign_val("Ml3", p_ml3);
    micromegas::Micromegas::assign_val("Mr1", p_mr1);
    micromegas::Micromegas::assign_val("Mr2", p_mr2);
    micromegas::Micromegas::assign_val("Mr3", p_mr3);
    micromegas::Micromegas::assign_val("Mq1", p_mq1);
    micromegas::Micromegas::assign_val("Mq2", p_mq2);
    micromegas::Micromegas::assign_val("Mq3", p_mq3);
    micromegas::Micromegas::assign_val("Mu1", p_mu1);
    micromegas::Micromegas::assign_val("Mu2", p_mu2);
    micromegas::Micromegas::assign_val("Mu3", p_mu3);
    micromegas::Micromegas::assign_val("Md1", p_md1);
    micromegas::Micromegas::assign_val("Md2", p_md2);
    micromegas::Micromegas::assign_val("Md3", p_md3);
    micromegas::Micromegas::assign_val("MH3", p_mh3);
    micromegas::Micromegas::assign_val("tb", p_tb);
    micromegas::Micromegas::assign_val("At", p_at);
    micromegas::Micromegas::assign_val("Ab", p_ab);
    micromegas::Micromegas::assign_val("Al", p_al);
  }

  // void initialize() const {
  //   const int err = micromegas::Micromegas::mssm_ewsb();
  //   if (err == 0) {
  //     // TODO Handle error
  //   }
  //   assign_all();
  // }

  [[nodiscard]] auto get_mu() const -> value_type { return p_mu; }
  [[nodiscard]] auto get_mg1() const -> value_type { return p_mg1; }
  [[nodiscard]] auto get_mg2() const -> value_type { return p_mg2; }
  [[nodiscard]] auto get_mg3() const -> value_type { return p_mg3; }
  [[nodiscard]] auto get_ml1() const -> value_type { return p_ml1; }
  [[nodiscard]] auto get_ml2() const -> value_type { return p_ml2; }
  [[nodiscard]] auto get_ml3() const -> value_type { return p_ml3; }
  [[nodiscard]] auto get_mr1() const -> value_type { return p_mr1; }
  [[nodiscard]] auto get_mr2() const -> value_type { return p_mr2; }
  [[nodiscard]] auto get_mr3() const -> value_type { return p_mr3; }
  [[nodiscard]] auto get_mq1() const -> value_type { return p_mq1; }
  [[nodiscard]] auto get_mq2() const -> value_type { return p_mq2; }
  [[nodiscard]] auto get_mq3() const -> value_type { return p_mq3; }
  [[nodiscard]] auto get_mu1() const -> value_type { return p_mu1; }
  [[nodiscard]] auto get_mu2() const -> value_type { return p_mu2; }
  [[nodiscard]] auto get_mu3() const -> value_type { return p_mu3; }
  [[nodiscard]] auto get_md1() const -> value_type { return p_md1; }
  [[nodiscard]] auto get_md2() const -> value_type { return p_md2; }
  [[nodiscard]] auto get_md3() const -> value_type { return p_md3; }
  [[nodiscard]] auto get_mh3() const -> value_type { return p_mh3; }
  [[nodiscard]] auto get_tb() const -> value_type { return p_tb; }
  [[nodiscard]] auto get_at() const -> value_type { return p_at; }
  [[nodiscard]] auto get_ab() const -> value_type { return p_ab; }
  [[nodiscard]] auto get_al() const -> value_type { return p_al; }
  [[nodiscard]] auto get_au() const -> value_type { return p_au; }
  [[nodiscard]] auto get_ad() const -> value_type { return p_ad; }
  [[nodiscard]] auto get_am() const -> value_type { return p_am; }
  [[nodiscard]] auto get_mtp() const -> value_type { return p_mtp; }

  [[nodiscard]] auto get_alfsmz() const -> value_type { return p_alfsmz; }
  [[nodiscard]] auto get_mbmb() const -> value_type { return p_mbmb; }
  [[nodiscard]] auto get_mcmc() const -> value_type { return p_mcmc; }

  [[nodiscard]] auto get_ee() const -> value_type { return p_mcmc; }
  [[nodiscard]] auto get_q() const -> value_type { return p_mcmc; }
  [[nodiscard]] auto get_ml() const -> value_type { return p_mcmc; }
  [[nodiscard]] auto get_mq() const -> value_type { return p_mcmc; }
  [[nodiscard]] auto get_gg() const -> value_type { return p_mcmc; }
  [[nodiscard]] auto get_mz() const -> value_type { return p_mcmc; }
  [[nodiscard]] auto get_mw() const -> value_type { return p_mcmc; }

  auto set_mu(double val) -> void { p_mu = val; }
  auto set_mg1(double val) -> void { p_mg1 = val; }
  auto set_mg2(double val) -> void { p_mg2 = val; }
  auto set_mg3(double val) -> void { p_mg3 = val; }
  auto set_ml1(double val) -> void { p_ml1 = val; }
  auto set_ml2(double val) -> void { p_ml2 = val; }
  auto set_ml3(double val) -> void { p_ml3 = val; }
  auto set_mr1(double val) -> void { p_mr1 = val; }
  auto set_mr2(double val) -> void { p_mr2 = val; }
  auto set_mr3(double val) -> void { p_mr3 = val; }
  auto set_mq1(double val) -> void { p_mq1 = val; }
  auto set_mq2(double val) -> void { p_mq2 = val; }
  auto set_mq3(double val) -> void { p_mq3 = val; }
  auto set_mu1(double val) -> void { p_mu1 = val; }
  auto set_mu2(double val) -> void { p_mu2 = val; }
  auto set_mu3(double val) -> void { p_mu3 = val; }
  auto set_md1(double val) -> void { p_md1 = val; }
  auto set_md2(double val) -> void { p_md2 = val; }
  auto set_md3(double val) -> void { p_md3 = val; }
  auto set_mh3(double val) -> void { p_mh3 = val; }
  auto set_tb(double val) -> void { p_tb = val; }
  auto set_at(double val) -> void { p_at = val; }
  auto set_ab(double val) -> void { p_ab = val; }
  auto set_al(double val) -> void { p_al = val; }
  auto set_au(double val) -> void { p_au = val; }
  auto set_ad(double val) -> void { p_ad = val; }
  auto set_am(double val) -> void { p_am = val; }
  auto set_mtp(double val) -> void { p_mtp = val; }
  auto set_alfsmz(double val) -> void { p_alfsmz = val; }
  auto set_mz(double val) -> void { p_mz = val; }
  auto set_mw(double val) -> void { p_mw = val; }
  auto set_mbmb(double val) -> void { p_mbmb = val; }
  auto set_mcmc(double val) -> void { p_mcmc = val; }
  auto set_ee(double val) -> void { p_ee = val; }
  auto set_q(double val) -> void { p_q = val; }
  auto set_ml(double val) -> void { p_ml = val; }
  auto set_mq(double val) -> void { p_mq = val; }
};

#endif // PYMICROMEGAS_EWSB_HPP
