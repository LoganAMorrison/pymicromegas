
#ifndef MICROMEGAS_HPP
#define MICROMEGAS_HPP

#include <array>
#include <string>

namespace micromegas {

struct alignas(64) DirectDetectionAmps {
  std::array<double, 2> proton_si;
  std::array<double, 2> proton_sd;
  std::array<double, 2> neutron_si;
  std::array<double, 2> neutron_sd;
};

enum class DirectDetectionExperiment : unsigned int {
  Unknown = 0,
  Xenon1T2018 = 1,
  Cresst2019 = 2,
  DarkSide2018 = 4,
  Pico2019 = 8,
};

struct alignas(32) DirectDetectionResults {
  double xenon1T;
  double cresst;
  double darkside;
  double pico;
};

enum class MssmParameter {
  // Higgs Superpotential term: ℒ ⊃ μ Hu Hd
  Mu,
  // U(1)_Y gaugino (Bino) mass: ℒ ⊃ -1/2 M₁ B B
  Mg1,
  // SU(2) gaugino (Wino) mass: ℒ ⊃ -1/2 M₂ Wᵃ Wₐ
  Mg2,
  // SU(3) gaugino (Gluino) mass: ℒ ⊃ -1/2 M₃ Gᵃ Gₐ
  Mg3,
  // Mass of the 1st gen slepton doublet: ℒ ⊃ -ML₁ * ⟨L₁, L₁⟩
  Ml1,
  // Mass of the 2nd gen slepton doublet: ℒ ⊃ -ML₂ * ⟨L₂, L₂⟩
  Ml2,
  // Mass of the 3rd gen slepton doublet: ℒ ⊃ -ML₃ * ⟨L₃, L₃⟩
  Ml3,
  // Mass of the 1st gen RH slepton : ℒ ⊃ -Ml₁ * |lR₁|²
  Mr1,
  // Mass of the 2nd gen RH slepton : ℒ ⊃ -Ml₂ * |lR₂|²
  Mr2,
  // Mass of the 3rd gen RH slepton : ℒ ⊃ -Ml₃ * |lR₃|²
  Mr3,
  // Mass of the 1st gen squark doublet: ℒ ⊃ -MQ₁ * ⟨Q₁, Q₁⟩
  Mq1,
  // Mass of the 2nd gen squark doublet: ℒ ⊃ -MQ₂ * ⟨Q₂, Q₂⟩
  Mq2,
  // Mass of the 3rd gen squark doublet: ℒ ⊃ -MQ₃ * ⟨Q₃, Q₃⟩
  Mq3,
  // Mass of the 1st gen RH up-squark : ℒ ⊃ -Mu₁ * |uR₁|²
  Mu1,
  // Mass of the 2nd gen RH up-squark : ℒ ⊃ -Mu₂ * |uR₂|²
  Mu2,
  // Mass of the 3rd gen RH up-squark : ℒ ⊃ -Mu₃ * |uR₃|²
  Mu3,
  // Mass of the 1st gen RH down-squark : ℒ ⊃ -Md₁ * |dR₁|²
  Md1,
  // Mass of the 2nd gen RH down-squark : ℒ ⊃ -Md₂ * |dR₂|²
  Md2,
  // Mass of the 3rd gen RH down-squark : ℒ ⊃ -Md₃ * |dR₃|²
  Md3,
  // Mass of the psuedo-scalar Higgs boson MA
  Mh3,
  // Tangent beta: tan(β) = v₂ / v₁
  Tb,
  // Trilinear top coupling: ℒ ⊃ At Yt tR ⟨Hu, Q₃⟩
  At,
  // Trilinear bottom coupling: ℒ ⊃ Ab Yb bR ⟨Hd, Q₃⟩
  Ab,
  // Trilinear tau coupling: ℒ ⊃ Al Yl lR ⟨Hd, Q₃⟩
  Al,
  // Trilinear 1st/2nd gen up-quark couplings:
  //  ℒ ⊃ Au [ Yu uR ⟨Hu, Q₁⟩ + Yc cR ⟨Hu, Q₂⟩ ]
  Au,
  // Trilinear 1st/2nd gen down-quark couplings:
  //  ℒ ⊃ Ad [ Yd dR ⟨Hd, Q₁⟩ + Ys sR ⟨Hd, Q₂⟩ ]
  Ad,
  // alpha_QCD(MZ)
  Alfsmz,
  // Z-boson mass.  PDG value   91.1876 +/- 0.0021
  Mz,
  // W-boson mass.  PDG valeu   80.385  +/- 0.015
  Mw,
  // top quark pole mass
  Mtp,
  // Mb(Mb) running mass
  Mbmb,
  // Mc(Mc) running mass
  Mcmc,
  // Drived by CalcHEP
  Gg,
  // QCD scale for running quark masses
  Q,
  // PDG value 4*pi/EE^2=128 corresponds to EE=0.3134
  Ee,
  // mass of tau-lepton
  Ml,
  // mass for light quarks
  Mq,
  // No idea what this represents
  Am,
};

static const std::string XENON1T_2018_NAME = "XENON1T_2018";   // NOLINT
static const std::string CRESST_2019_NAME = "CRESST_2019";     // NOLINT
static const std::string DARKSIDE_2018_NAME = "DarkSide_2018"; // NOLINT
static const std::string PICO_2019_NAME = "PICO_2019";         // NOLINT

class Micromegas {
public:
  /**
   * Set the value of VZdecay and VWdecay.
   *
   * @param flag If 0, off-shell Z/W decays are turned off. If 1, 3-body
   * final-states are included. If 2, Z's and W's are included in
   * coannihilations as well.
   */
  static void set_v_decay(int);

  /**
   * Return the value of micromegas CDM1.
   */
  static std::string get_cdm1();

  /**
   * Return the value of micromegas CDM2.
   */
  static std::string get_cdm2();

  /**
   * Return the mass of micromegas first CDM.
   */
  static double get_mcdm1();

  /**
   * Return the mass of micromegas second CDM.
   */
  static double get_mcdm2();

  /**
   * Calculates the masses of Higgs and supersymmetric particles in the MSSM
   * including one-loop corrections starting from weak scale input parameters.
   */
  static int mssm_ewsb();

  /**
   * Derives the parameters at the electroweak symmetry breaking scale assuming
   * thate all input parameters except `tb` and `signMu` are defined at the GUT
   * scale.
   */
  static int mssm_sugra(double tb, double gMG1, double gMG2, double gMG3,
                        double gAl, double gAt, double gAb, double sgn,
                        double gMHu, double gMHd, double gMl1, double gMl3,
                        double gMr1, double gMr3, double gMq1, double gMq3,
                        double gMu1, double gMu3, double gMd1, double gMd3);

  /**
   * Assign value val to parameter name. Returns a non-zero
   * value if it cannot recognize a parameter name.
   */
  static int assign_val(const std::string &name, double val);

  /**
   * Assign value val to parameter name. Writes an error message if it cannot
   * recognize a parameter name.
   */
  static int assign_val_w(const std::string &name, double val);

  /**
   * Assign value val to parameter name. Returns a non-zero
   * value if it cannot recognize a parameter name.
   */
  static int mssm_assign_val(MssmParameter param, double val);

  /**
   * Assign value val to parameter name. Writes an error message if it cannot
   * recognize a parameter name.
   */
  static int mssm_assign_val_w(MssmParameter param, double val);

  /**
   * Finds the value of variable name and assigns it to parameter val. It
   * returns a non-zero value if it cannot recognize a parameter name.
   */
  static double find_val(const std::string &name);

  /**
   * Returns the value of variable name and writes an error message if it cannot
   * recognize a parameter name.
   */
  static double find_val_w(const std::string &fname);

  /**
   * Finds the value of variable name and assigns it to parameter val. It
   * returns a non-zero value if it cannot recognize a parameter name.
   */
  static double mssm_find_val(MssmParameter param);

  /**
   * Returns the value of variable name and writes an error message if it cannot
   * recognize a parameter name.
   */
  static double mssm_find_val_w(MssmParameter param);

  /**
   * Reads parameters from a file. The file should contain two columns with the
   * following format:
   *  name value
   *
   * Returns zero when the file has been read successfully, a negative value
   * when the file cannot be opened for reading and a positive value
   * corresponding to the line where a wrong file record was found.
   */
  static void read_var(const std::string &fname);

  // ===========================
  // ---- Utility Functions ----
  // ===========================

  /**
   * Sorts the odd particles with increasing masses. This routine fills the text
   * parameters CDM1 and CDM2 with the names of the lightest odd particle
   * starting with one and two tildes respectively and assigns the value of the
   * mass of the lightest odd particle in each sector to the global parameters
   * Mcdm1 and Mcdm2. For models with only one DM candidate, micrOMEGAs will set
   * CDM2=NULL and Mcdm2=0. This routine returns a non zero error code for a
   * wrong set of parameters, for example parameters for which some constraint
   * cannot be calculated. The name of the corresponding constraint is written
   * in txt. This routine has to be called after a reassignment of any input
   * parameter.
   */
  static void sort_odd_particles();

  /**
   * The sortOddParticles command which must be used to recompute the particle
   * spectrum after changing the model parameters also clears the decay table.
   */
  static void clean_decay_table();

  // ============================
  // ---- Particle Functions ----
  // ============================

  /**
   * Returns the PDG code.
   */
  static int name_to_pdg(const std::string &name);

  /**
   * Returns the numerical value of the particle mass.
   */
  static double name_to_mass(const std::string &name);

  /**
   * Returns the name of the particle which PDG code is nPDG. If this particle
   * does not exist in the model the return value is NULL.
   */
  static std::string pdg_to_name(int pdg);

  // =========================
  // ---- Pheno Functions ----
  // =========================

  /**
   * Returns the value of the supersymmetric contribution to the anomalous
   * magnetic moment of the muon.
   */
  static double mssm_gmuon();

  /**
   * Calculates the ∆ρ parameter in the MSSM. It contains for example the
   * stop/sbottom contributions, as well as the two-loop QCD corrections due to
   * gluon exchange and the correction due to gluino exchange in the heavy
   * gluino limit
   */
  static double mssm_deltarho();

  /**
   * Returns the value of the branching ratio for b → sγ, see Appendix A. We
   * have included some new contributions beyond the leading order that are
   * especially important for high tan β. SMbsg gives the SM contribution.
   */
  static std::pair<double, double> mssm_bsgnlo();

  /**
   * Returns the value of the branching ratio Bs → μ+μ− in the MSSM. It includes
   * the loop contributions due to chargino, sneutrino, stop and Higgs exchange.
   * The ∆mb effect relevant for high tan β is taken into account.
   */
  static double mssm_bsmumu();

  /**
   * Computes the ratio between the MSSM and SM branching fractions for  ̄B+ → τ
   * +ντ .
   */
  static double mssm_btaunu();

  /**
   * Computes the ratio of the MSSM to SM value for Rl23 in K+ → μν due to a
   * charged higgs contribution.
   */
  static double mssm_rl23();

  /**
   * Computes the branching ratio for D+s → τ +ντ . `dmunu` gives the branching
   * ratio for D+s →μ+νμ.
   */
  static std::pair<double, double> mssm_d_taunu_and_munu();

  /**
   * Returns a positive value and prints a WARNING when the choice of parameters
   * conflicts with a direct accelerator limits on sparticle masses from LEP.
   * The constraint on the light Higgs mass from the LHC is included.
   */
  static double mssm_masslimits();

  // =================================
  // ---- Relic Density Functions ----
  // =================================

  static std::pair<double, double> relic_density(bool fast = true,
                                                 double beps = 1e-4);

  // =================================
  // ---- Neutrino Flux Functions ----
  // =================================

  /**
   * Calculates the muon neutrino/anti-neutrino fluxes near the surface of the
   * Earth. Here f is the DM velocity distribution normalized such that
   *  ∫_{0}^{∞} v f (v)dv = 1.
   * The units are km/s for v and s2/km2 for f(v). For example, one can
   * use the same Maxwell function introduced for direct detection. This routine
   * implicitly depends on the WIMPSIM switch.
   */
  static std::pair<double, double> neutrino_flux(double (*vfv)(double),
                                                 int forSun);

  /**
   * calculates the exclusion confidence level for number of signal events
   * generated by given νμ and ̄νμ fluxes. The fluxes are assumed to be in [GeV
   * km2 Year]−1. This function uses the IC22BGdCos(cs) and IC22sigma(E) angular
   * distribution for background and signal as well as the event files
   * distributed by IceCube22 with φ < φcut = 8◦. The returned parameter B is a
   * Bayesian factor representing the ratio of likelihood functions for the
   * model with given fluxes and the model with null signal.
   */
  static double exclusion_confidence_icecube(double *nu, double *NU, double *B);

  // ====================================
  // ---- Direct Detection Functions ----
  // ====================================

  /**
   * Returns
   *  F(v) = cnorm / (2pi v Rot^2)^(3/2) * exp( - v^2 / vRot^2) theta(vEsc - v)
   * which corresponds to the isothermal model. Here vRot is the orbital
   * velocity of stars in the Milky Way, it is also a global parameter of
   * micrOMEGAs. cnorm is the normalization factor.
   */
  static double maxwell(double v);

  /**
   * Returns
   *    F(v) = (1 - eta) FM(v) + eta FS(v)
   * where FM is the standard Maxwell velocity distribution and the second
   * component is the velocity distribution from the Gaia sausage, which is not
   * spherically symmetric.
   */
  static double shmpp(double v);

  /**
   * Calculates the amplitudes for CDM-nucleon elastic scattering at zero
   * momentum. pAsi(nAsi) are spin independent amplitudes for protons(neutrons)
   * whereas pAsd(nAsd) are the cor- responding spin dependent amplitudes. Each
   * of these four parameters is an array of dimension 2. The zeroth (first)
   * element of these arrays gives the χ-nucleon amplitudes whereas the second
   * element gives χ-nucleon amplitudes. Amplitudes (in GeV−2) are normalized
   * such that the total cross section for either χ or χ cross sections is:
   *  sigma_tot = 4 mx^2 mn^2 / (pi * (mx + mn)^2) * (|Asi|^2 + 3 |Asd|^2)
   */
  static DirectDetectionAmps nucleon_amplitudes(char *WINP);

  /**
   * Same as `nucleon_amplitudes` but uses first dark-matter candidate.
   */
  static DirectDetectionAmps nucleon_amplitudes_cdm1();

  /**
   * Same as `nucleon_amplitudes` but uses second dark-matter candidate.
   */
  static DirectDetectionAmps nucleon_amplitudes_cdm2();

  /**
   * Returns the overall factor which should be applied to the cross sections,
   * σSIP , σSIN , σSDP , σSDN to reach the exclusion level α. All parameters
   * are the same as in DD_pvalCS.
   */
  static DirectDetectionResults
  direct_detection_factor_cs(double pval, double (*vfv)(double), double cs_SI_P,
                             double cs_SI_N, double cs_SD_P, double cs_SD_N);

  /**
   * Same as `direct_detection_factor_cs` but uses Maxwell distribution.
   */
  static DirectDetectionResults
  direct_detection_factor_cs_maxwell(double pval, double cs_SI_P,
                                     double cs_SI_N, double cs_SD_P,
                                     double cs_SD_N);

  /**
   * Same as `direct_detection_factor_cs` but uses SHM++ distribution.
   */
  static DirectDetectionResults
  direct_detection_factor_cs_shmpp(double pval, double cs_SI_P, double cs_SI_N,
                                   double cs_SD_P, double cs_SD_N);

  /**
   * calculates the value α = 1−C.L. for a model with DM-nucleon cross sections
   * σSIP , σSIN , σSDP , σSDN . Cross sections are specified in [pb] units. The
   * return value 0.1 corresponds to a 90% exclusion. The expCode parameter can
   * be any of the codes XENON1T_2018,DarkSide_2018, CRESST_2019,PICO_2019 or
   * their combination concatenated with the symbol |.
   */
  static DirectDetectionResults
  direct_detection_pval_cs(double (*vfv)(double), double cs_SI_P,
                           double cs_SI_N, double cs_SD_P, double cs_SD_N);

  /**
   * Same as `direct_detection_pval_cs` but uses Maxwell distribution.
   */
  static DirectDetectionResults
  direct_detection_pval_cs_maxwell(double cs_SI_P, double cs_SI_N,
                                   double cs_SD_P, double cs_SD_N);

  /**
   * Same as `direct_detection_pval_cs` but uses SHM++ distribution.
   */
  static DirectDetectionResults direct_detection_pval_cs_shmpp(double cs_SI_P,
                                                               double cs_SI_N,
                                                               double cs_SD_P,
                                                               double cs_SD_N);

  /**
   * Similar to DD_factorCS but use the cross section calculated
   * from the DM model under consideration in micrOMEGAs. The necessary
   * corrections for a light mediator are implemented automatically, these
   * functions do not use dNdEFact.
   */
  static DirectDetectionResults direct_detection_factor(double pval,
                                                        double (*vfv)(double));

  /**
   * Same as `direct_detection_factor` but uses Maxwell distribution.
   */
  static DirectDetectionResults direct_detection_factor_maxwell(double pval);

  /**
   * Same as `direct_detection_factor` but uses SMH++ distribution.
   */
  static DirectDetectionResults direct_detection_factor_shmpp(double pval);

  /**
   * Similar to DD_pvalCS but use the cross section calculated
   * from the DM model under consideration in micrOMEGAs. The necessary
   * corrections for a light mediator are implemented automatically, these
   * functions do not use dNdEFact.
   */
  static DirectDetectionResults direct_detection_pval(double (*vfv)(double));

  /**
   * Same as `direct_detection_factor` but uses Maxwell distribution.
   */
  static DirectDetectionResults direct_detection_pval_maxwell();

  /**
   * Same as `direct_detection_factor` but uses SHM++ distribution.
   */
  static DirectDetectionResults direct_detection_pval_shmpp();

  /**
   * Returns the excluded 90% SI cross section in cm^2 from the XENON1T
   * experiment. For a DM mass outside 6 < MDM < 1000 GeV, NaN is returned.
   */
  static double xenon1T_90(double M);

  /**
   * Returns the excluded 90% SD proton cross section in cm^2 from the XENON1T
   * experiment. For a DM mass outside 6 < MDM < 1000 GeV, NaN is returned.
   */
  static double xenon1T_sdp_90(double M);

  /**
   * Returns the excluded 90% SD neutron cross section in cm^2 from the XENON1T
   * experiment. For a DM mass outside 6 < MDM < 1000 GeV, NaN is returned.
   */
  static double xenon1T_sdn_90(double M);

  /**
   * Returns the excluded 90% SI cross section in cm^2 from the DarkSide
   * experiment. For a DM mass outside 0.7 < MDM < 15 GeV, NaN is returned.
   */
  static double darkside50_90(double M);

  static double darkside50_90_nob(double M);

  /**
   * Returns the excluded 90% SI cross section in cm^2 from the CRESST III
   * experiment. For a DM mass outside 0.35 < MDM < 12 GeV, NaN is returned.
   */
  static double cresst3_90(double M);

  /**
   * Returns the excluded 90% SD neutron cross section in cm^2 from the CRESST
   * III experiment. For a DM mass outside 0.35 < MDM < 12 GeV, NaN is returned.
   */
  static double cresst3_sdn_90(double M);

  /**
   * Returns the excluded 90% SI cross section in cm^2 from the PICO-60
   * experiment. For a DM mass outside 3 < MDM < 10000 GeV, NaN is returned.
   */
  static double pico60_90(double M);

  /**
   * Returns the excluded 90% SD proton cross section in cm^2 from the PICO-60
   * experiment. For a DM mass outside 3 < MDM < 10000 GeV, NaN is returned.
   */
  static double pico60_sdp_90(double M);

  // ==============================
  // ---- Constraint Functions ----
  // ==============================

  /**
   * Returns 1 and prints a WARNING if the invisible width of the Z boson of the
   * Standard Model is larger than 0.5 MeV and returns 0 if this
   * constraint is fulfilled. This routine can be used in any model with one DM
   * where the Z boson is uniquely defined by its PDG=23 and whether the neutral
   * LSP is its own antiparticle or not.
   */
  static bool z_invisible();

  /**
   * Returns 0 if the scenario considered is not excluded by Z′ constraints, 1
   * if the point is excluded and 2 if both subroutines dealing with Z′
   * constraints cannot test the given scenario. The routine can be used for any
   * Z′ uniquely defined by the PDG code 32. Currently two types of searches
   * defined in different subroutines of Zprimelimits() are implemented.
   * The 3.2/fb Z′ search in the dilepton final state at √s = 13 TeV from ATLAS
   * is considered in the first subroutine Zprime_dilepton. If the scenario
   * considered is allowed or not tested by Zprime_dilepton, a second subroutine
   * called Zprime_dijet analyses the point using constraints from LHC dijet
   * searches at √s = 8 TeV and at √s = 13 TeV. This subroutine uses the
   * recasting performed in [88] for a combination of ATLAS and CMS searches.
   * Zprimelimits() returns 1 if MZ′ < 0.5 TeV and 2 for points for which the
   * narrow-width approximation is not valid, i.e. Γ/MZ′ > 0.3.
   */
  static std::pair<double, double> lsp_nlsp_lep();

  /**
   * Checks the compatibility with the upper limit on the cross section for
   * the production of neutralinos σ(e+ e− →  ̃χ01  ̃χ0i ), i != 1, when the
   * heavier neutralino decays into quark pairs and the LSP,  ̃χ0i →  ̃χ01 + q +
   * q. The function returns σ × BR = ∑ i σ(e+e− → ̃χ01  ̃χ0i) × BR(χ0i →  ̃χ01 + q
   * + ̄q) in pb as well as a flag greater than one if σ × BR > 0.1(0.5) pb if
   * mNLSP > (<)100 GeV. This function can also be applied for non-SUSY models
   * which feature the same signature, in this case the function will compute
   * the cross section for production of the LSP and any other neutral particle
   * from the odd sector which can decay into the LSP and a Z boson.
   */
  static bool z_prime_limits();

  /**
   * Computes the cross section for p, p → pName1, pName2 + jet at √s = 8 TeV
   * where pName1, pName2 are the names of neutral outgoing particles and jet
   * includes light quarks (u,d,s) and gluons. The function returns the
   * resulting confidence level obtained with the CLs technique for each signal
   * region of the 8 TeV CMS mono-jet analysis and chooses the most
   * constraining one.
   */
  static double monojet();
};

} // namespace micromegas

#endif // MICROMEGAS_HPP
