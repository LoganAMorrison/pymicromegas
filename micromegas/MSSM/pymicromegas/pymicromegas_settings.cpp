#include "pymicromegas.hpp"

namespace py = pybind11;

void define_settings(pybind11::class_<MicromegasSettings> *settings) {

  // clang-format off
  static const std::string settings_doc_relic_density = "If true, the DM relic density and DM freeze-out are computed.";
  static const std::string settings_doc_masses = "If true, the SUSY particle masses are queried.";
  static const std::string settings_doc_gmuon = "If true, the g-2 of the muon is computed.";
  static const std::string settings_doc_bsg = "If true, micrOMEGAs will compute the branching ratio for b → s + γ.";
  static const std::string settings_doc_bsmumu = "If true, micrOMEGAs will compute the branching ratio for Bs → μ⁺ + μ⁻.";
  static const std::string settings_doc_btaunu = "If true, micrOMEGAs will compute the branching ratio for B⁺ → τ⁺ + ντ.";
  static const std::string settings_doc_deltarho = "If true, micrOMEGAs will compute the value of Δρ.";
  static const std::string settings_doc_rl23 = "If true, micrOMEGAs will compute the value of Rl23 in K⁺ → μ⁺ + νμ.";
  static const std::string settings_doc_d_taunu_and_munu = "If true, micrOMEGAs will compute the branching ratios for Ds⁺ → τ⁺ + ντ\nand Ds⁺ → μ⁺ + νμ.";
  static const std::string settings_doc_masslimits = "If true, micrOMEGAs will determine if sparticle masses are in conflict\nwith LEP and if the Higgs mass is in conflict with LHC.";
  static const std::string settings_doc_nucleon_amplitudes = "If true, micrOMEGAs will compute the DM-proton and DM-neutron\nspin-independent and spin-dependent scattering amplitudes are computed.";
  static const std::string settings_doc_direct_detection_pvalues = "If true, micrOMEGAs will compute the p-values for the Xenon1T,\nCRESST-III, PICO and DarkSide direct-detection experiments are\ncomputed.";
  static const std::string settings_doc_z_invisible = "If true, micrOMEGAs will determine if the invisible width of the Z boson\nis larger that 0.5 MeV.";
  static const std::string settings_doc_lsp_nlsp_lep = "If true, micrOMEGAs will check compatibility with the upper limit of\nσ(e⁺ + e⁻ → χ01 + χ0i) (where χ0i decays to χ01 + q + q) is violated\nby LEP.";
  static const std::string settings_doc_z_prime_limits = "If true, micrOMEGAs will check if model is excluded by Z' experiments.";
  static const std::string settings_doc_monojet = "If true, micrOMEGAs will compute confidence-limit on\nσ(proton + proton → N1 + N2 + jet) from an 8 TeV CMS mono-jet analysis.";
  static const std::string settings_doc_fast = "If true, micrOMEGAs will compute relic density using a approximation\naccurate to about 2% rather than use a more accurate (but much slower)\nmethod.";
  static const std::string settings_doc_beps = "The criteria for including a given coannihilation channel.\nThe recommended value is 1e-4 <= beps <= 1e-6. beps=1 means only\nannihilation with the lightest odd particle are computed.";
  static const std::string settings_doc_debug = "If true, debug/error messages will be printed.";
  // clang-format on

  static const std::string settings_init_doc = R"pbdoc(
    Construct a settings object used to pass to one of micrOMEGAs' RGEs. This
    object is used to tell micrOMEGAs what to compute.

    Parameters
    ----------
    masses: bool, optional
        If true, the SUSY particle masses are queried.
        Default is true.
    gmuon: bool, optional = true;
        If true, the g-2 of the muon is computed.
        Default is true.
    bsg: bool, optional = false;
        If true, micrOMEGAs will compute the branching ratio for b → s + γ.
        Default is false.
    bsmumu: bool, optional = false;
        If true, micrOMEGAs will compute the branching ratio for Bs → μ⁺ + μ⁻.
        Default is false.
    btaunu: bool, optional = false;
        If true, micrOMEGAs will compute the branching ratio for B⁺ → τ⁺ + ντ.
        Default is false.
    deltarho: bool, optional = false;
        If true, micrOMEGAs will compute the value of Δρ. Default is false.
    rl23: bool, optional = false;
        If true, micrOMEGAs will compute the value of Rl23 in K⁺ → μ⁺ + νμ.
        Default is false.
    d_taunu_and_munu: bool, optional = false;
        If true, micrOMEGAs will compute the branching ratios for Ds⁺ → τ⁺ + ντ
        and Ds⁺ → μ⁺ + νμ. Default is false.
    masslimits: bool, optional = false;
        If true, micrOMEGAs will determine if sparticle masses are in conflict
        with LEP and if the Higgs mass is in conflict with LHC. Default is false.
    nucleon_amplitudes: bool, optional = false;
        If true, micrOMEGAs will compute the DM-proton and DM-neutron
        spin-independent and spin-dependent scattering amplitudes are computed.
        Default is false.
    direct_detection_pvalues: bool, optional = false;
        If true, micrOMEGAs will compute the p-values for the Xenon1T,
        CRESST-III, PICO and DarkSide direct-detection experiments are
        computed. Default is false.
    z_invisible: bool, optional = false;
        If true, micrOMEGAs will determine if the invisible width of the Z boson
        is larger that 0.5 MeV. Default is false.
    lsp_nlsp_lep: bool, optional = false;
        If true, micrOMEGAs will check compatibility with the upper limit of
        σ(e⁺ + e⁻ → χ01 + χ0i) (where χ0i decays to χ01 + q + q) is violated
        by LEP. Default is false.
    z_prime_limits: bool, optional = false;
        If true, micrOMEGAs will check if model is excluded by Z' experiments.
        Default is false.
    monojet: bool, optional
        If true, micrOMEGAs will compute confidence-limit one
        σ(proton + proton → N1 + N2 + jet) from an 8 TeV CMS mono-jet analysis.
        Default is false.
    fast: bool, optional
        If true, micrOMEGAs will compute relic density using a approximation
        accurate to about 2% rather than use a more accurate (but much slower)
        method. Default is true.
    beps: float, optional
        The criteria for including a given coannihilation channel.
        The recommended value is 1e-4 <= beps <= 1e-6. beps=1 means only
        annihilation with the lightest odd particle are computed.
        Default is 1e-4.
    debug: bool, optional
        If true, debug/error messages will be printed. Default is False.
)pbdoc";

  settings->def(
      py::init<bool, bool, bool, bool, bool, bool, bool, bool, bool, bool, bool,
               bool, bool, bool, bool, bool, bool, double, bool>(),
      settings_init_doc.c_str(), py::arg("relic_density") = true,
      py::arg("masses") = true, py::arg("gmuon") = true, py::arg("bsg") = false,
      py::arg("bsmumu") = false, py::arg("btaunu") = false,
      py::arg("deltarho") = false, py::arg("rl23") = false,
      py::arg("d_taunu_and_munu") = false, py::arg("masslimits") = false,
      py::arg("nucleon_amplitudes") = false,
      py::arg("direct_detection_pvalues") = false,
      py::arg("z_invisible") = false, py::arg("lsp_nlsp_lep") = false,
      py::arg("z_prime_limits") = false, py::arg("monojet") = false,
      py::arg("fast") = true, py::arg("beps") = 1e-4, py::arg("debug") = false);

  // clang-format off
  settings->def_property("relic_density", &MicromegasSettings::get_relic_density, &MicromegasSettings::set_relic_density,settings_doc_relic_density.c_str());
  settings->def_property("masses", &MicromegasSettings::get_masses, &MicromegasSettings::set_masses, settings_doc_masses.c_str());
  settings->def_property("gmuon;", &MicromegasSettings::get_gmuon, &MicromegasSettings::set_gmuon,settings_doc_gmuon.c_str());
  settings->def_property("bsg;", &MicromegasSettings::get_bsg, &MicromegasSettings::set_bsg,settings_doc_bsg.c_str());
  settings->def_property("bsmumu", &MicromegasSettings::get_bsmumu,&MicromegasSettings::set_bsmumu,settings_doc_bsmumu.c_str());
  settings->def_property("btaunu", &MicromegasSettings::get_btaunu,&MicromegasSettings::set_btaunu,settings_doc_btaunu.c_str());
  settings->def_property("deltarho", &MicromegasSettings::get_deltarho, &MicromegasSettings::set_deltarho,settings_doc_deltarho.c_str());
  settings->def_property("rl23", &MicromegasSettings::get_rl23, &MicromegasSettings::set_rl23,settings_doc_rl23.c_str());
  settings->def_property("d_taunu_and_munu", &MicromegasSettings::get_d_taunu_and_munu, &MicromegasSettings::set_d_taunu_and_munu,settings_doc_d_taunu_and_munu.c_str());
  settings->def_property("masslimits", &MicromegasSettings::get_masslimits, &MicromegasSettings::set_masslimits,settings_doc_masslimits.c_str());
  settings->def_property("nucleon_amplitudes", &MicromegasSettings::get_nucleon_amplitudes, &MicromegasSettings::set_nucleon_amplitudes,settings_doc_nucleon_amplitudes.c_str());
  settings->def_property("direct_detection_pvalues", &MicromegasSettings::get_direct_detection_pvalues, &MicromegasSettings::set_direct_detection_pvalues,settings_doc_direct_detection_pvalues.c_str());
  settings->def_property("z_invisible", &MicromegasSettings::get_z_invisible, &MicromegasSettings::set_z_invisible,settings_doc_z_invisible.c_str());
  settings->def_property("lsp_nlsp_lep", &MicromegasSettings::get_lsp_nlsp_lep, &MicromegasSettings::set_lsp_nlsp_lep,settings_doc_lsp_nlsp_lep.c_str());
  settings->def_property("z_prime_limits", &MicromegasSettings::get_z_prime_limits, &MicromegasSettings::set_z_prime_limits,settings_doc_z_prime_limits.c_str());
  settings->def_property("monojet", &MicromegasSettings::get_monojet, &MicromegasSettings::set_monojet,settings_doc_monojet.c_str());
  settings->def_property("fast", &MicromegasSettings::get_fast, &MicromegasSettings::set_fast,settings_doc_fast.c_str());
  settings->def_property("beps", &MicromegasSettings::get_beps, &MicromegasSettings::set_beps,settings_doc_beps.c_str());
  settings->def_property("debug", &MicromegasSettings::get_debug, &MicromegasSettings::set_debug,settings_doc_debug.c_str());
  // clang-format on
}
