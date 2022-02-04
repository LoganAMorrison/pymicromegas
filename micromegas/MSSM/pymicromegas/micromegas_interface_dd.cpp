#include "micromegas.hpp"
#include "micromegas_interface.hpp"
#include <fmt/format.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void define_micromegas_direct_detection(py::module_ *ll) {

  using micromegas::DirectDetectionAmps;
  using micromegas::DirectDetectionResults;
  using micromegas::Micromegas;
  using DDExp = micromegas::DirectDetectionExperiment;

  // ========================================================================
  // ---- CDM1 --------------------------------------------------------------
  // ========================================================================

  py::enum_<DDExp>(*ll, "DirectDetectionExperiment", py::arithmetic(),
                   "Enumeration of direct detection experiments.")
      .value("Xenon1T2018", DDExp::Xenon1T2018, "XENON1T-2018 experiments.")
      .value("Cresst2019", DDExp::Cresst2019, "CRESST-III-2019 experiment.")
      .value("DarkSide2018", DDExp::DarkSide2018, "DarkSide-2018 experiment.")
      .value("Pico2019", DDExp::Pico2019, "PICO-2019 experiment.")
      .export_values();

  // ========================================================================
  // ---- CDM1 --------------------------------------------------------------
  // ========================================================================

  // clang-format off
  py::class_<DirectDetectionAmps>(*ll, "DirectDetectionAmps",
                   "Class holding the direct-detection amplitudes.")
      .def_readonly("proton_si", &DirectDetectionAmps::proton_si, "DM-Proton spin-independent amplitudes.")
      .def_readonly("proton_sd", &DirectDetectionAmps::proton_sd, "DM-Proton spin-dependent amplitudes.")
      .def_readonly("neutron_si", &DirectDetectionAmps::neutron_si, "DM-Proton spin-independent amplitudes.")
      .def_readonly("neutron_sd", &DirectDetectionAmps::neutron_sd, "DM-Proton spin-dependent amplitudes.")
      .def("__repr__", [](const DirectDetectionAmps& amps){
          const auto str1 = fmt::format("\tproton_si = [\u03C7 = {:.4e}, \u03C7\u203E = {:.4e}],\n", amps.proton_si[0], amps.proton_si[1]);
          const auto str2 = fmt::format("\tproton_sd = [\u03C7 = {:.4e}, \u03C7\u203E = {:.4e}],\n", amps.proton_sd[0], amps.proton_sd[1]);
          const auto str3 = fmt::format("\tneutron_si = [\u03C7 = {:.4e}, \u03C7\u203E = {:.4e}],\n", amps.neutron_si[0], amps.neutron_si[1]);
          const auto str4 = fmt::format("\tneutron_sd = [\u03C7 = {:.4e}, \u03C7\u203E = {:.4e}],\n", amps.neutron_sd[0], amps.neutron_sd[1]);

          return "DirectDetectionAmps(\n" + str1 + str2 + str3 + str4 + "\n)";
      });
  // clang-format on

  // ========================================================================
  // ---- CDM1 --------------------------------------------------------------
  // ========================================================================

  // clang-format off
  py::class_<DirectDetectionResults>(*ll, "DirectDetectionResults", 
                   "Class holding either scaling factors or p-values from a direct-detection exclusions.")
      .def_readonly("xenon1T", &DirectDetectionResults::xenon1T, "Results from XENON1T-2018 experiment.")
      .def_readonly("cresst", &DirectDetectionResults::cresst, "Results from CRESST-III-2019 experiment.")
      .def_readonly("pico", &DirectDetectionResults::pico, "Results from PICO-2019 experiment.")
      .def_readonly("darkside", &DirectDetectionResults::darkside, "Results from DarkSide 2018 experiment.")
      .def("__repr__", [](const DirectDetectionResults& results){
          const auto str1 = fmt::format("\tXENON1T (2018) = {:.4e},\n", results.xenon1T);
          const auto str2 = fmt::format("\tCRESST-III (2019) = {:.4e},\n", results.cresst);
          const auto str3 = fmt::format("\tDarkSide (2018) = {:.4e},\n", results.darkside);
          const auto str4 = fmt::format("\tPICO (2019) = {:.4e},\n", results.pico);

          return "DirectDetectionResults(\n" + str1 + str2 + str3 + str4 + "\n)";
      });
  // clang-format on

  // ========================================================================
  // ---- nucleonAmplitudes (CDM1) ------------------------------------------
  // ========================================================================

  static const char *nucleon_amplitudes_cdm1_doc = R"pbdoc(
  Calculates the amplitudes for CDM-nucleon elastic scattering at zero
  momentum for the 1st dark-matter candidate. 

  Notes
  -----
  The resulting object consists of four two-element arrays: `proton_si`, 
  `proton_sd`, `neutron_si`, and `neutron_sd`. The first element of each
  of these arrays corresponds to the DM cross-section while the second
  to the anti-DM cross-section. The `proton_*` arrays correspond to the
  scattering off protons, while the `neutron_*` arrays correspond to
  scattering off neutrons. The `*_si` arrays are the spin-independent
  parts while the `*_sd` are spin-dependent. The cross-sections (in Gev^-2) 
  are normalized such that the total cross-section is given by:
    σTOT = 4 mx^2 mn^2 / (pi * (mx + mn)^2) * (|Asi|^2 + 3 |Asd|^2)

  Returns
  -------
  cs: DirectDetectionCS 
    Object containing the proton and neutron spin-independent and
    spin-dependent cross sections in units of GeV^-2.

  Examples
  --------
  Compute the proton SI cross section:

  >>> import numpy as np
  >>> import micromegas.lowlevel as ll
  >>> ll->mssm_ewsb()
  >>> ll->sort_odd_particles()
  >>> amps = ll->nucleon_amplitudes_cdm1()
  >>> mx = ll->mcdm1()
  >>> nmass = 0.939 # nucleon mass 
  >>> pre = 4.0 / np.pi * (mx * nmass / (mx + nmass))**2
  >>> proton_si_cs = pre * amps.proton_si[0]**2

  Computing the other cross sections:
  
  >>> proton_sd_cs = pre * amps.proton_sd[0]**2 # Proton spin-dependent
  >>> neutron_si_cs = pre * amps.neutron_si[0]**2 # Neutron spin-independent
  >>> neutron_sd_cs = pre * amps.neutron_sd[0]**2 # Neutron spin-dependent
)pbdoc";

  ll->def(
      "nucleon_amplitudes_cdm1",
      []() { return Micromegas::nucleon_amplitudes_cdm1(); },
      nucleon_amplitudes_cdm1_doc);

  // ========================================================================
  // ---- nucleonAmplitudes (CDM2) ------------------------------------------
  // ========================================================================

  static const char *nucleon_amplitudes_cdm2_doc = R"pbdoc(
  Calculates the amplitudes for CDM-nucleon elastic scattering at zero
  momentum for the 2nd dark-matter candidate. 

  Notes
  -----
  The resulting object consists of four two-element arrays: `proton_si`, 
  `proton_sd`, `neutron_si`, and `neutron_sd`. The first element of each
  of these arrays corresponds to the DM cross-section while the second
  to the anti-DM cross-section. The `proton_*` arrays correspond to the
  scattering off protons, while the `neutron_*` arrays correspond to
  scattering off neutrons. The `*_si` arrays are the spin-independent
  parts while the `*_sd` are spin-dependent. The cross-sections (in Gev^-2) 
  are normalized such that the total cross-section is given by:
    σTOT = 4 mx^2 mn^2 / (pi * (mx + mn)^2) * (|Asi|^2 + 3 |Asd|^2)

  Returns
  -------
  cs: DirectDetectionCS 
    Object containing the proton and neutron spin-independent and
    spin-dependent cross sections in units of GeV^-2.

  Examples
  --------
  Compute the proton SI cross section:

  >>> import numpy as np
  >>> import micromegas.lowlevel as ll
  >>> ll->mssm_ewsb()
  >>> ll->sort_odd_particles()
  >>> amps = ll->nucleon_amplitudes_cdm2()
  >>> mx = ll->mcdm1()
  >>> nmass = 0.939 # nucleon mass 
  >>> pre = 4.0 / np.pi * (mx * nmass / (mx + nmass))**2
  >>> proton_si_cs = pre * amps.proton_si[0]**2

  Computing the other cross sections:
  
  >>> proton_sd_cs = pre * amps.proton_sd[0]**2 # Proton spin-dependent
  >>> neutron_si_cs = pre * amps.neutron_si[0]**2 # Neutron spin-independent
  >>> neutron_sd_cs = pre * amps.neutron_sd[0]**2 # Neutron spin-dependent
)pbdoc";

  ll->def(
      "nucleon_amplitudes_cdm2",
      []() { return Micromegas::nucleon_amplitudes_cdm2(); },
      nucleon_amplitudes_cdm2_doc);

  // ========================================================================
  // ---- DD_factorCS (Maxwell) ---------------------------------------------
  // ========================================================================

  //   static const char *direct_detection_factor_cs_maxwell_doc = R"pbdoc(
  //   Returns the overall factor which should be applied to the cross sections,
  //   σSIP , σSIN , σSDP , σSDN to reach the exclusion level α. Assumes a
  //   Maxwell velocity distribution for DM.

  //   Parameters
  //   ----------
  //   pval:
  //     p-value for exclusion level.
  //   proton_si_cs:
  //     Spin-independent DM-proton cross section.
  //   proton_sd_cs:
  //     Spin-dependent DM-proton cross section.
  //   neutron_si_cs:
  //     Spin-independent DM-neutron cross section.
  //   neutron_sd_cs:
  //     Spin-dependent DM-neutron cross section.

  //   Returns
  //   -------
  //   factor: float
  //     Factor needed to rescale cross-section to reach exclusion level.

  //   See Also
  //   --------
  //   micromegas.lowlevel.direct_detection_factor_cs_shmpp:
  //     Same but uses the SHM++ velocity distribution.
  //   micromegas.lowlevel.nucleon_amplitudes_cdm1
  //   micromegas.lowlevel.nucleon_amplitudes_cdm2:
  //     Computing amplitudes needed for computing cross-sections.
  // )pbdoc";

  //   ll->def(
  //       "direct_detection_factor_cs_maxwell",
  //       [](double pval, double proton_si_cs, double proton_sd_cs,
  //          double neutron_si_cs, double neutron_sd_cs) {
  //         if (Micromegas::get_cdm1() == "null") {
  //           py::print(py::str("Dark matter not yet determined. Make sure "
  //                             "'sort_odd_particles' has been called."));
  //         }
  //         return Micromegas::direct_detection_factor_cs_maxwell(
  //             pval, proton_si_cs, proton_sd_cs, neutron_si_cs,
  //             neutron_sd_cs);
  //       },
  //       direct_detection_factor_cs_maxwell_doc, py::arg("pval"),
  //       py::arg("proton_si_cs"), py::arg("proton_sd_cs"),
  //       py::arg("neutron_si_cs"), py::arg("neutron_sd_cs"));

  // ========================================================================
  // ---- DD_factorCS (SHM++) -----------------------------------------------
  // ========================================================================

  //   static const char *direct_detection_factor_cs_shmpp_doc = R"pbdoc(
  //   Returns the overall factor which should be applied to the cross sections,
  //   σSIP , σSIN , σSDP , σSDN to reach the exclusion level α. Assumes a SHM++
  //   velocity distribution for DM.

  //   Parameters
  //   ----------
  //   pval:
  //     p-value for exclusion level.
  //   proton_si_cs:
  //     Spin-independent DM-proton cross section.
  //   proton_sd_cs:
  //     Spin-dependent DM-proton cross section.
  //   neutron_si_cs:
  //     Spin-independent DM-neutron cross section.
  //   neutron_sd_cs:
  //     Spin-dependent DM-neutron cross section.

  //   Returns
  //   -------
  //   factor: float
  //     Factor needed to rescale cross-section to reach exclusion level.

  //   See Also
  //   --------
  //   micromegas.lowlevel.direct_detection_factor_cs_shmpp:
  //     Same but uses the Maxwell velocity distribution.
  //   micromegas.lowlevel.nucleon_amplitudes_cdm1
  //   micromegas.lowlevel.nucleon_amplitudes_cdm2:
  //     Computing amplitudes needed for computing cross-sections.
  // )pbdoc";

  //   ll->def(
  //       "direct_detection_factor_cs_shmpp",
  //       [](double pval, double proton_si_cs, double proton_sd_cs,
  //          double neutron_si_cs, double neutron_sd_cs) {
  //         if (Micromegas::get_cdm1() == "null") {
  //           py::print(py::str("Dark matter not yet determined. Make sure "
  //                             "'sort_odd_particles' has been called."));
  //         }
  //         return Micromegas::direct_detection_factor_cs_shmpp(
  //             pval, proton_si_cs, proton_sd_cs, neutron_si_cs,
  //             neutron_sd_cs);
  //       },
  //       direct_detection_factor_cs_shmpp_doc, py::arg("pval"),
  //       py::arg("proton_si_cs"), py::arg("proton_sd_cs"),
  //       py::arg("neutron_si_cs"), py::arg("neutron_sd_cs"));

  // ========================================================================
  // ---- DD_pvalCS (Maxwell) -----------------------------------------------
  // ========================================================================

  //   static const char *direct_detection_pval_cs_maxwell_doc = R"pbdoc(
  //   Calculate α = 1 − CL for a model with the given cross sections. Assumes a
  //   Maxwell velocity distribution for DM.

  //   Parameters
  //   ----------
  //   proton_si_cs:
  //     Spin-independent DM-proton cross section.
  //   proton_sd_cs:
  //     Spin-dependent DM-proton cross section.
  //   neutron_si_cs:
  //     Spin-independent DM-neutron cross section.
  //   neutron_sd_cs:
  //     Spin-dependent DM-neutron cross section.

  //   Returns
  //   -------
  //   factor: float
  //     Factor needed to rescale cross-section to reach exclusion level.

  //   See Also
  //   --------
  //   micromegas.lowlevel.nucleon_amplitudes_cdm1
  //   micromegas.lowlevel.nucleon_amplitudes_cdm2:
  //     Computing amplitudes needed for computing cross-sections.
  // )pbdoc";

  //   ll->def(
  //       "direct_detection_pval_cs_maxwell",
  //       [](double proton_si_cs, double proton_sd_cs, double neutron_si_cs,
  //          double neutron_sd_cs) {
  //         if (Micromegas::get_cdm1() == "null") {
  //           py::print(py::str("Dark matter not yet determined. Make sure "
  //                             "'sort_odd_particles' has been called."));
  //         }
  //         return Micromegas::direct_detection_pval_cs_maxwell(
  //             proton_si_cs, proton_sd_cs, neutron_si_cs, neutron_sd_cs);
  //       },
  //       direct_detection_pval_cs_maxwell_doc, py::arg("proton_si_cs"),
  //       py::arg("proton_sd_cs"), py::arg("neutron_si_cs"),
  //       py::arg("neutron_sd_cs"));

  // ========================================================================
  // ---- DD_factorCS (SHM++) -----------------------------------------------
  // ========================================================================

  //   static const char *direct_detection_pval_cs_shmpp_doc = R"pbdoc(
  //   Calculate α = 1 − CL for a model with the given cross sections. Assumes a
  //   SHM++ velocity distribution for DM.

  //   Parameters
  //   ----------
  //   proton_si_cs:
  //     Spin-independent DM-proton cross section.
  //   proton_sd_cs:
  //     Spin-dependent DM-proton cross section.
  //   neutron_si_cs:
  //     Spin-independent DM-neutron cross section.
  //   neutron_sd_cs:
  //     Spin-dependent DM-neutron cross section.

  //   Returns
  //   -------
  //   factor: float
  //     Factor needed to rescale cross-section to reach exclusion level.

  //   See Also
  //   --------
  //   micromegas.lowlevel.nucleon_amplitudes_cdm1
  //   micromegas.lowlevel.nucleon_amplitudes_cdm2:
  //     Computing amplitudes needed for computing cross-sections.
  // )pbdoc";

  //   ll->def(
  //       "direct_detection_pval_cs_shmpp",
  //       [](double proton_si_cs, double proton_sd_cs, double neutron_si_cs,
  //          double neutron_sd_cs) {
  //         if (Micromegas::get_cdm1() == "null") {
  //           py::print(py::str("Dark matter not yet determined. Make sure "
  //                             "'sort_odd_particles' has been called."));
  //         }
  //         return Micromegas::direct_detection_pval_cs_shmpp(
  //             proton_si_cs, proton_sd_cs, neutron_si_cs, neutron_sd_cs);
  //       },
  //       direct_detection_pval_cs_shmpp_doc, py::arg("proton_si_cs"),
  //       py::arg("proton_sd_cs"), py::arg("neutron_si_cs"),
  //       py::arg("neutron_sd_cs"));

  // ========================================================================
  // ---- DD_factor (Maxwell) -----------------------------------------------
  // ========================================================================

  static const char *direct_detection_factor_maxwell_doc = R"pbdoc(
  Returns the overall factor to reach the exclusion level α which should be
  applied to the cross sections calculated by micrOMEGAs using DM model under
  consideration. Assumes a Maxwell velocity distribution for DM.

  Parameters
  ----------
  pval:
    p-value for exclusion level.

  Returns
  -------
  factor: float 
    Factor needed to rescale cross-section to reach exclusion level.
)pbdoc";

  ll->def(
      "direct_detection_factor_maxwell",
      [](double pval) {
        if (Micromegas::get_cdm1() == "null") {
          py::print(py::str("Dark matter not yet determined. Make sure "
                            "'sort_odd_particles' has been called."));
        }
        return Micromegas::direct_detection_factor_maxwell(pval);
      },
      direct_detection_factor_maxwell_doc);

  // ========================================================================
  // ---- DD_factor (SHM++) -------------------------------------------------
  // ========================================================================

  static const char *direct_detection_factor_shmpp_doc = R"pbdoc(
  Returns the overall factor to reach the exclusion level α which should be
  applied to the cross sections calculated by micrOMEGAs using DM model under
  consideration. Assumes a SHM++ velocity distribution for DM.

  Parameters
  ----------
  pval:
    p-value for exclusion level.

  Returns
  -------
  factor: float 
    Factor needed to rescale cross-section to reach exclusion level.
)pbdoc";

  ll->def(
      "direct_detection_factor_shmpp",
      [](double pval) {
        if (Micromegas::get_cdm1() == "null") {
          py::print(py::str("Dark matter not yet determined. Make sure "
                            "'sort_odd_particles' has been called."));
        }
        return Micromegas::direct_detection_factor_shmpp(pval);
      },
      direct_detection_factor_shmpp_doc);

  // ========================================================================
  // ---- DD_pval (Maxwell) -------------------------------------------------
  // ========================================================================

  static const char *direct_detection_pval_maxwell_doc = R"pbdoc(
  Calculates the value α = 1 − C.L. using cross sections calculated by
  micrOMEGAs using DM model under consideration. For exmaple, a return 
  value of 0.1 corresponds to a 90% exclusion. Assumes a Maxwell velocity
  distribution for DM.

  Returns
  -------
  alpha: float 
    Value of 1 - C.L.
)pbdoc";

  ll->def(
      "direct_detection_pval_maxwell",
      []() {
        if (Micromegas::get_cdm1() == "null") {
          py::print(py::str("Dark matter not yet determined. Make sure "
                            "'sort_odd_particles' has been called."));
        }
        return Micromegas::direct_detection_pval_maxwell();
      },
      direct_detection_pval_maxwell_doc);

  // ========================================================================
  // ---- DD_pval (SHM++) ---------------------------------------------------
  // ========================================================================

  static const char *direct_detection_pval_shmpp_doc = R"pbdoc(
  Calculates the value α = 1 − C.L. using cross sections calculated by
  micrOMEGAs using DM model under consideration. For exmaple, a return 
  value of 0.1 corresponds to a 90% exclusion. Assumes a SHM++ velocity
  distribution for DM.

  Returns
  -------
  alpha: float 
    Value of 1 - C.L.
)pbdoc";

  ll->def(
      "direct_detection_pval_shmpp",
      []() {
        if (Micromegas::get_cdm1() == "null") {
          py::print(py::str("Dark matter not yet determined. Make sure "
                            "'sort_odd_particles' has been called."));
        }
        return Micromegas::direct_detection_pval_maxwell();
      },
      direct_detection_pval_shmpp_doc);

  // ========================================================================
  // ---- XENON1T_90 --------------------------------------------------------
  // ========================================================================

  static const char *xenon1T_90_doc = R"pbdoc(
  Returns the excluded 90% SI cross section in cm^2 from the XENON1T
  experiment. For a DM mass outside 6 < MDM < 1000 GeV, NaN is returned.

  Parameters
  ----------
  mass: float
    Dark matter mass.

  Returns
  -------
  cs: float
    Constraint on cross-section in cm^2.
)pbdoc";

  ll->def(
      "xenon1T_90", [](double mass) { return Micromegas::xenon1T_90(mass); },
      xenon1T_90_doc, py::arg("mass"));

  // ========================================================================
  // ---- XENON1T_SDp_90 ----------------------------------------------------
  // ========================================================================

  static const char *xenon1T_sdp_90_doc = R"pbdoc(
  Returns the excluded 90% SD proton cross section in cm^2 from the XENON1T
  experiment. For a DM mass outside 6 < MDM < 1000 GeV, NaN is returned.

  Parameters
  ----------
  mass: float
    Dark matter mass.

  Returns
  -------
  cs: float
    Constraint on cross-section in cm^2.
)pbdoc";

  ll->def(
      "xenon1T_sdp_90",
      [](double mass) { return Micromegas::xenon1T_sdp_90(mass); },
      xenon1T_sdp_90_doc, py::arg("mass"));

  // ========================================================================
  // ---- XENON1T_SDn_90 ----------------------------------------------------
  // ========================================================================

  static const char *xenon1T_sdn_90_doc = R"pbdoc(
  Returns the excluded 90% SD neutron cross section in cm^2 from the XENON1T
  experiment. For a DM mass outside 6 < MDM < 1000 GeV, NaN is returned.

  Parameters
  ----------
  mass: float
    Dark matter mass.

  Returns
  -------
  cs: float
    Constraint on cross-section in cm^2.
)pbdoc";

  ll->def(
      "xenon1T_sdn_90",
      [](double mass) { return Micromegas::xenon1T_sdn_90(mass); },
      xenon1T_sdn_90_doc, py::arg("mass"));

  // ========================================================================
  // ---- DS50_90 -----------------------------------------------------------
  // ========================================================================

  static const char *darkside50_90_doc = R"pbdoc(
  Returns the excluded 90% SI cross section in cm^2 from the DarkSide
  experiment. For a DM mass outside 0.7 < MDM < 15 GeV, NaN is returned.

  Parameters
  ----------
  mass: float
    Dark matter mass.

  Returns
  -------
  cs: float
    Constraint on cross-section in cm^2.
)pbdoc";

  ll->def(
      "darkside50_90",
      [](double mass) { return Micromegas::darkside50_90(mass); },
      darkside50_90_doc, py::arg("mass"));

  // ========================================================================
  // ---- DS50_90_noB -------------------------------------------------------
  // ========================================================================

  ll->def("darkside50_90_nob",
          [](double mass) { return Micromegas::darkside50_90_nob(mass); });

  // ========================================================================
  // ---- CRESST_III_90 -----------------------------------------------------
  // ========================================================================

  static const char *cresst3_90_doc = R"pbdoc(
  Returns the excluded 90% SI cross section in cm^2 from the CRESST III
  experiment. For a DM mass outside 0.35 < MDM < 12 GeV, NaN is returned.

  Parameters
  ----------
  mass: float
    Dark matter mass.

  Returns
  -------
  cs: float
    Constraint on cross-section in cm^2.
)pbdoc";

  ll->def(
      "cresst3_90", [](double mass) { return Micromegas::cresst3_90(mass); },
      cresst3_90_doc, py::arg("mass"));

  // ========================================================================
  // ---- CRESST_III_SDn_90 -------------------------------------------------
  // ========================================================================

  static const char *cresst3_sdn_90_doc = R"pbdoc(
  Returns the excluded 90% SD neutron cross section in cm^2 from the CRESST
  III experiment. For a DM mass outside 0.35 < MDM < 12 GeV, NaN is returned.

  Parameters
  ----------
  mass: float
    Dark matter mass.

  Returns
  -------
  cs: float
    Constraint on cross-section in cm^2.
)pbdoc";

  ll->def(
      "cresst3_sdn_90",
      [](double mass) { return Micromegas::cresst3_sdn_90(mass); },
      cresst3_sdn_90_doc, py::arg("mass"));

  // ========================================================================
  // ---- PICO60_90 ---------------------------------------------------------
  // ========================================================================

  static const char *pico60_90_doc = R"pbdoc(
  Returns the excluded 90% SI cross section in cm^2 from the PICO-60
  experiment. For a DM mass outside 3 < MDM < 10000 GeV, NaN is returned.

  Parameters
  ----------
  mass: float
    Dark matter mass.

  Returns
  -------
  cs: float
    Constraint on cross-section in cm^2.
)pbdoc";

  ll->def(
      "pico60_90", [](double mass) { return Micromegas::pico60_90(mass); },
      pico60_90_doc, py::arg("mass"));

  // ========================================================================
  // ---- PICO60_SDp_90 -----------------------------------------------------
  // ========================================================================

  static const char *pico60_sdp_90_doc = R"pbdoc(
  Returns the excluded 90% SD proton cross section in cm^2 from the PICO-60
  experiment. For a DM mass outside 3 < MDM < 10000 GeV, NaN is returned.

  Parameters
  ----------
  mass: float
    Dark matter mass.

  Returns
  -------
  cs: float
    Constraint on cross-section in cm^2.
)pbdoc";

  ll->def(
      "pico60_sdp_90",
      [](double mass) { return Micromegas::pico60_sdp_90(mass); },
      pico60_sdp_90_doc, py::arg("mass"));
}
