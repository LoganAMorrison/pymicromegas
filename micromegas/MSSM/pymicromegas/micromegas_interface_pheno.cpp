
#include "micromegas.hpp"
#include "micromegas_interface.hpp"
#include <fmt/format.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void define_micromegas_pheno(py::module_ *ll) {

  using micromegas::Micromegas;
  using micromegas::MssmParameter;

  // -- gmuon ------

  ll->def(
      "mssm_gmuon", []() { return Micromegas::mssm_gmuon(); },
      R"pbdoc(
  Returns the value of the supersymmetric contribution to the anomalous
  magnetic moment of the muon.
)pbdoc");

  // -- deltarho ------

  ll->def(
      "mssm_deltarho", []() { return Micromegas::mssm_deltarho(); },
      R"pbdoc(
  Calculates the ∆ρ parameter in the MSSM. It contains for example the
  stop/sbottom contributions, as well as the two-loop QCD corrections due to
  gluon exchange and the correction due to gluino exchange in the heavy
  gluino limit
)pbdoc");

  // -- bsgnlo ------

  ll->def(
      "mssm_bsgnlo", []() { return Micromegas::mssm_bsgnlo(); },
      R"pbdoc(
  Returns the value of the branching ratio for b → sγ, see Appendix A. We
  have included some new contributions beyond the leading order that are
  especially important for high tan β. SMbsg gives the SM contribution.
)pbdoc");

  // -- bsmumu ------

  ll->def(
      "mssm_bsmumu", []() { return Micromegas::mssm_bsmumu(); },
      R"pbdoc(
  Returns the value of the branching ratio Bs → μ+μ− in the MSSM. It includes
  the loop contributions due to chargino, sneutrino, stop and Higgs exchange.
  The ∆mb effect relevant for high tan β is taken into account.
)pbdoc");

  // -- btaunu ------

  ll->def(
      "mssm_btaunu", []() { return Micromegas::mssm_btaunu(); },
      R"pbdoc(
  Computes the ratio between the MSSM and SM branching fractions for  ̄B+ → τ
  +ντ .
)pbdoc");

  // -- Rl23 ------

  ll->def(
      "mssm_rl23", []() { return Micromegas::mssm_rl23(); },
      R"pbdoc(
  Computes the ratio of the MSSM to SM value for Rl23 in K+ → μν due to a
  charged higgs contribution.
)pbdoc");

  // -- dtaunu ------

  ll->def(
      "mssm_d_taunu_and_munu",
      []() { return Micromegas::mssm_d_taunu_and_munu(); },
      R"pbdoc(
  Computes the branching ratio for D+s → τ +ντ . `dmunu` gives the branching
  ratio for D+s →μ+νμ.
)pbdoc");

  // -- masslimits ------

  ll->def(
      "mssm_masslimits", []() { return Micromegas::mssm_masslimits(); },
      R"pbdoc(
  Returns a positive value and prints a WARNING when the choice of parameters
  conflicts with a direct accelerator limits on sparticle masses from LEP.
  The constraint on the light Higgs mass from the LHC is included.
)pbdoc");

  ll->def(
      "relic_density",
      [](bool fast = true, double beps = 1e-4) {
        return Micromegas::relic_density(fast, beps);
      },
      R"pbdoc(
  Computes the dark matter relic density Ωh² and the scaled dark matter
  freeze-out temperature: xᶠ = mᵡ/Tᶠ.

  Parameters
  ----------
  fast: bool
    If true, an approximate scheme is used which is accurate to about 2%.
  beps:
    Parameter defining which coannihilation channels are included. Recommended
    values are between 1e-4 and 1e-6.

  Returns
  -------
  rd: float
    Dark matter relic density Ωh².
  xf: float
    Scaled dark matter freeze-out temperature: xᶠ = mᵡ/Tᶠ.
)pbdoc",
      py::arg("fast") = true, py::arg("beps") = 1e-4);

  // -- Zinvisible ------

  ll->def(
      "z_invisible", []() { return Micromegas::z_invisible(); },
      R"pbdoc(
  Returns 1 and prints a WARNING if the invisible width of the Z boson of the
  Standard Model is larger than 0.5 MeV and returns 0 if this
  constraint is fulfilled. This routine can be used in any model with one DM
  where the Z boson is uniquely defined by its PDG=23 and whether the neutral
  LSP is its own antiparticle or not.
)pbdoc");

  // -- LspNlsp_LEP ------

  ll->def(
      "lsp_nlsp_lep", []() { return Micromegas::lsp_nlsp_lep(); },
      R"pbdoc(
  Returns 0 if the scenario considered is not excluded by Z′ constraints, 1
  if the point is excluded and 2 if both subroutines dealing with Z′
  constraints cannot test the given scenario. The routine can be used for any
  Z′ uniquely defined by the PDG code 32. Currently two types of searches
  defined in different subroutines of Zprimelimits() are implemented.
  The 3.2/fb Z′ search in the dilepton final state at √s = 13 TeV from ATLAS
  is considered in the first subroutine Zprime_dilepton. If the scenario
  considered is allowed or not tested by Zprime_dilepton, a second subroutine
  called Zprime_dijet analyses the point using constraints from LHC dijet
  searches at √s = 8 TeV and at √s = 13 TeV. This subroutine uses the
  recasting performed in [88] for a combination of ATLAS and CMS searches.
  Zprimelimits() returns 1 if MZ′ < 0.5 TeV and 2 for points for which the
  narrow-width approximation is not valid, i.e. Γ/MZ′ > 0.3.
)pbdoc");

  // -- Zprimelimits ------

  ll->def(
      "z_prime_limits", []() { return Micromegas::z_prime_limits(); },
      R"pbdoc(
  Checks the compatibility with the upper limit on the cross section for
  the production of neutralinos σ(e+ e− →  ̃χ01  ̃χ0i ), i != 1, when the
  heavier neutralino decays into quark pairs and the LSP,  ̃χ0i →  ̃χ01 + q +
  q. The function returns σ × BR = ∑ i σ(e+e− → ̃χ01  ̃χ0i) × BR(χ0i →  ̃χ01 + q
  + ̄q) in pb as well as a flag greater than one if σ × BR > 0.1(0.5) pb if
  mNLSP > (<)100 GeV. This function can also be applied for non-SUSY models
  which feature the same signature, in this case the function will compute
  the cross section for production of the LSP and any other neutral particle
  from the odd sector which can decay into the LSP and a Z boson.
)pbdoc");

  // -- monoJet ------

  ll->def(
      "monojet", []() { return Micromegas::monojet(); },
      R"pbdoc(
  Computes the cross section for p, p → pName1, pName2 + jet at √s = 8 TeV
  where pName1, pName2 are the names of neutral outgoing particles and jet
  includes light quarks (u,d,s) and gluons. The function returns the
  resulting confidence level obtained with the CLs technique for each signal
  region of the 8 TeV CMS mono-jet analysis and chooses the most
  constraining one.
)pbdoc");
}
