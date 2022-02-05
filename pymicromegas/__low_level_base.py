from abc import ABCMeta, abstractmethod
from typing import Tuple

from . import micromegas as low_level
from .micromegas import DirectDetectionAmps, DirectDetectionResults


class LowLevelMicromegasBase:

    __metaclass__ = ABCMeta

    def __init__(self) -> None:
        pass

    @abstractmethod
    def _initialize(self):
        pass

    def nucleon_amplitudes_cdm1(self) -> DirectDetectionAmps:
        """
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
        parts while the `*_sd` are spin-dependent. The cross-sections
        (in Gev^-2) are normalized such that the total cross-section is given
        by:
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

        # Proton spin-dependent
        >>> proton_sd_cs = pre * amps.proton_sd[0]**2
        # Neutron spin-independent
        >>> neutron_si_cs = pre * amps.neutron_si[0]**2
        # Neutron spin-dependent
        >>> neutron_sd_cs = pre * amps.neutron_sd[0]**2
        """
        self._initialize()
        return low_level.nucleon_amplitudes_cdm1()

    def nucleon_amplitudes_cdm2(self) -> DirectDetectionAmps:
        """
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
        parts while the `*_sd` are spin-dependent. The cross-sections
        (in Gev^-2) are normalized such that the total cross-section is given
        by:
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

        # Proton spin-dependent
        >>> proton_sd_cs = pre * amps.proton_sd[0]**2
        # Neutron spin-independent
        >>> neutron_si_cs = pre * amps.neutron_si[0]**2
        # Neutron spin-dependent
        >>> neutron_sd_cs = pre * amps.neutron_sd[0]**2
        """
        self._initialize()
        return low_level.nucleon_amplitudes_cdm2()

    def direct_detection_factor(
        self, pval: float, maxwell: bool = True
    ) -> DirectDetectionResults:
        """
        Returns the overall factor to reach the exclusion level α which should
        be applied to the cross sections calculated by micrOMEGAs using DM
        model under consideration. Assumes a Maxwell velocity distribution for
        DM.

        Parameters
        ----------
        pval:
          p-value for exclusion level.

        Returns
        -------
        factor: float
            Factor needed to rescale cross-section to reach exclusion level.
        """
        self._initialize()
        if maxwell:
            return low_level.direct_detection_factor_maxwell(pval)
        return low_level.direct_detection_factor_shmpp(pval)

    def direct_detection_pval(self, maxwell: bool = True) -> DirectDetectionResults:
        """
        Calculates the value α = 1 − C.L. using cross sections calculated by
        micrOMEGAs using DM model under consideration. For exmaple, a return
        value of 0.1 corresponds to a 90% exclusion. Assumes a Maxwell velocity
        distribution for DM.

        Returns
        -------
        alpha: float
            Value of 1 - C.L.
        """
        self._initialize()
        if maxwell:
            return low_level.direct_detection_pval_maxwell()
        return low_level.direct_detection_pval_shmpp()

    def xenon1T_90(self, mass: float) -> float:
        """
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
        """
        return low_level.xenon1T_90(mass)

    def xenon1T_sdp_90(self, mass: float) -> float:
        """
        Returns the excluded 90% SD proton cross section in cm^2 from the
        XENON1T experiment. For a DM mass outside 6 < MDM < 1000 GeV, NaN is
        returned.

        Parameters
        ----------
        mass: float
          Dark matter mass.

        Returns
        -------
        cs: float
            Constraint on cross-section in cm^2.
        """
        return low_level.xenon1T_sdp_90(mass)

    def xenon1T_sdn_90(self, mass: float) -> float:
        """
        Returns the excluded 90% SD neutron cross section in cm^2 from the
        XENON1T experiment. For a DM mass outside 6 < MDM < 1000 GeV, NaN is
        returned.

        Parameters
        ----------
        mass: float
            Dark matter mass.

        Returns
        -------
        cs: float
            Constraint on cross-section in cm^2.
        """
        return low_level.xenon1T_sdn_90(mass)

    def darkside50_90(self, mass: float) -> float:
        """
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
        """
        return low_level.darkside50_90(mass)

    def darkside50_90_nob(self, mass: float) -> float:
        return low_level.darkside50_90_nob(mass)

    def cresst3_90(self, mass: float) -> float:
        """
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
        """
        return low_level.cresst3_90(mass)

    def cresst3_sdn_90(self, mass: float) -> float:
        """
        Returns the excluded 90% SD neutron cross section in cm^2 from the
        CRESST III experiment. For a DM mass outside 0.35 < MDM < 12 GeV, NaN
        is returned.

        Parameters
        ----------
        mass: float
            Dark matter mass.

        Returns
        -------
        cs: float
            Constraint on cross-section in cm^2.
        """
        return low_level.cresst3_sdn_90(mass)

    def pico60_90(self, mass: float) -> float:
        """
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
        """
        return low_level.pico60_90(mass)

    def pico60_sdp_90(self, mass: float) -> float:
        """
        Returns the excluded 90% SD proton cross section in cm^2 from the
        PICO-60 experiment. For a DM mass outside 3 < MDM < 10000 GeV, NaN is
        returned.

        Parameters
        ----------
        mass: float
            Dark matter mass.

        Returns
        -------
        cs: float
            Constraint on cross-section in cm^2.
        """
        return low_level.pico60_sdp_90(mass)

    def cdm1(self) -> str:
        """
        Name of the 1st DM particle.
        """
        self._initialize()
        return low_level.cdm1()

    def cdm2(self) -> str:
        """
        Name of the 2nd DM particle.
        """
        self._initialize()
        return low_level.cdm2()

    def mcdm1(self) -> float:
        """
        Mass of the 1st DM particle.
        """
        self._initialize()
        return low_level.mcdm1()

    def mcdm2(self) -> float:
        """
        Mass of the 2nd DM particle.
        """
        self._initialize()
        return low_level.mcdm2()

    def gmuon(self) -> float:
        """
        Returns the value of the supersymmetric contribution to the anomalous
        magnetic moment of the muon.
        """
        self._initialize()
        return low_level.mssm_gmuon()

    def deltarho(self) -> float:
        """
        Calculates the ∆ρ parameter in the MSSM. It contains for example the
        stop/sbottom contributions, as well as the two-loop QCD corrections
        due to gluon exchange and the correction due to gluino exchange in the
        heavy gluino limit
        """
        self._initialize()
        return low_level.mssm_deltarho()

    def bsgnlo(self) -> Tuple[float, float]:
        """
        Returns the value of the branching ratio for b → sγ, see Appendix A. We
        have included some new contributions beyond the leading order that are
        especially important for high tan β. SMbsg gives the SM contribution.
        """
        self._initialize()
        return low_level.mssm_bsgnlo()

    def bsmumu(self) -> float:
        """
        Returns the value of the branching ratio Bs → μ+μ− in the MSSM.
        It includes the loop contributions due to chargino, sneutrino, stop and
        Higgs exchange. The ∆mb effect relevant for high tan β is taken into
        account.
        """
        self._initialize()
        return low_level.mssm_bsmumu()

    def btaunu(self) -> float:
        """
        Computes the ratio between the MSSM and SM branching fractions for
        B+ → τ+ντ.
        """
        self._initialize()
        return low_level.mssm_btaunu()

    def rl23(self) -> float:
        """
        Computes the ratio of the MSSM to SM value for Rl23 in K+ → μν due to a
        charged higgs contribution.
        """
        self._initialize()
        return low_level.mssm_rl23()

    def d_taunu_and_munu(self) -> Tuple[float, float]:
        """
        Computes the branching ratio for D+s → τ +ντ . `dmunu` gives the
        branching ratio for D+s →μ+νμ.
        """
        self._initialize()
        return low_level.mssm_d_taunu_and_munu()

    def masslimits(self) -> float:
        """
        Returns a positive value and prints a WARNING when the choice of
        parameters conflicts with a direct accelerator limits on sparticle
        masses from LEP. The constraint on the light Higgs mass from the LHC
        is included.
        """
        self._initialize()
        return low_level.mssm_masslimits()

    def relic_density(
        self, fast: bool = True, beps: float = 1e-4
    ) -> Tuple[float, float]:
        """
        Computes the dark matter relic density Ωh² and the scaled dark matter
        freeze-out temperature: xᶠ = mᵡ/Tᶠ.

        Parameters
        ----------
        fast: bool
          If true, an approximate scheme is used which is accurate to about 2%.
        beps:
          Parameter defining which coannihilation channels are included.
          Recommended values are between 1e-4 and 1e-6.

        Returns
        -------
        rd: float
          Dark matter relic density Ωh².
        xf: float
          Scaled dark matter freeze-out temperature: xᶠ = mᵡ/Tᶠ.
        """
        self._initialize()
        return low_level.relic_density(fast, beps)

    def z_invisible(self) -> bool:
        """
        Returns 1 and prints a WARNING if the invisible width of the Z boson
        of the Standard Model is larger than 0.5 MeV and returns 0 if this
        constraint is fulfilled. This routine can be used in any model with
        one DM where the Z boson is uniquely defined by its PDG=23 and whether
        the neutral LSP is its own antiparticle or not.
        """
        self._initialize()
        return low_level.z_invisible()

    def lsp_nlsp_lep(self) -> Tuple[float, float]:
        """
        Returns 0 if the scenario considered is not excluded by Z′
        constraints, 1 if the point is excluded and 2 if both subroutines
        dealing with Z′ constraints cannot test the given scenario.
        """
        self._initialize()
        return low_level.lsp_nlsp_lep()

    def z_prime_limits(self) -> bool:
        """
        Checks the compatibility with the upper limit on the cross section for
        the production of neutralinos σ(e+ e− →  ̃χ01  ̃χ0i ), i != 1, when the
        heavier neutralino decays into quark pairs and the LSP,  ̃
        χ0i →  ̃χ01 + q + q. The function returns
            σ × BR = ∑ i σ(e+e− → ̃χ01  ̃χ0i) × BR(χ0i →  ̃χ01 + q + ̄q)
        in pb as well as a flag greater than one if σ × BR > 0.1(0.5)
        pb if mNLSP > (<)100 GeV. This function can also be applied for
        non-SUSY models which feature the same signature, in this case the
        function will compute the cross section for production of the LSP and
        any other neutral particle from the odd sector which can decay into the
        LSP and a Z boson.
        """
        self._initialize()
        return low_level.z_prime_limits()

    def monojet(self) -> float:
        """
        Computes the cross section for p, p → pName1, pName2 + jet at √s =
        8 TeV where pName1, pName2 are the names of neutral outgoing particles
        and jet includes light quarks (u,d,s) and gluons. The function returns
        the resulting confidence level obtained with the CLs technique for each
        signal region of the 8 TeV CMS mono-jet analysis and chooses the most
        constraining one.
        """
        self._initialize()
        return low_level.monojet()

    def mass_left_selectron(self) -> float:
        """
        Compute the left-handed Selectron mass.
        """
        return low_level.find_val("MSeL")

    def mass_right_selectron(self) -> float:
        """
        Compute the right-handed Selectron mass.
        """
        return low_level.find_val("MSeR")

    def mass_electron_sneutrino(self) -> float:
        """
        Compute the electron Snuetrino mass.
        """
        return low_level.find_val("MSne")

    def mass_left_smuon(self) -> float:
        """
        Compute the left-handed smuon mass.
        """
        return low_level.find_val("MSmL")

    def mass_right_smuon(self) -> float:
        """
        Compute the right-handed smuon mass.
        """
        return low_level.find_val("MSmR")

    def mass_muon_sneutrino(self) -> float:
        """
        Compute the muon snuetrino mass.
        """
        return low_level.find_val("MSnm")

    def mass_stau1(self) -> float:
        """
        Compute the light stau mass.
        """
        return low_level.find_val("MSl1")

    def mass_stau2(self) -> float:
        """
        Compute the heavy stau mass.
        """
        return low_level.find_val("MSl2")

    def mass_stau_neutrino(self) -> float:
        """
        Compute the tau snuetrino mass.
        """
        return low_level.find_val("MSnl")

    def mass_left_up_squark(self) -> float:
        """
        Compute the left-handed up-squark mass.
        """
        return low_level.find_val("MSuL")

    def mass_right_up_squark(self) -> float:
        """
        Compute the right-handed up-squark mass.
        """
        return low_level.find_val("MSuR")

    def mass_left_down_squark(self) -> float:
        """
        Compute the left-handed down-squark mass.
        """
        return low_level.find_val("MSdL")

    def mass_right_down_squark(self) -> float:
        """
        Compute the right-handed down-squark mass.
        """
        return low_level.find_val("MSdR")

    def mass_left_charm_squark(self) -> float:
        """
        Compute the left-handed charm-squark mass.
        """
        return low_level.find_val("MScL")

    def mass_right_charm_squark(self) -> float:
        """
        Compute the right-handed charm-squark mass.
        """
        return low_level.find_val("MScR")

    def mass_left_strange_squark(self) -> float:
        """
        Compute the left-handed strange-squark mass.
        """
        return low_level.find_val("MSsL")

    def mass_right_strange_squark(self) -> float:
        """
        Compute the right-handed strange-squark mass.
        """
        return low_level.find_val("MSsR")

    def mass_top_squark1(self) -> float:
        """
        Compute the light top-squark mass.
        """
        return low_level.find_val("MSt1")

    def mass_top_squark2(self) -> float:
        """
        Compute the heavy top-squark mass.
        """
        return low_level.find_val("MSt2")

    def mass_bottom_squark1(self) -> float:
        """
        Compute the light bottom-squark mass.
        """
        return low_level.find_val("MSb1")

    def mass_bottom_squark2(self) -> float:
        """
        Compute the heavy bottom-squark mass.
        """
        return low_level.find_val("MSb2")

    def mass_neutralino1(self) -> float:
        """
        Compute the 1st neutralino mass.
        """
        return low_level.find_val("MNE1")

    def mass_neutralino2(self) -> float:
        """
        Compute the 2nd neutralino mass.
        """
        return low_level.find_val("MNE2")

    def mass_neutralino3(self) -> float:
        """
        Compute the 3rd neutralino mass.
        """
        return low_level.find_val("MNE3")

    def mass_neutralino4(self) -> float:
        """
        Compute the 4th neutralino mass.
        """
        return low_level.find_val("MNE4")

    def mass_chargino1(self) -> float:
        """
        Compute the light chargino mass.
        """
        return low_level.find_val("MC1")

    def mass_chargino2(self) -> float:
        """
        Compute the heavy chargino mass.
        """
        return low_level.find_val("MC2")

    def mass_higgs(self) -> float:
        """
        Compute the SM Higgs boson mass.
        """
        return low_level.find_val("Mh")

    def mass_higgs2(self) -> float:
        """
        Compute the 2nd Higgs boson mass.
        """
        return low_level.find_val("MH")

    def mass_pseudo_scal_higgs(self) -> float:
        """
        Compute the pseudo-scalar Higgs boson mass.
        """
        return low_level.find_val("MH3")

    def mass_charged_higgs(self) -> float:
        """
        Compute the charged Higgs boson mass.
        """
        return low_level.find_val("MHc")

    def mass_glunio(self) -> float:
        """
        Compute the gluino mass.
        """
        return low_level.find_val("MSG")
