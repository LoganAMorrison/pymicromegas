from typing import List, overload

class SugraParameters:
    def __repr__(self) -> str: ...
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(
        self, m0: float, mhf: float, a0: float, tb: float, sgn: float
    ) -> None: ...
    @property
    def m0(self) -> float: ...
    @m0.setter
    def m0(self, val) -> None: ...
    @property
    def mhf(self) -> float: ...
    @mhf.setter
    def mhf(self, val) -> None: ...
    @property
    def a0(self) -> float: ...
    @a0.setter
    def a0(self, val) -> None: ...
    @property
    def tb(self) -> float: ...
    @tb.setter
    def tb(self, val) -> None: ...
    @property
    def sgn(self) -> float: ...
    @sgn.setter
    def sgn(self, val) -> None: ...

class EwsbParameters:
    def __repr__(self) -> str: ...
    @overload
    def __init__(self) -> None: ...
    @overload
    def __init__(
        self,
        mu: float,
        mg1: float,
        mg2: float,
        mg3: float,
        ml1: float,
        ml2: float,
        ml3: float,
        mr1: float,
        mr2: float,
        mr3: float,
        mq1: float,
        mq2: float,
        mq3: float,
        mu1: float,
        mu2: float,
        mu3: float,
        md1: float,
        md2: float,
        md3: float,
        mh3: float,
        tb: float,
        at: float,
        ab: float,
        al: float,
    ) -> None: ...
    @property
    def mg1(self) -> float: ...
    @mg1.setter
    def mg1(self, val: float) -> None: ...
    @property
    def mg2(self) -> float: ...
    @mg2.setter
    def mg2(self, val: float) -> None: ...
    @property
    def mg3(self) -> float: ...
    @mg3.setter
    def mg3(self, val: float) -> None: ...
    @property
    def ml1(self) -> float: ...
    @ml1.setter
    def ml1(self, val: float) -> None: ...
    @property
    def ml2(self) -> float: ...
    @ml2.setter
    def ml2(self, val: float) -> None: ...
    @property
    def ml3(self) -> float: ...
    @ml3.setter
    def ml3(self, val: float) -> None: ...
    @property
    def mr1(self) -> float: ...
    @mr1.setter
    def mr1(self, val: float) -> None: ...
    @property
    def mr2(self) -> float: ...
    @mr2.setter
    def mr2(self, val: float) -> None: ...
    @property
    def mr3(self) -> float: ...
    @mr3.setter
    def mr3(self, val: float) -> None: ...
    @property
    def mq1(self) -> float: ...
    @mq1.setter
    def mq1(self, val: float) -> None: ...
    @property
    def mq2(self) -> float: ...
    @mq2.setter
    def mq2(self, val: float) -> None: ...
    @property
    def mq3(self) -> float: ...
    @mq3.setter
    def mq3(self, val: float) -> None: ...
    @property
    def mu1(self) -> float: ...
    @mu1.setter
    def mu1(self, val: float) -> None: ...
    @property
    def mu2(self) -> float: ...
    @mu2.setter
    def mu2(self, val: float) -> None: ...
    @property
    def mu3(self) -> float: ...
    @mu3.setter
    def mu3(self, val: float) -> None: ...
    @property
    def md1(self) -> float: ...
    @md1.setter
    def md1(self, val: float) -> None: ...
    @property
    def md2(self) -> float: ...
    @md2.setter
    def md2(self, val: float) -> None: ...
    @property
    def md3(self) -> float: ...
    @md3.setter
    def md3(self, val: float) -> None: ...
    @property
    def mu(self) -> float: ...
    @mu.setter
    def mu(self, val: float) -> None: ...
    @property
    def mh3(self) -> float: ...
    @mh3.setter
    def mh3(self, val: float) -> None: ...
    @property
    def tb(self) -> float: ...
    @tb.setter
    def tb(self, val: float) -> None: ...
    @property
    def at(self) -> float: ...
    @at.setter
    def at(self, val: float) -> None: ...
    @property
    def ab(self) -> float: ...
    @ab.setter
    def ab(self, val: float) -> None: ...
    @property
    def al(self) -> float: ...
    @al.setter
    def al(self, val: float) -> None: ...
    @property
    def au(self) -> float: ...
    @au.setter
    def au(self, val: float) -> None: ...
    @property
    def ad(self) -> float: ...
    @ad.setter
    def ad(self, val: float) -> None: ...
    @property
    def am(self) -> float: ...
    @am.setter
    def am(self, val: float) -> None: ...

class MicromegasResults:
    def __repr__(self) -> str: ...
    def omega(self) -> List[float]: ...
    def xf(self) -> List[float]: ...
    def bsgsm(self) -> List[float]: ...
    def bsgnlo(self) -> List[float]: ...
    def deltarho(self) -> List[float]: ...
    def bsmumu(self) -> List[float]: ...
    def bstaunu(self) -> List[float]: ...
    def gmuon(self) -> List[float]: ...
    def rl23(self) -> List[float]: ...
    def dtaunu(self) -> List[float]: ...
    def dmunu(self) -> List[float]: ...
    def masslimits(self) -> List[float]: ...
    def proton_si_amp(self) -> List[float]: ...
    def proton_sd_amp(self) -> List[float]: ...
    def neutron_si_amp(self) -> List[float]: ...
    def neutron_sd_amp(self) -> List[float]: ...
    def pval_xenon1T(self) -> List[float]: ...
    def pval_cresst(self) -> List[float]: ...
    def pval_darkside(self) -> List[float]: ...
    def pval_pico(self) -> List[float]: ...
    def msel(self) -> List[float]: ...
    def mser(self) -> List[float]: ...
    def msne(self) -> List[float]: ...
    def msml(self) -> List[float]: ...
    def msmr(self) -> List[float]: ...
    def msnm(self) -> List[float]: ...
    def msl1(self) -> List[float]: ...
    def msl2(self) -> List[float]: ...
    def msnl(self) -> List[float]: ...
    def msul(self) -> List[float]: ...
    def msur(self) -> List[float]: ...
    def msdl(self) -> List[float]: ...
    def msdr(self) -> List[float]: ...
    def mscl(self) -> List[float]: ...
    def mscr(self) -> List[float]: ...
    def mssl(self) -> List[float]: ...
    def mssr(self) -> List[float]: ...
    def mst1(self) -> List[float]: ...
    def mst2(self) -> List[float]: ...
    def msb1(self) -> List[float]: ...
    def msb2(self) -> List[float]: ...
    def mneut1(self) -> List[float]: ...
    def mneut2(self) -> List[float]: ...
    def mneut3(self) -> List[float]: ...
    def mneut4(self) -> List[float]: ...
    def mchg1(self) -> List[float]: ...
    def mchg2(self) -> List[float]: ...
    def mhsm(self) -> List[float]: ...
    def mh(self) -> List[float]: ...
    def mhc(self) -> List[float]: ...
    def mg(self) -> List[float]: ...

class MicromegasSettings:
    def __init__(
        self: MicromegasSettings,
        relic_density: bool = True,
        masses: bool = True,
        gmuon: bool = True,
        bsg: bool = False,
        bsmumu: bool = False,
        btaunu: bool = False,
        deltarho: bool = False,
        rl23: bool = False,
        d_taunu_and_munu: bool = False,
        masslimits: bool = False,
        nucleon_amplitudes: bool = False,
        direct_detection_pvalues: bool = False,
        z_invisible: bool = False,
        lsp_nlsp_lep: bool = False,
        z_prime_limits: bool = False,
        monojet: bool = False,
        fast: bool = True,
        beps: float = 0.0001,
    ) -> None: ...
    def __repr__(self) -> str: ...
    @property
    def relic_density(self) -> bool: ...
    @relic_density.setter
    def relic_density(self, val: bool) -> None: ...
    @property
    def masses(self) -> bool: ...
    @masses.setter
    def masses(self, val: bool) -> None: ...
    @property
    def gmuon(self) -> bool: ...
    @gmuon.setter
    def gmuon(self, val: bool) -> None: ...
    @property
    def bsg(self) -> bool: ...
    @bsg.setter
    def bsg(self, val: bool) -> None: ...
    @property
    def bsmumu(self) -> bool: ...
    @bsmumu.setter
    def bsmumu(self, val: bool) -> None: ...
    @property
    def btaunu(self) -> bool: ...
    @btaunu.setter
    def btaunu(self, val: bool) -> None: ...
    @property
    def deltarho(self) -> bool: ...
    @deltarho.setter
    def deltarho(self, val: bool) -> None: ...
    @property
    def rl23(self) -> bool: ...
    @rl23.setter
    def rl23(self, val: bool) -> None: ...
    @property
    def d_taunu_and_munu(self) -> bool: ...
    @d_taunu_and_munu.setter
    def d_taunu_and_munu(self, val: bool) -> None: ...
    @property
    def masslimits(self) -> bool: ...
    @masslimits.setter
    def masslimits(self, val: bool) -> None: ...
    @property
    def nucleon_amplitudes(self) -> bool: ...
    @nucleon_amplitudes.setter
    def nucleon_amplitudes(self, val: bool) -> None: ...
    @property
    def direct_detection_pvalues(self) -> bool: ...
    @direct_detection_pvalues.setter
    def direct_detection_pvalues(self, val: bool) -> None: ...
    @property
    def z_invisible(self) -> bool: ...
    @z_invisible.setter
    def z_invisible(self, val: bool) -> None: ...
    @property
    def lsp_nlsp_lep(self) -> bool: ...
    @lsp_nlsp_lep.setter
    def lsp_nlsp_lep(self, val: bool) -> None: ...
    @property
    def z_prime_limits(self) -> bool: ...
    @z_prime_limits.setter
    def z_prime_limits(self, val: bool) -> None: ...
    @property
    def monojet(self) -> bool: ...
    @monojet.setter
    def monojet(self, val: bool) -> None: ...
    @property
    def fast(self) -> bool: ...
    @fast.setter
    def fast(self, val: bool) -> None: ...
    @property
    def beps(self) -> float: ...
    @beps.setter
    def beps(self, val: float) -> None: ...
