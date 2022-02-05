from typing import Union

from .pymicromegas import (
    EwsbParameters,
    SugraParameters,
    MicromegasResults,
    MicromegasSettings,
)
from .spheno import spheno
from .suspect import suspect
from . import micromegas as _low_level
from .__low_level_base import LowLevelMicromegasBase as _LowLevelMicromegasBase

try:
    from .softsusy import softsusy
except ImportError:
    pass


EwsbOrSugra = Union[EwsbParameters, SugraParameters]


class LowLevelMicromegasEwsb(_LowLevelMicromegasBase):
    def __init__(self, params: EwsbParameters) -> None:
        self._params: EwsbParameters = params
        self._initialized = False
        _low_level.assign_val("MG1", self._params.mg1)
        _low_level.assign_val("MG2", self._params.mg2)
        _low_level.assign_val("MG3", self._params.mg3)
        _low_level.assign_val("Ml1", self._params.ml1)
        _low_level.assign_val("Ml2", self._params.ml2)
        _low_level.assign_val("Ml3", self._params.ml3)
        _low_level.assign_val("Mr1", self._params.mr1)
        _low_level.assign_val("Mr2", self._params.mr2)
        _low_level.assign_val("Mr3", self._params.mr3)
        _low_level.assign_val("Mq1", self._params.mq1)
        _low_level.assign_val("Mq2", self._params.mq2)
        _low_level.assign_val("Mq3", self._params.mq3)
        _low_level.assign_val("Mu1", self._params.mu1)
        _low_level.assign_val("Mu2", self._params.mu2)
        _low_level.assign_val("Mu3", self._params.mu3)
        _low_level.assign_val("Md1", self._params.md1)
        _low_level.assign_val("Md2", self._params.md2)
        _low_level.assign_val("Md3", self._params.md3)
        _low_level.assign_val("mu", self._params.mu)
        _low_level.assign_val("MH3", self._params.mh3)
        _low_level.assign_val("tb", self._params.tb)
        _low_level.assign_val("At", self._params.at)
        _low_level.assign_val("Ab", self._params.ab)
        _low_level.assign_val("Al", self._params.al)
        super().__init__()

    def _initialize(self):
        if not self._initialized:
            _low_level.mssm_ewsb()
            _low_level.sort_odd_particles()
            self._initialized = True

    @property
    def mg1(self) -> float:
        return self._params.mg1

    @mg1.setter
    def mg1(self, val: float) -> None:
        self._params.mg1 = val
        _low_level.assign_val("MG1", self._params.mg1)
        self._initialized = False

    @property
    def mg2(self) -> float:
        return self._params.mg2

    @mg2.setter
    def mg2(self, val: float) -> None:
        self._params.mg2 = val
        _low_level.assign_val("MG2", self._params.mg2)
        self._initialized = False

    @property
    def mg3(self) -> float:
        return self._params.mg3

    @mg3.setter
    def mg3(self, val: float) -> None:
        self._params.mg3 = val
        _low_level.assign_val("MG3", self._params.mg3)
        self._initialized = False

    @property
    def ml1(self) -> float:
        return self._params.ml1

    @ml1.setter
    def ml1(self, val: float) -> None:
        self._params.ml1 = val
        _low_level.assign_val("Ml1", self._params.ml1)
        self._initialized = False

    @property
    def ml2(self) -> float:
        return self._params.ml2

    @ml2.setter
    def ml2(self, val: float) -> None:
        self._params.ml2 = val
        _low_level.assign_val("Ml2", self._params.ml2)
        self._initialized = False

    @property
    def ml3(self) -> float:
        return self._params.ml3

    @ml3.setter
    def ml3(self, val: float) -> None:
        self._params.ml3 = val
        _low_level.assign_val("Ml3", self._params.ml3)
        self._initialized = False

    @property
    def mr1(self) -> float:
        return self._params.mr1

    @mr1.setter
    def mr1(self, val: float) -> None:
        self._params.mr1 = val
        _low_level.assign_val("Mr1", self._params.mr1)
        self._initialized = False

    @property
    def mr2(self) -> float:
        return self._params.mr2

    @mr2.setter
    def mr2(self, val: float) -> None:
        self._params.mr2 = val
        _low_level.assign_val("Mr2", self._params.mr2)
        self._initialized = False

    @property
    def mr3(self) -> float:
        return self._params.mr3

    @mr3.setter
    def mr3(self, val: float) -> None:
        self._params.mr3 = val
        _low_level.assign_val("Mr3", self._params.mr3)
        self._initialized = False

    @property
    def mq1(self) -> float:
        return self._params.mq1

    @mq1.setter
    def mq1(self, val: float) -> None:
        self._params.mq1 = val
        _low_level.assign_val("Mq1", self._params.mq1)
        self._initialized = False

    @property
    def mq2(self) -> float:
        return self._params.mq2

    @mq2.setter
    def mq2(self, val: float) -> None:
        self._params.mq2 = val
        _low_level.assign_val("Mq2", self._params.mq2)
        self._initialized = False

    @property
    def mq3(self) -> float:
        return self._params.mq3

    @mq3.setter
    def mq3(self, val: float) -> None:
        self._params.mq3 = val
        _low_level.assign_val("Mq3", self._params.mq3)
        self._initialized = False

    @property
    def mu1(self) -> float:
        return self._params.mu1

    @mu1.setter
    def mu1(self, val: float) -> None:
        self._params.mu1 = val
        _low_level.assign_val("Mu1", self._params.mu1)
        self._initialized = False

    @property
    def mu2(self) -> float:
        return self._params.mu2

    @mu2.setter
    def mu2(self, val: float) -> None:
        self._params.mu2 = val
        _low_level.assign_val("Mu2", self._params.mu2)
        self._initialized = False

    @property
    def mu3(self) -> float:
        return self._params.mu3

    @mu3.setter
    def mu3(self, val: float) -> None:
        self._params.mu3 = val
        _low_level.assign_val("Mu3", self._params.mu3)
        self._initialized = False

    @property
    def md1(self) -> float:
        return self._params.md1

    @md1.setter
    def md1(self, val: float) -> None:
        self._params.md1 = val
        _low_level.assign_val("Md1", self._params.md1)
        self._initialized = False

    @property
    def md2(self) -> float:
        return self._params.md2

    @md2.setter
    def md2(self, val: float) -> None:
        self._params.md2 = val
        _low_level.assign_val("Md2", self._params.md2)
        self._initialized = False

    @property
    def md3(self) -> float:
        return self._params.md3

    @md3.setter
    def md3(self, val: float) -> None:
        self._params.md3 = val
        _low_level.assign_val("Md3", self._params.md3)
        self._initialized = False

    @property
    def mu(self) -> float:
        return self._params.mu

    @mu.setter
    def mu(self, val: float) -> None:
        self._params.mu = val
        _low_level.assign_val("mu", self._params.mu)
        self._initialized = False

    @property
    def mh3(self) -> float:
        return self._params.mh3

    @mh3.setter
    def mh3(self, val: float) -> None:
        self._params.mh3 = val
        _low_level.assign_val("MH3", self._params.mh3)
        self._initialized = False

    @property
    def tb(self) -> float:
        return self._params.tb

    @tb.setter
    def tb(self, val: float) -> None:
        self._params.tb = val
        _low_level.assign_val("tb", self._params.tb)
        self._initialized = False

    @property
    def at(self) -> float:
        return self._params.at

    @at.setter
    def at(self, val: float) -> None:
        self._params.at = val
        _low_level.assign_val("At", self._params.at)
        self._initialized = False

    @property
    def ab(self) -> float:
        return self._params.ab

    @ab.setter
    def ab(self, val: float) -> None:
        self._params.ab = val
        _low_level.assign_val("Ab", self._params.ab)
        self._initialized = False

    @property
    def al(self) -> float:
        return self._params.al

    @al.setter
    def al(self, val: float) -> None:
        self._params.al = val
        _low_level.assign_val("Al", self._params.al)
        self._initialized = False

    @property
    def au(self) -> float:
        return self._params.au

    @au.setter
    def au(self, val: float) -> None:
        self._params.au = val
        _low_level.assign_val("Au", self._params.au)
        self._initialized = False

    @property
    def ad(self) -> float:
        return self._params.ad

    @ad.setter
    def ad(self, val: float) -> None:
        self._params.ad = val
        _low_level.assign_val("Ad", self._params.ad)
        self._initialized = False

    @property
    def am(self) -> float:
        return self._params.am

    @am.setter
    def am(self, val: float) -> None:
        self._params.am = val
        _low_level.assign_val("Am", self._params.am)
        self._initialized = False

    def verify_parameter_state(self, rtol=1e-5):
        """
        Verify that the state of micrOMEGAs is the same as the class state.
        """

        def frac_diff(val1: float, val2: float) -> float:
            return (val1 - val2) / val1

        def verify(val1: float, val2: float) -> None:
            fdiff = frac_diff(val1, val2)
            assert fdiff <= rtol

        assert verify(self._params.mg1, _low_level.find_val("MG1"))
        assert verify(self._params.mg2, _low_level.find_val("MG2"))
        assert verify(self._params.mg3, _low_level.find_val("MG3"))
        assert verify(self._params.ml1, _low_level.find_val("Ml1"))
        assert verify(self._params.ml2, _low_level.find_val("Ml2"))
        assert verify(self._params.ml3, _low_level.find_val("Ml3"))
        assert verify(self._params.mr1, _low_level.find_val("Mr1"))
        assert verify(self._params.mr2, _low_level.find_val("Mr2"))
        assert verify(self._params.mr3, _low_level.find_val("Mr3"))
        assert verify(self._params.mq1, _low_level.find_val("Mq1"))
        assert verify(self._params.mq2, _low_level.find_val("Mq2"))
        assert verify(self._params.mq3, _low_level.find_val("Mq3"))
        assert verify(self._params.mu1, _low_level.find_val("Mu1"))
        assert verify(self._params.mu2, _low_level.find_val("Mu2"))
        assert verify(self._params.mu3, _low_level.find_val("Mu3"))
        assert verify(self._params.md1, _low_level.find_val("Md1"))
        assert verify(self._params.md2, _low_level.find_val("Md2"))
        assert verify(self._params.md3, _low_level.find_val("Md3"))
        assert verify(self._params.mu, _low_level.find_val("mu"))
        assert verify(self._params.mh3, _low_level.find_val("MH3"))
        assert verify(self._params.tb, _low_level.find_val("tb"))
        assert verify(self._params.at, _low_level.find_val("At"))
        assert verify(self._params.ab, _low_level.find_val("Ab"))
        assert verify(self._params.al, _low_level.find_val("Al"))


class LowLevelMicromegasSugra(_LowLevelMicromegasBase):
    def __init__(self, params: SugraParameters) -> None:
        self._params: SugraParameters = params
        self._initialized = False
        super().__init__()

    def _initialize(self):
        if not self._initialized:
            mhf = self._params.mhf
            a0 = self._params.a0
            m0 = self._params.m0
            tb = self._params.tb
            sgn = self._params.sgn
            _low_level.mssm_sugra(
                tb,
                mhf,
                mhf,
                mhf,
                a0,
                a0,
                a0,
                sgn,
                m0,
                m0,
                m0,
                m0,
                m0,
                m0,
                m0,
                m0,
                m0,
                m0,
                m0,
                m0,
            )
            _low_level.sort_odd_particles()
            self._initialized = True

    @property
    def m0(self) -> float:
        return self._params.m0

    @m0.setter
    def m0(self, val: float) -> None:
        self._params.m0 = val
        self._initialized = False

    @property
    def mhf(self) -> float:
        return self._params.mhf

    @mhf.setter
    def mhf(self, val: float) -> None:
        self._params.mhf = val
        self._initialized = False

    @property
    def a0(self) -> float:
        return self._params.a0

    @a0.setter
    def a0(self, val: float) -> None:
        self._params.a0 = val
        self._initialized = False

    @property
    def tb(self) -> float:
        return self._params.tb

    @tb.setter
    def tb(self, val: float) -> None:
        self._params.tb = val
        self._initialized = False

    @property
    def sgn(self) -> float:
        return self._params.sgn

    @sgn.setter
    def sgn(self, val: float) -> None:
        self._params.sgn = val
        self._initialized = False
