from typing import Union, Optional

from .pymicromegas import (
    EwsbParameters,
    SugraParameters,
    MicromegasResults,
    MicromegasSettings,
)
from .spheno import spheno
from .softsusy import softsusy
from .suspect import suspect
from . import micromegas as low_level

EwsbOrSugra = Union[EwsbParameters, SugraParameters]


class LowLevelMicromegasEwsb:
    def __init__(self, params: EwsbParameters) -> None:
        self._params: EwsbParameters = params

    @property
    def mg1(self) -> float:
        return low_level.find_val("MG1")

    @mg1.setter
    def mg1(self, val: float) -> None:
        low_level.assign_val("MG1", val)

    @property
    def mg2(self) -> float:
        ...

    @mg2.setter
    def mg2(self, val: float) -> None:
        ...

    @property
    def mg3(self) -> float:
        ...

    @mg3.setter
    def mg3(self, val: float) -> None:
        ...

    @property
    def ml1(self) -> float:
        ...

    @ml1.setter
    def ml1(self, val: float) -> None:
        ...

    @property
    def ml2(self) -> float:
        ...

    @ml2.setter
    def ml2(self, val: float) -> None:
        ...

    @property
    def ml3(self) -> float:
        ...

    @ml3.setter
    def ml3(self, val: float) -> None:
        ...

    @property
    def mr1(self) -> float:
        ...

    @mr1.setter
    def mr1(self, val: float) -> None:
        ...

    @property
    def mr2(self) -> float:
        ...

    @mr2.setter
    def mr2(self, val: float) -> None:
        ...

    @property
    def mr3(self) -> float:
        ...

    @mr3.setter
    def mr3(self, val: float) -> None:
        ...

    @property
    def mq1(self) -> float:
        ...

    @mq1.setter
    def mq1(self, val: float) -> None:
        ...

    @property
    def mq2(self) -> float:
        ...

    @mq2.setter
    def mq2(self, val: float) -> None:
        ...

    @property
    def mq3(self) -> float:
        ...

    @mq3.setter
    def mq3(self, val: float) -> None:
        ...

    @property
    def mu1(self) -> float:
        ...

    @mu1.setter
    def mu1(self, val: float) -> None:
        ...

    @property
    def mu2(self) -> float:
        ...

    @mu2.setter
    def mu2(self, val: float) -> None:
        ...

    @property
    def mu3(self) -> float:
        ...

    @mu3.setter
    def mu3(self, val: float) -> None:
        ...

    @property
    def md1(self) -> float:
        ...

    @md1.setter
    def md1(self, val: float) -> None:
        ...

    @property
    def md2(self) -> float:
        ...

    @md2.setter
    def md2(self, val: float) -> None:
        ...

    @property
    def md3(self) -> float:
        ...

    @md3.setter
    def md3(self, val: float) -> None:
        ...

    @property
    def mu(self) -> float:
        ...

    @mu.setter
    def mu(self, val: float) -> None:
        ...

    @property
    def mh3(self) -> float:
        ...

    @mh3.setter
    def mh3(self, val: float) -> None:
        ...

    @property
    def tb(self) -> float:
        ...

    @tb.setter
    def tb(self, val: float) -> None:
        ...

    @property
    def at(self) -> float:
        ...

    @at.setter
    def at(self, val: float) -> None:
        ...

    @property
    def ab(self) -> float:
        ...

    @ab.setter
    def ab(self, val: float) -> None:
        ...

    @property
    def al(self) -> float:
        ...

    @al.setter
    def al(self, val: float) -> None:
        ...

    @property
    def au(self) -> float:
        ...

    @au.setter
    def au(self, val: float) -> None:
        ...

    @property
    def ad(self) -> float:
        ...

    @ad.setter
    def ad(self, val: float) -> None:
        ...

    @property
    def am(self) -> float:
        ...

    @am.setter
    def am(self, val: float) -> None:
        ...
