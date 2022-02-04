from typing import List, overload
from pymicromegas import (
    SugraParameters,
    MicromegasResults,
    MicromegasSettings,
    EwsbParameters,
)

@overload
def spheno(
    params: SugraParameters, settings: MicromegasSettings
) -> MicromegasResults: ...
@overload
def spheno(params: SugraParameters) -> MicromegasResults: ...
@overload
def spheno(
    params: List[SugraParameters], settings: MicromegasSettings
) -> MicromegasResults: ...
@overload
def spheno(params: List[SugraParameters]) -> MicromegasResults: ...
@overload
def spheno(
    params: EwsbParameters, settings: MicromegasSettings
) -> MicromegasResults: ...
@overload
def spheno(params: EwsbParameters) -> MicromegasResults: ...
@overload
def spheno(
    params: List[EwsbParameters], settings: MicromegasSettings
) -> MicromegasResults: ...
@overload
def spheno(params: List[EwsbParameters]) -> MicromegasResults: ...
