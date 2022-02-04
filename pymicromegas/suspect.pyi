from typing import List, overload
from pymicromegas import (
    SugraParameters,
    MicromegasResults,
    MicromegasSettings,
    EwsbParameters,
)

@overload
def suspect(
    params: SugraParameters, settings: MicromegasSettings
) -> MicromegasResults: ...
@overload
def suspect(params: SugraParameters) -> MicromegasResults: ...
@overload
def suspect(
    params: List[SugraParameters], settings: MicromegasSettings
) -> MicromegasResults: ...
@overload
def suspect(params: List[SugraParameters]) -> MicromegasResults: ...
@overload
def suspect(
    params: EwsbParameters, settings: MicromegasSettings
) -> MicromegasResults: ...
@overload
def suspect(params: EwsbParameters) -> MicromegasResults: ...
@overload
def suspect(
    params: List[EwsbParameters], settings: MicromegasSettings
) -> MicromegasResults: ...
@overload
def suspect(params: List[EwsbParameters]) -> MicromegasResults: ...
