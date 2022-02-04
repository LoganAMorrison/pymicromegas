from typing import List, overload
from pymicromegas import (
    SugraParameters,
    MicromegasResults,
    MicromegasSettings,
    EwsbParameters,
)

@overload
def softsusy(
    params: SugraParameters, settings: MicromegasSettings
) -> MicromegasResults: ...
@overload
def softsusy(params: SugraParameters) -> MicromegasResults: ...
@overload
def softsusy(
    params: List[SugraParameters], settings: MicromegasSettings
) -> MicromegasResults: ...
@overload
def softsusy(params: List[SugraParameters]) -> MicromegasResults: ...
@overload
def softsusy(
    params: EwsbParameters, settings: MicromegasSettings
) -> MicromegasResults: ...
@overload
def softsusy(params: EwsbParameters) -> MicromegasResults: ...
@overload
def softsusy(
    params: List[EwsbParameters], settings: MicromegasSettings
) -> MicromegasResults: ...
@overload
def softsusy(params: List[EwsbParameters]) -> MicromegasResults: ...
