//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[SeasonsInfo](index.md)

# SeasonsInfo

[jvm]\
class [SeasonsInfo](index.md)(mar_equinox: [AstroTime](../-astro-time/index.md), jun_solstice: [AstroTime](../-astro-time/index.md), sep_equinox: [AstroTime](../-astro-time/index.md), dec_solstice: [AstroTime](../-astro-time/index.md))

The dates and times of changes of season for a given calendar year.

Call #Astronomy.seasons to calculate this data structure for a given year.

## Constructors

| | |
|---|---|
| [SeasonsInfo](-seasons-info.md) | [jvm]<br>fun [SeasonsInfo](-seasons-info.md)(mar_equinox: [AstroTime](../-astro-time/index.md), jun_solstice: [AstroTime](../-astro-time/index.md), sep_equinox: [AstroTime](../-astro-time/index.md), dec_solstice: [AstroTime](../-astro-time/index.md)) |

## Properties

| Name | Summary |
|---|---|
| [dec_solstice](dec_solstice.md) | [jvm]<br>val [dec_solstice](dec_solstice.md): [AstroTime](../-astro-time/index.md)<br>The date and time of the December solstice for the specified year. |
| [jun_solstice](jun_solstice.md) | [jvm]<br>val [jun_solstice](jun_solstice.md): [AstroTime](../-astro-time/index.md)<br>The date and time of the June soltice for the specified year. |
| [mar_equinox](mar_equinox.md) | [jvm]<br>val [mar_equinox](mar_equinox.md): [AstroTime](../-astro-time/index.md)<br>The date and time of the March equinox for the specified year. |
| [sep_equinox](sep_equinox.md) | [jvm]<br>val [sep_equinox](sep_equinox.md): [AstroTime](../-astro-time/index.md)<br>The date and time of the September equinox for the specified year. |
