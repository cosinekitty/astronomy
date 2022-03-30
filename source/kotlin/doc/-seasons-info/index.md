//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[SeasonsInfo](index.md)

# SeasonsInfo

[jvm]\
class [SeasonsInfo](index.md)(marEquinox: [AstroTime](../-astro-time/index.md), junSolstice: [AstroTime](../-astro-time/index.md), sepEquinox: [AstroTime](../-astro-time/index.md), decSolstice: [AstroTime](../-astro-time/index.md))

The dates and times of changes of season for a given calendar year.

Call Astronomy.seasons to calculate this data structure for a given year.

## Constructors

| | |
|---|---|
| [SeasonsInfo](-seasons-info.md) | [jvm]<br>fun [SeasonsInfo](-seasons-info.md)(marEquinox: [AstroTime](../-astro-time/index.md), junSolstice: [AstroTime](../-astro-time/index.md), sepEquinox: [AstroTime](../-astro-time/index.md), decSolstice: [AstroTime](../-astro-time/index.md)) |

## Properties

| Name | Summary |
|---|---|
| [decSolstice](dec-solstice.md) | [jvm]<br>val [decSolstice](dec-solstice.md): [AstroTime](../-astro-time/index.md)<br>The date and time of the December solstice for the specified year. |
| [junSolstice](jun-solstice.md) | [jvm]<br>val [junSolstice](jun-solstice.md): [AstroTime](../-astro-time/index.md)<br>The date and time of the June soltice for the specified year. |
| [marEquinox](mar-equinox.md) | [jvm]<br>val [marEquinox](mar-equinox.md): [AstroTime](../-astro-time/index.md)<br>The date and time of the March equinox for the specified year. |
| [sepEquinox](sep-equinox.md) | [jvm]<br>val [sepEquinox](sep-equinox.md): [AstroTime](../-astro-time/index.md)<br>The date and time of the September equinox for the specified year. |
