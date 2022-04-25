//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[searchMoonQuarter](search-moon-quarter.md)

# searchMoonQuarter

fun [searchMoonQuarter](search-moon-quarter.md)(startTime: [Time](-time/index.md)): [MoonQuarterInfo](-moon-quarter-info/index.md)

Finds the first lunar quarter after the specified date and time. A lunar quarter is one of the following four lunar phase events: new moon, first quarter, full moon, third quarter. This function finds the lunar quarter that happens soonest after the specified date and time.

To continue iterating through consecutive lunar quarters, call this function once, followed by calls to #NextMoonQuarter as many times as desired.

#### Return

A [MoonQuarterInfo](-moon-quarter-info/index.md) object reporting the next quarter phase and the time it will occur.

## Parameters

| | |
|---|---|
| startTime | The date and time at which to start the search. |
