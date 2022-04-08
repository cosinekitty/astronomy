//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[nextMoonQuarter](next-moon-quarter.md)

# nextMoonQuarter

[jvm]\
fun [nextMoonQuarter](next-moon-quarter.md)(mq: [MoonQuarterInfo](-moon-quarter-info/index.md)): [MoonQuarterInfo](-moon-quarter-info/index.md)

Continues searching for lunar quarters from a previous search.

After calling [searchMoonQuarter](search-moon-quarter.md), this function can be called one or more times to continue finding consecutive lunar quarters. This function finds the next consecutive moon quarter event after the one passed in as the parameter mq.

#### Return

The moon quarter that occurs next in time after the one passed in mq.

## Parameters

jvm

| | |
|---|---|
| The | previous moon quarter found by a call to [searchMoonQuarter](search-moon-quarter.md) or nextMoonQuarter. |
