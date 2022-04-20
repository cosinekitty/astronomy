//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[nextLunarApsis](next-lunar-apsis.md)

# nextLunarApsis

[jvm]\
fun [nextLunarApsis](next-lunar-apsis.md)(apsis: [ApsisInfo](-apsis-info/index.md)): [ApsisInfo](-apsis-info/index.md)

Finds the next lunar perigee or apogee event in a series.

Finds the next consecutive time the Moon is closest or farthest from the Earth in its orbit. See [searchLunarApsis](search-lunar-apsis.md) for more details.

## Parameters

jvm

| | |
|---|---|
| apsis | An [ApsisInfo] value obtained from a call     to [searchLunarApsis] or `nextLunarApsis`. |