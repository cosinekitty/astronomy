//[astronomy](../../../index.md)/[io.github.cosinekitty.astronomy](../index.md)/[SearchContext](index.md)

# SearchContext

[jvm]\
interface [SearchContext](index.md)

Represents a function whose ascending root is to be found.

This interface must be implemented for callers of [Astronomy.search](../-astronomy/search.md) in order to find the ascending root of a smooth function. A class that implements SearchContext can hold state information needed to evaluate the scalar function eval.

## Functions

| Name | Summary |
|---|---|
| [eval](eval.md) | [jvm]<br>abstract fun [eval](eval.md)(time: [AstroTime](../-astro-time/index.md)): [Double](https://kotlinlang.org/api/latest/jvm/stdlib/kotlin/-double/index.html)<br>Evaluates a scalar function at a given time. |
