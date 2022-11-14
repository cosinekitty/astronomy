//[astronomy](../../index.md)/[io.github.cosinekitty.astronomy](index.md)/[massProduct](mass-product.md)

# massProduct

fun [massProduct](mass-product.md)(body: [Body](-body/index.md)): Double

Returns the product of mass and universal gravitational constant of a Solar System body.

For problems involving the gravitational interactions of Solar System bodies, it is helpful to know the product GM, where G = the universal gravitational constant and M = the mass of the body. In practice, GM is known to a higher precision than either G or M alone, and thus using the product results in the most accurate results. This function returns the product GM in the units au^3/day^2. The values come from page 10 of a [JPL memorandum regarding the DE405/LE405 ephemeris](https://web.archive.org/web/20120220062549/http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf).

#### Return

The mass product of the given body in au^3/day^2.

## Parameters

| | |
|---|---|
| body | The body for which to find the GM product. Allowed to be the Sun, Moon, EMB (Earth/Moon Barycenter), or any planet. Any other value will cause an exception to be thrown. |
