package io.github.cosinekitty.astronomy

/*
    Astronomy Engine for Kotlin / JVM.
    https://github.com/cosinekitty/astronomy

    MIT License

    Copyright (c) 2019-2022 Don Cross <cosinekitty@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/

/**
 * The enumeration of celestial bodies supported by Astronomy Engine.
 */
enum class Body
{
    /**
     * A placeholder value representing an invalid or unknown celestial body.
     */
    Invalid,

    /**
     * The planet Mercury.
     */
    Mercury,

    /**
     * The planet Venus.
     */
    Venus,

    /**
     * The planet Earth.
     * Some functions that accept a `Body` parameter will fail if passed this value
     * because they assume that an observation is being made from the Earth,
     * and therefore the Earth is not a target of observation.
     */
    Earth,

    /**
     * The planet Mars.
     */
    Mars,

    /**
     * The planet Jupiter.
     */
    Jupiter,

    /**
     * The planet Saturn.
     */
    Saturn,

    /**
     * The planet Uranus.
     */
    Uranus,

    /**
     * The planet Neptune.
     */
    Neptune,

    /**
     * The planet Pluto.
     */
    Pluto,

    /**
     * The Sun.
     */
    Sun,

    /**
     * The Earth's natural satellite, the Moon.
     */
    Moon,

    /**
     * The Earth/Moon Barycenter.
     */
    EMB,

    /**
     * The Solar System Barycenter.
     */
    SSB,
};
