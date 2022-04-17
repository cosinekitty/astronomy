package io.github.cosinekitty.astronomy

import java.io.File
import java.io.PrintWriter
import org.junit.jupiter.api.Test
import org.junit.jupiter.params.ParameterizedTest
import org.junit.jupiter.params.provider.CsvSource
import kotlin.test.fail
import kotlin.test.assertEquals
import kotlin.test.assertTrue
import kotlin.math.abs
import kotlin.math.cos
import kotlin.math.PI
import kotlin.math.sqrt
import kotlin.text.Regex

private const val dataRootDir = "../../generate/"
private const val MINUTES_PER_DAY = 24.0 * 60.0

private val regexTokenize = Regex("""\s+""")
private fun tokenize(s: String): List<String> = s.split(regexTokenize)

private fun tokenize(s: String, expectedTokenCount: Int, filename: String, lnum: Int): List<String> {
    val token = tokenize(s)
    assertEquals(expectedTokenCount, token.size, "$filename line $lnum: expected $expectedTokenCount tokens, found ${token.size}")
    return token
}

private fun stripLine(line: String): String {
    // Ignore anything after a '#' character and trim leading and trailing whitespace.
    val index = line.indexOf('#')
    return (if (index < 0) line else line.substring(0, index)).trim()
}

private fun groupString(match: MatchResult, groupIndex: Int, filename: String, lnum: Int): String =
    (match.groups[groupIndex] ?: fail("$filename line $lnum: cannot extract group $groupIndex")).value

private fun groupInt(match: MatchResult, groupIndex: Int, filename: String, lnum: Int): Int =
    groupString(match, groupIndex, filename, lnum).toInt()

private fun groupDouble(match: MatchResult, groupIndex: Int, filename: String, lnum: Int): Double =
    groupString(match, groupIndex, filename, lnum).toDouble()

private fun groupBody(match: MatchResult, groupIndex: Int, filename: String, lnum: Int): Body =
    Body.valueOf(groupString(match, groupIndex, filename, lnum))

private fun groupDirection(match: MatchResult, groupIndex: Int, filename: String, lnum: Int): Direction =
    when (val text = groupString(match, groupIndex, filename, lnum)) {
        "r" -> Direction.Rise
        "s" -> Direction.Set
        else -> fail("Invalid direction symbol: $text")
    }

private val regexDate = Regex("""^(\d+)-(\d+)-(\d+)T(\d+):(\d+)(:(\d+))?Z$""")

private val epoch2000 = Time(0.0)

private fun parseDate(text: String): Time {
    val m = regexDate.matchEntire(text) ?: fail("parseDate failed for string: '$text'")
    val year   = (m.groups[1] ?: fail("Cannot parse year from string: '$text'")).value.toInt()
    val month  = (m.groups[2] ?: fail("Cannot parse month from string: '$text'")).value.toInt()
    val day    = (m.groups[3] ?: fail("Cannot parse day from string: '$text'")).value.toInt()
    val hour   = (m.groups[4] ?: fail("Cannot parse hour from string: '$text'")).value.toInt()
    val minute = (m.groups[5] ?: fail("Cannot parse minute from string: '$text'")).value.toInt()
    val secondText = m.groups[7]?.value
    var second = (
        if (secondText != null && secondText != "")
            secondText.toDouble()
        else
            0.0
    )
    return Time(year, month, day, hour, minute, second)
}

private fun optionalParseDate(text: String): Time? = (
    if (text == "-")
        null
    else
        parseDate(text)
)

class Tests {
    private fun checkVector(
        vec: Vector,
        x: Double, y: Double, z: Double,
        tolerance: Double,
        message: String
    ) {
        assertTrue(vec.x.isFinite(), "$message: x is not a finite value")
        assertTrue(vec.y.isFinite(), "$message: y is not a finite value")
        assertTrue(vec.z.isFinite(), "$message: z is not a finite value")
        val dx = vec.x - x
        val dy = vec.y - y
        val dz = vec.z - z
        val diff = sqrt(dx*dx + dy*dy + dz*dz)
        assertTrue(
            diff < tolerance,
            "$message: EXCESSIVE ERROR: calc=(${vec.x}, ${vec.y}, ${vec.z}), expected=($x, $y, $z), diff=$diff"
        )
    }

    private fun checkScalar(
        calculated: Double,
        expected: Double,
        tolerance: Double,
        message: String
    ) {
        assertTrue(calculated.isFinite(), "$message: scalar is not finite")
        val diff = abs(calculated - expected)
        assertTrue(diff < tolerance, "$message: EXCESSIVE ERROR: calc=$calculated, expected=$expected, diff=$diff")
    }

    @Test
    fun `test deltaT calculation`() {
        val time = Time(0.0)
        assertEquals(0.0007389709440951036, time.tt)
    }

    @ParameterizedTest
    @CsvSource(
        /*
            These test data were made using the Python version of Astronmy Engine like this:

            don@spearmint:~/github/astronomy/demo/python $ python3
            Python 3.7.3 (default, Jan 22 2021, 20:04:44)
            [GCC 8.3.0] on linux
            Type "help", "copyright", "credits" or "license" for more information.
            >>> from astronomy import Time
            >>> t = Time.Make(1970, 12, 13, 23, 45, 12.345)
            >>> t
            Time('1970-12-13T23:45:12.345Z')
            >>> t.ut
            -10610.510273784721
        */
        value = [
            "2000,  1,  1, 12,  0,  0.0,        0.0,            '2000-01-01T12:00:00.000Z'",
            "2022,  1,  1, 12,  0,  0.0,     8036.0,            '2022-01-01T12:00:00.000Z'",
            "2022,  1,  1, 18,  0,  0.0,     8036.25,           '2022-01-01T18:00:00.000Z'",
            "1970, 12, 13, 23, 45, 12.345, -10610.510273784721, '1970-12-13T23:45:12.344Z'",    // we tolerate the 1 ms error in string
            "2022,  1,  1, 18, 59, 59.9999, 8036.291666665509,  '2022-01-01T18:59:59.999Z'",
        ]
    )
    fun `universal time calculation should match expectations`(
        year: Int, month: Int, day: Int, hour: Int, minute: Int, second: Double, ut: Double, expectedToString: String
    ) {
        val time = Time(year, month, day, hour, minute, second)
        assertEquals(ut, time.ut)
        assertEquals(expectedToString, time.toString())
    }

    @Test
    fun `Time should be able to add fractional days`() {
        val time = Time(2000, 1, 1, 12, 0, 0.0)
        assertEquals("2000-01-02T18:00:00.000Z", time.addDays(1.25).toString())
    }

    @Test
    fun `TerseVector methods should work as expected`() {
        val ones = TerseVector(1.0, 1.0, 1.0)
        assertEquals(ones, ones + TerseVector.zero())
        assertEquals(ones, ones + TerseVector.zero())
        assertEquals(TerseVector(6.0, 8.0, 4.0), TerseVector(3.0, 4.0, 2.0) * 2.0)
        assertEquals(TerseVector(-1.5, 2.0, -1.0), TerseVector(-3.0, 4.0, -2.0) / 2.0)
        assertEquals(29.0, TerseVector(-3.0, 4.0, -2.0).quadrature())
        assertEquals(5.744562646538029, TerseVector(-2.0, -2.0, 5.0).magnitude())
    }

    private fun testTimeRoundTrip(millis: Long, expectedFormat: String) {
        val time = Time.fromMillisecondsSince1970(millis)
        assertEquals(millis, time.toMillisecondsSince1970())
        assertEquals(expectedFormat, time.toString())
    }

    @Test
    fun `Convert milliseconds since 1970 to Time`() {
        // I created this test data using Node.js, like this:
        // $ node
        // Welcome to Node.js v14.19.1.
        // Type ".help" for more information.
        // > date = new Date('1066-07-15T11:54:37.123Z')
        // 1066-07-15T11:54:37.123Z
        // > date.getTime()
        // -28510574722877
        testTimeRoundTrip(-28510574722877L, "1066-07-15T11:54:37.123Z")

        // The answer should end with ".762Z", but we are one millsecond off.
        // We tolerate 1ms round-trip formatting error, due to the limitations of floating point math.
        testTimeRoundTrip(1649942537762L, "2022-04-14T13:22:17.761Z")
    }

    //----------------------------------------------------------------------------------------
    // RotationMatrix tests

    private fun compareMatrices(testName: String, a: RotationMatrix, b: RotationMatrix, tolerance: Double = 1.0e-15) {
        for (i: Int in 0..2) {
            for (j: Int in 0..2) {
                assertTrue(abs(a.rot[i][j] - b.rot[i][j]) < tolerance, "Matrix mismatch at i=$i, j=$j: $testName")
            }
        }
    }

    private fun compareVectors(testName: String, a: Vector, b: Vector, tolerance: Double = 1.0e-15) {
        assertTrue(a.x.isFinite(), "a.x is not finite")
        assertTrue(a.y.isFinite(), "a.y is not finite")
        assertTrue(a.z.isFinite(), "a.z is not finite")

        assertTrue(b.x.isFinite(), "b.x is not finite")
        assertTrue(b.y.isFinite(), "b.y is not finite")
        assertTrue(b.z.isFinite(), "b.z is not finite")

        assertTrue(abs(a.x - b.x) < tolerance, "a.x=${a.x} but b.x=${b.x}: $testName")
        assertTrue(abs(a.y - b.y) < tolerance, "a.x=${a.x} but b.x=${b.x}: $testName")
        assertTrue(abs(a.z - b.z) < tolerance, "a.x=${a.x} but b.x=${b.x}: $testName")

        assertTrue(a.t.tt == b.t.tt, "a.t=${a.t}, but b.t=${b.t}")
    }

    @Test
    fun `Pivot rotation matrix`() {
        // Pivot an identity matrix 90 degrees counterclockwise around the z-axis.
        var r = RotationMatrix.identity().pivot(2, +90.0)

        compareMatrices("Pivot #1", r, RotationMatrix(
             0.0, +1.0,  0.0,
            -1.0,  0.0,  0.0,
             0.0,  0.0, +1.0
        ))

        // Pivot again, -30 degrees around the x-axis.
        r = r.pivot(0, -30.0)

        // Pivot a third time, 180 degrees around the y-axis.
        r = r.pivot(1, +180.0)

        // Use the resulting matrix to rotate a vector.
        val time = Time(2000, 1, 1, 0, 0, 0.0)
        val v1 = Vector(1.0, 2.0, 3.0, time)

        val v2 = r.rotate(v1)

        compareVectors("Pivot #2", v2, Vector(
            +2.0,
            +2.3660254037844390,
            -2.0980762113533156,
            time
        ))
    }

    @Test
    fun `Combine rotation matrices`() {
        // I derived these calculation answers by running the same steps
        // in the Python version of Astronomy Engine, which has already been validated.

        val a = RotationMatrix
            .identity()
            .pivot(0, +37.0)
            .pivot(1, -72.0)
            .pivot(2, +23.0)

        compareMatrices("matrix a", a, RotationMatrix(
            +0.28445164312142457, +0.12074255893448678, +0.9510565162951535,
            -0.8389120034878212,  +0.5115089556077241,  +0.18597106962413537,
            -0.46401930253985274, -0.8507525038429493,  +0.24679194491591758
        ))

        val b = RotationMatrix
            .identity()
            .pivot(0, -94.0)
            .pivot(1, +63.0)
            .pivot(2, +112.0)

        compareMatrices("matrix b", b, RotationMatrix(
            -0.17006783455061913, +0.42093266148521513, -0.8910065241883678,
            +0.3976409311461465,  -0.7979832250245378,  -0.452884601699664,
            -0.9016421804288488,  -0.4313217674479356,  -0.031668776375164034
        ))

        val c = a combine b

        compareMatrices("c = a combine b", c, RotationMatrix(
            -0.8578765624797756,  -0.38682692692461523, -0.33824951168322753,
            +0.17838948449655287, -0.8415143988492727,  +0.5099320624851914,
            -0.4818972871566489,  +0.3771186088426366,  +0.7909213358455168
        ))
    }

    //----------------------------------------------------------------------------------------

    @Test
    fun `Sanity check sidereal time`() {
        // An exhaustive test is not needed, because sidereal time calculations
        // are indirectly tested through their use in thousands of other
        // verified calculations. This is just to help isolate a problem
        // in sidereal time in case it is broken.
        val correct = 9.398368460418821
        val time = Time(2022, 3, 15, 21, 50, 0.0)
        val gast = siderealTime(time)
        assertTrue(gast.isFinite())
        val diff = abs(gast - correct)
        assertTrue(diff < 1.0e-15, "correct=$correct, gast=$gast, diff=$diff")
    }

    //----------------------------------------------------------------------------------------

    @Test
    fun `Verify body axis calculations`() {
        axisTestBody(Body.Sun,      "Sun.txt",       0.0)
        axisTestBody(Body.Mercury,  "Mercury.txt",   0.074340)
        axisTestBody(Body.Venus,    "Venus.txt",     0.0)
        axisTestBody(Body.Earth,    "Earth.txt",     0.000591)
        axisTestBody(Body.Moon,     "Moon.txt",      0.264845)
        axisTestBody(Body.Mars,     "Mars.txt",      0.075323)
        axisTestBody(Body.Jupiter,  "Jupiter.txt",   0.000324)
        axisTestBody(Body.Saturn,   "Saturn.txt",    0.000304)
        axisTestBody(Body.Uranus,   "Uranus.txt",    0.0)
        axisTestBody(Body.Neptune,  "Neptune.txt",   0.000462)
        axisTestBody(Body.Pluto,    "Pluto.txt",     0.0)
    }

    private fun axisTestBody(body: Body, shortName: String, arcminTolerance: Double) {
        val filename = dataRootDir + "axis/" + shortName
        val file = File(filename)
        var lnum = 0
        var foundData = false
        var maxArcmin = 0.0
        for (line in file.readLines()) {
            ++lnum
            if (!foundData) {
                if (line == "\$\$SEO")
                    foundData = true
            } else {
                if (line == "\$\$EOE")
                    break

                assertTrue(line.length >= 61, "axisBodyTest($filename line $lnum): line is too short")
                val token = tokenize(line.substring(19), 3, filename, lnum)
                val jd = token[0].toDouble()
                val ra = token[1].toDouble()
                val dec = token[2].toDouble()
                val time = Time(jd - 2451545.0)
                val axis = rotationAxis(body, time)

                // Convert the reference angles to a reference north pole vector.
                // tricky: `ra` is in degrees, not sidereal hours; so don't multiply by 15.
                val sphere = Spherical(dec, ra, 1.0);
                val north = sphere.toVector(time)

                // Find angle between two versions of the north pole. Use that as the measure of error.
                val arcmin = 60.0 * north.angleWith(axis.north)
                if (arcmin > maxArcmin)
                    maxArcmin = arcmin
            }
        }
        assertTrue(maxArcmin <= arcminTolerance, "axisBodyTest($body): excessive error = $maxArcmin arcmin.")
        assertTrue(lnum > 900, "axisBodyTest($filename): only processed $lnum lines")
    }

    //----------------------------------------------------------------------------------------

    @Test
    fun `Sanity check geocentric moon`() {
        val time = Time(2019, 6, 24, 15, 45, 37.0)

        val eclSphere = eclipticGeoMoon(time)
        val dlat = eclSphere.lat  - (-4.851798346972171)
        val dlon = eclSphere.lon  - (+354.5951298193645)
        var drad = eclSphere.dist - 0.0026968810499258147
        assertTrue(abs(dlat) < 1.0e-15, "eclipticGeoMoon: excessive latitude error $dlat")
        assertTrue(abs(dlon) < 1.0e-15, "eclipticGeoMoon: excessive longitude error $dlon")
        assertTrue(abs(drad) < 1.0e-17, "eclipticGeoMoon: excessive distance error $drad")

        val equVec = geoMoon(time)
        checkVector(
            equVec,
            +0.002674037026701135, -0.0001531610316600666, -0.0003150159927069429,
            5.43e-20,
            "geoMoon"
        )
    }

    //----------------------------------------------------------------------------------------

    @Test
    fun `Heliocentric vectors and distances`() {
        val time = Time(2022, 3, 28, 15, 21, 41.0)

        verifyHelio(
            Body.Sun,
            time,
            0.0, 0.0, 0.0,
            0.0, 0.0, 0.0
        )

        // I used the Python version to calculate reference state vectors.
        // This is just a sanity check.
        // Exhaustive tests referenced against JPL Horizons data
        // are done separately by the test `Heliocentric state vectors`.

        verifyHelio(
            Body.Mercury,
            time,
            0.3601225052845418, -0.05702761412960074, -0.06778598705660548,
            0.0005993053269706367, 0.025455638276077892, 0.013536525517872652
        )

        verifyHelio(
            Body.Venus,
            time,
            -0.41031775375890733, -0.553797703820094, -0.2232212930914451,
            0.01652666689705999, -0.010155161581682825, -0.005615571915795434
        )

        verifyHelio(
            Body.Earth,
            time,
            -0.989323089873633, -0.12148064775730126, -0.05265719320840268,
            0.001998233870734938, -0.01571097012164241, -0.00681060160402074
        )

        verifyHelio(
            Body.Mars,
            time,
            0.33465883914714956, -1.2596467270049816, -0.5868032433717035,
            0.014131155779203354, 0.0042119182043996745, 0.0015503545422968985
        )

        verifyHelio(
            Body.Jupiter,
            time,
            4.843202980566162, -1.0038257126069707, -0.5481395157798546,
            0.0016395952415457383, 0.007099637074264691, 0.0030031754891337277
        )

        verifyHelio(
            Body.Saturn,
            time,
            7.2700610693589, -6.096752520060888, -2.8318384086244084,
            0.0034786432886192487, 0.003838295762587206, 0.00143581155860176
        )

        verifyHelio(
            Body.Uranus,
            time,
            14.159758260480487, 12.631575584016622, 5.332054883582301,
            -0.002763132653154497, 0.0024129846977416014, 0.0010960328116071437
        )

        verifyHelio(
            Body.Neptune,
            time,
            29.66860888675136, -3.261664651496572, -2.0721405654341702,
            0.00038196533635370553, 0.0029105280185167605, 0.0011819337466936254
        )

        verifyHelio(
            Body.Pluto,
            time,
            15.377665594383952, -27.85223298004125, -13.32288901996256,
            0.0028887491551056084, 0.0010313534744949828, -0.0005469417328328084
        )

        verifyHelio(
            Body.Moon,
            time,
            -0.9873532819712121, -0.12279866865704744, -0.05347300317660862,
            0.002381799204808397, -0.015280061848635642, -0.006626397121575363
        )

        verifyHelio(
            Body.EMB,
            time,
            -0.9892991555540802, -0.1214966624830789, -0.05266710577726255,
            0.0020028944141634816, -0.015705734333781345, -0.006808363411687111
        )

        verifyHelio(
            Body.SSB,
            time,
            0.008844144038360908, -0.002316519421832873, -0.0012061496455500675,
            2.4572258206466927e-06, 8.124488621688556e-06, 3.38376990631319e-06
        )
    }

    @Test
    fun `Another Pluto test`() {
        // This time offset caused an excessive discrepancy compared to the C code's output.
        val time = Time.fromTerrestrialTime(-58402.247295546535)
        verifyHelio(
            Body.Pluto,
            time,
            42.937184148112564, 20.340210822720223, -6.589767515046609,
            -0.0005050118478242109, 0.0019911149712207757, 0.0007743560643446314
        )
    }

    @Test
    fun `Venus horizontal coords`() {
        val body = Body.Venus
        val time = Time.fromTerrestrialTime(10373.141119633277)  // UT 2028-05-26T15:21:56.183Z
        val pos = helioVector(body, time)
        checkVector(
            pos,
            -0.34191041083994594, -0.5908265808125958, -0.24422546639998666,
            1.0e-16,
            "Venus heliocentric position"
        )

        val observer = Observer(29.0, -81.0, 10.0)

        val ofdate = equator(body, time, observer, EquatorEpoch.OfDate, Aberration.Corrected)
        checkScalar(ofdate.ra,   4.893134408107621,   8.9e-16, "Venus EQD RA")
        checkScalar(ofdate.dec, 24.6998405830952,     1.0e-16, "Venus EQD DEC")
        checkScalar(ofdate.dist, 0.29355004763155124, 1.0e-16, "Venus EQD distance")

        val hor = horizon(time, observer, ofdate.ra, ofdate.dec, Refraction.None)
        checkScalar(hor.azimuth,  87.5963687042015,   1.0e-16, "Venus azimuth")
        checkScalar(hor.altitude, 54.929061963517746, 7.2e-15, "Venus altitude")
    }

    private fun verifyHelio(
        body: Body,
        time: Time,
        x: Double, y: Double, z: Double,
        vx: Double, vy: Double, vz: Double
    ) {
        val tolerance = 1.0e-11
        val calcPos = helioVector(body, time)
        checkVector(calcPos, x, y, z, tolerance, "helioVector($body)")

        val calcDist = helioDistance(body, time)
        val expectedDist = sqrt(x*x + y*y + z*z)
        checkScalar(calcDist, expectedDist, tolerance, "helioDistance($body)")

        val calcState = helioState(body, time)
        checkVector(calcState.position(), x, y, z, tolerance, "helioState($body).position")
        checkVector(calcState.velocity(), vx, vy, vz, tolerance, "helioState($body).velocity")
    }

    //----------------------------------------------------------------------------------------

    @Test
    fun `Write numeric check output file`() {
        // The primary purpose of this test is to generate an output
        // file full of various calculations in a format that is understood
        // by the `generate check` command, which runs later in the build process.
        // That command compares the output of one Astronomy Engine language
        // implementation with another.
        // Even if this test passes, the numeric output may deviate from expectations
        // in a later build step, causing the script `unit_test_kotlin` to fail.
        // For comparison, see the functions named `AstroCheck` in:
        // - generate/ctest.c
        // - generate/test.js
        // - generate/test.py
        // - generate/dotnet/csharp_test/csharp_test.cs
        val filename = dataRootDir + "temp/k_check.txt"
        File(filename).printWriter().use { outfile -> AstroCheck(outfile) }
    }

    private fun AstroCheck(outfile: PrintWriter) {
        val bodylist = arrayOf (
            Body.Sun, Body.Mercury, Body.Venus, Body.Earth, Body.Mars,
            Body.Jupiter, Body.Saturn, Body.Uranus, Body.Neptune, Body.Pluto,
            Body.SSB, Body.EMB
        )

        val vectorOnlyList = arrayOf(
            Body.Earth, Body.SSB, Body.EMB
        )

        val observer = Observer(29.0, -81.0, 10.0)
        var time = Time(1700, 1, 1, 0, 0, 0.0)
        val stop = Time(2200, 1, 1, 0, 0, 0.0)
        var pos: Vector
        var j2000: Equatorial
        var ofdate: Equatorial
        var hor: Topocentric

        outfile.println("o ${observer.latitude} ${observer.longitude} ${observer.height}")
        while (time.tt < stop.tt) {
            for (body in bodylist) {
                pos = helioVector(body, time)
                outfile.println("v ${body} ${pos.t.tt} ${pos.x} ${pos.y} ${pos.z}")
                if (body !in vectorOnlyList) {
                    j2000 = equator(body, time, observer, EquatorEpoch.J2000, Aberration.None)
                    ofdate = equator(body, time, observer, EquatorEpoch.OfDate, Aberration.Corrected)
                    hor = horizon(time, observer, ofdate.ra, ofdate.dec, Refraction.None)
                    outfile.println("s ${body} ${time.tt} ${time.ut} ${j2000.ra} ${j2000.dec} ${j2000.dist} ${hor.azimuth} ${hor.altitude}")
                }
            }
            pos = geoVector(Body.Moon, time, Aberration.None)
            outfile.println("v GM ${pos.t.tt} ${pos.x} ${pos.y} ${pos.z}")
            j2000 = equator(Body.Moon, time, observer, EquatorEpoch.J2000, Aberration.None)
            ofdate = equator(Body.Moon, time, observer, EquatorEpoch.OfDate, Aberration.Corrected)
            hor = horizon(time, observer, ofdate.ra, ofdate.dec, Refraction.None)
            outfile.println("s GM ${time.tt} ${time.ut} ${j2000.ra} ${j2000.dec} ${j2000.dist} ${hor.azimuth} ${hor.altitude}")

            val jm = jupiterMoons(time)
            var mindex = 0
            for (moon in jm.moon) {
                outfile.println("j ${mindex} ${time.tt} ${time.ut} ${moon.x} ${moon.y} ${moon.z} ${moon.vx} ${moon.vy} ${moon.vz}")
                ++mindex
            }

            time = time.addDays(10.0 + PI/100.0)
        }
    }

    //----------------------------------------------------------------------------------------

    private fun compareMatrices(a: RotationMatrix, b: RotationMatrix, comment: String, tolerance: Double = 1.0e-99) {
        for (i in 0..2) {
            for (j in 0..2) {
                val diff = abs(a.rot[i][j] - b.rot[i][j])
                assertTrue(diff < tolerance, "Matrix mismatch ($comment): a[$i,$j]=${a.rot[i][j]}, b=${b.rot[i][j]}, diff=$diff")
            }
        }
    }

    @Test
    fun `Orientation conversion using rotation matrices`() {
        val time = Time(8126.418466077072)
        assertEquals("2022-04-01T22:02:35.469Z", time.toString(), "Unexpected time string.")

        val observer = Observer(-30.0, +150.0, 200.0)

        compareMatrices(
            rotationEqjEcl(),
            RotationMatrix(
                1.0, 0.0, 0.0,
                0.0, +0.9174821430670688, -0.3977769691083922,
                0.0, +0.3977769691083922, +0.9174821430670688
            ),
            "EQJ ECL"
        )

        compareMatrices(
            rotationEclEqj(),
            RotationMatrix(
                1.0, 0.0, 0.0,
                0.0, +0.9174821430670688, +0.3977769691083922,
                0.0, -0.3977769691083922, +0.9174821430670688
            ),
            "EQJ ECL"
        )

        compareMatrices(
            rotationEqjEqd(time),
            RotationMatrix(
                0.9999856608656787, 0.004911515527973243, 0.002134262929010771,
                -0.004911577234051216, 0.9999879378516134, 2.367175135753887e-05,
                -0.002134120921040258, -3.4154009138639524e-05, 0.9999977221781049
            ),
            "EQJ EQD"
        )

        compareMatrices(
            rotationEqdEqj(time),
            RotationMatrix(
                0.9999856608656787, -0.004911577234051216, -0.002134120921040258,
                0.004911515527973243, 0.9999879378516134, -3.4154009138639524e-05,
                0.002134262929010771, 2.367175135753887e-05, 0.9999977221781049,
            ),
            "EQD EQJ"
        )

        compareMatrices(
            rotationEqdHor(time, observer),
            RotationMatrix(
                0.3272894142412824, -0.7559937548038297, 0.5668818942453582,
                -0.3779968774019148, -0.6545788284825649, -0.6547097967625005,
                0.8660254037844387, 0.0, -0.49999999999999994,
            ),
            "EQD HOR"
        )

        compareMatrices(
            rotationHorEqd(time, observer),
            RotationMatrix(
                0.3272894142412824, -0.3779968774019148, 0.8660254037844387,
                -0.7559937548038297, -0.6545788284825649, 0.0,
                0.5668818942453582, -0.6547097967625005, -0.49999999999999994,
            ),
            "HOR EQD"
        )

        compareMatrices(
            rotationHorEqj(time, observer),
            RotationMatrix(
                0.3272765095764035, -0.37957932484539564, 0.86533786605545,
                -0.7591978885882081, -0.6508578111404256, 0.0016357385795925856,
                0.5625910168521117, -0.657498019637632, -0.5011862946349386,
            ),
            "HOR EQJ"
        )

        compareMatrices(
            rotationEqjHor(time, observer),
            RotationMatrix(
                0.3272765095764035, -0.7591978885882081, 0.5625910168521117,
                -0.37957932484539564, -0.6508578111404256, -0.657498019637632,
                0.86533786605545, 0.0016357385795925856, -0.5011862946349386,
            ),
            "EQJ HOR"
        )

        compareMatrices(
            rotationEqdEcl(time),
            RotationMatrix(
                0.9999856608656787, -0.00535518855821894, -4.305530497609837e-06,
                0.004911515527973243, 0.917457490583079, -0.39780350675706494,
                0.002134262929010771, 0.39779778145246825, 0.9174706371286464,
            ),
            "EQD ECL"
        )

        compareMatrices(
            rotationEclEqd(time),
            RotationMatrix(
                0.9999856608656787, 0.004911515527973243, 0.002134262929010771,
                -0.00535518855821894, 0.917457490583079, 0.39779778145246825,
                -4.305530497609837e-06, -0.39780350675706494, 0.9174706371286464,
            ),
            "ECL EQD"
        )

        compareMatrices(
            rotationEclHor(time, observer),
            RotationMatrix(
                0.3272765095764035, -0.7591978885882081, 0.5625910168521117,
                -0.004045778808843881, -0.5964997602626152, -0.8026030573580397,
                0.9449199531988498, 0.26039700837346297, -0.1982919062312794,
            ),
            "ECL HOR"
        )

        compareMatrices(
            rotationHorEcl(time, observer),
            RotationMatrix(
                0.3272765095764035, -0.004045778808843881, 0.9449199531988498,
                -0.7591978885882081, -0.5964997602626152, 0.26039700837346297,
                0.5625910168521117, -0.8026030573580397, -0.1982919062312794,
            ),
            "HOR ECL"
        )

        compareMatrices(
            rotationGalEqj(),
            RotationMatrix(
                -0.0548624779711344, -0.8734572784246782, -0.483800052994852,
                0.4941095946388765, -0.4447938112296831, 0.7470034631630423,
                -0.8676668813529025, -0.1980677870294097, 0.4559861124470794,
            ),
            "GAL EQJ"
        )

        compareMatrices(
            rotationEqjGal(),
            RotationMatrix(
                -0.0548624779711344, 0.4941095946388765, -0.8676668813529025,
                -0.8734572784246782, -0.4447938112296831, -0.1980677870294097,
                -0.483800052994852, 0.7470034631630423, 0.4559861124470794,
            ),
            "EQJ GAL"
        )
    }

    //----------------------------------------------------------------------------------------

    @Test
    fun `Verify constellations`() {
        val filename = dataRootDir + "constellation/test_input.txt"
        val infile = File(filename)
        var lnum = 0
        val regex = Regex("""^\s*(\d+)\s+(\S+)\s+(\S+)\s+([A-Z][a-zA-Z]{2})\s*$""")
        for (line in infile.readLines()) {
            val m = regex.matchEntire(line) ?: fail("$filename line $lnum: syntax error")
            val id: Int = groupInt(m, 1, filename, lnum)
            val ra: Double = groupDouble(m, 2, filename, lnum)
            val dec: Double = groupDouble(m, 3, filename, lnum)
            val symbol: String = groupString(m, 4, filename, lnum)
            val constel = constellation(ra, dec)
            assertTrue(constel.symbol == symbol, "$filename line $lnum: expected constellation $symbol, but found id=$id, symbol=${constel.symbol}")

            ++lnum;
        }
    }

    //----------------------------------------------------------------------------------------

    @Test
    fun `Basic search`() {
        val toleranceDays = 1.0e-9
        val toleranceSeconds = toleranceDays / 86400.0
        val time1 = Time(0.0)
        val time2 = Time(1.0)
        val tsolve = search(time1, time2, toleranceSeconds) { time -> time.ut - cos(time.ut) }
        if (tsolve == null)
            fail("Basic search failed")

        val utCorrect = 0.7390851332151607
        val diff = abs(tsolve.ut - utCorrect)
        assertTrue(diff <= toleranceDays, "search found a solution outside tolerance: diff = $diff")
    }

    //----------------------------------------------------------------------------------------

    @Test
    fun `Seasons test`() {
        val filename = dataRootDir + "seasons/seasons.txt"
        val regex = Regex("""^(\d+)-(\d+)-(\d+)T(\d+):(\d+)Z\s+([A-Za-z]+)\s*$""")
        var lnum = 0
        var currentYear = 0
        var marchCount = 0
        var juneCount = 0
        var septemberCount = 0
        var decemberCount = 0
        var seasons: SeasonsInfo? = null
        val infile = File(filename)
        var maxMinutes = 0.0
        for (line in infile.readLines()) {
            ++lnum
            // 2019-01-03T05:20Z Perihelion
            // 2019-03-20T21:58Z Equinox
            // 2019-06-21T15:54Z Solstice
            // 2019-07-04T22:11Z Aphelion
            // 2019-09-23T07:50Z Equinox
            // 2019-12-22T04:19Z Solstice
            val m = regex.matchEntire(line) ?: fail("$filename line $lnum: syntax error")
            val year = groupInt(m, 1, filename, lnum)
            val month = groupInt(m, 2, filename, lnum)
            val day = groupInt(m, 3, filename, lnum)
            val hour = groupInt(m, 4, filename, lnum)
            val minute = groupInt(m, 5, filename, lnum)
            val name = groupString(m, 6, filename, lnum)
            val correctTime = Time(year, month, day, hour, minute, 0.0)
            if (year != currentYear) {
                currentYear = year
                seasons = seasons(year)
            }

            if (seasons == null)
                fail("internal error: seasons == null")

            var calcTime: Time
            when (name) {
                "Equinox" -> when (month) {
                    3 -> {
                        calcTime = seasons.marchEquinox
                        ++marchCount
                    }
                    9 -> {
                        calcTime = seasons.septemberEquinox
                        ++septemberCount
                    }
                    else -> fail("$filename line $lnum: Invalid equinox date in test data.")
                }
                "Solstice" -> when (month) {
                    6 -> {
                        calcTime = seasons.juneSolstice
                        ++juneCount
                    }
                    12 -> {
                        calcTime = seasons.decemberSolstice
                        ++decemberCount
                    }
                    else -> fail("$filename line $lnum: Invalid solstice date in test data.")
                }
                "Aphelion", "Perihelion" -> {
                    // not yet calculated
                    continue
                }
                else -> fail("$filename line $lnum: unknown event type: $name")
            }
            // Verify that the calculated time matches the current time for this event.
            val diffMinutes = (24.0 * 60.0) * abs(calcTime.tt - correctTime.tt)
            if (diffMinutes > maxMinutes)
                maxMinutes = diffMinutes

            if (diffMinutes > 2.38)
                fail("$filename line $lnum: excessive error ($name): $diffMinutes minutes.")
        }
        assertTrue(marchCount == 301, "marchCount = $marchCount")
        assertTrue(juneCount == 301, "juneCount = $juneCount")
        assertTrue(septemberCount == 301, "septemberCount = $septemberCount")
        assertTrue(decemberCount == 301, "decemberCount = $decemberCount")
    }

    @Test
    fun `Seasons issue 187`() {
        // This is a regression test for:
        // https://github.com/cosinekitty/astronomy/issues/187
        // For years far from the present, the seasons search was sometimes failing.
        for (year in -2000 .. 9999)
            seasons(year)
    }

    //----------------------------------------------------------------------------------------

    @Test
    fun `Moon phase search`() {
        var thresholdSeconds = 90.0
        val filename = dataRootDir + "moonphase/moonphases.txt"
        val infile = File(filename)
        var lnum = 0
        var maxArcmin = 0.0
        var prevYear = 0
        var mq: MoonQuarterInfo? = null
        // 0 1800-01-25T03:21:00.000Z
        // 1 1800-02-01T20:40:00.000Z
        // 2 1800-02-09T17:26:00.000Z
        // 3 1800-02-16T15:49:00.000Z
        var re = Regex("""^([0-3])\s+(\d+)-(\d+)-(\d+)T(\d+):(\d+):(\d+)\.000Z$""")
        for (line in infile.readLines()) {
            ++lnum
            val m = re.matchEntire(line) ?: fail("$filename line $lnum: syntax error")
            val quarter = groupInt(m, 1, filename, lnum)
            val year = groupInt(m, 2, filename, lnum)
            val month = groupInt(m, 3, filename, lnum)
            val day = groupInt(m, 4, filename, lnum)
            val hour = groupInt(m, 5, filename, lnum)
            val minute = groupInt(m, 6, filename, lnum)
            val second = groupDouble(m, 7, filename, lnum)

            val expectedTime = Time(year, month, day, hour, minute, second)
            val expectedElong = 90.0 * quarter
            val calcElong = moonPhase(expectedTime)
            var degreeError = abs(calcElong - expectedElong)
            if (degreeError > 180.0)
                degreeError = 360.0 - degreeError
            val arcmin = 60.0 * degreeError
            assertTrue(arcmin < 1.0, "$filename line $lnum: excessive angular error = $arcmin arcmin.")
            if (arcmin > maxArcmin)
                maxArcmin = arcmin
            if (year != prevYear) {
                prevYear = year
                // The test data contains a single year's worth of data for every 10 years.
                // Every time we see the year value change, it breaks continuity of the phases.
                // Start the search over again.
                val startTime = Time(year, 1, 1, 0, 0, 0.0)
                mq = searchMoonQuarter(startTime)
            } else {
                if (mq == null)
                    fail("mq == null")  // should not be possible

                // Yet another lunar quarter in the same year.
                val expectedQuarter = (1 + mq.quarter) % 4
                mq = nextMoonQuarter(mq)

                // Make sure we find the next expected quarter.
                assertTrue(expectedQuarter == mq.quarter, "$filename line $lnum: expected quarter $expectedQuarter, but found ${mq.quarter}")
            }

            // Make sure the time matches what we expect.
            val diffSeconds = abs(mq.time.tt - expectedTime.tt) * 86400
            assertTrue(diffSeconds <= thresholdSeconds, "$filename line $lnum: excessive time error = $diffSeconds seconds.")
        }
    }

    //----------------------------------------------------------------------------------------

    @Test
    fun `Rise set test`() {
        val filename = dataRootDir + "riseset/riseset.txt"
        val infile = File(filename)
        var lnum = 0
        // Moon  103 -61 1944-01-02T17:08Z s
        // Moon  103 -61 1944-01-03T05:47Z r
        var re = Regex("""^([A-Za-z]+)\s+([\-\+]?\d+\.?\d*)\s+([\-\+]?\d+\.?\d*)\s+(\d+)-(\d+)-(\d+)T(\d+):(\d+)Z\s+([rs])\s*$""")
        var currentBody: Body? = null
        var observer: Observer? = null
        var rSearchDate: Time? = null
        var sSearchDate: Time? = null
        var rEvt: Time?            // rise event search result
        var sEvt: Time?            // set event search result
        var aEvt: Time?            // chronologically first event (whether rise or set)
        var bEvt: Time? = null     // chronologically second event (whether rise or set)
        var aDir: Direction
        var bDir = Direction.Rise
        val nudgeDays = 0.01
        for (line in infile.readLines()) {
            ++lnum
            var m: MatchResult = re.matchEntire(line) ?: fail("$filename line $lnum: syntax error")
            val body = groupBody(m, 1, filename, lnum)
            val longitude = groupDouble(m, 2, filename, lnum)
            val latitude = groupDouble(m, 3, filename, lnum)
            val year = groupInt(m, 4, filename, lnum)
            val month = groupInt(m, 5, filename, lnum)
            val day = groupInt(m, 6, filename, lnum)
            val hour = groupInt(m, 7, filename, lnum)
            val minute = groupInt(m, 8, filename, lnum)
            val direction = groupDirection(m, 9, filename, lnum)

            val correctDate = Time(year, month, day, hour, minute, 0.0)

            // Every time we see a new geographic location or body, start a new iteration
            // of finding all rise/set times for that UTC calendar year.
            if (observer == null || observer.latitude != latitude || observer.longitude != longitude
                || currentBody == null || currentBody != body
            ) {
                currentBody = body
                observer = Observer(latitude, longitude, 0.0)
                rSearchDate = Time(year, 1, 1, 0, 0, 0.0)
                sSearchDate = rSearchDate
                bEvt = null
            }

            if (bEvt != null) {
                // The previous iteration found two events.
                // We already processed the earlier event (aEvt).
                // Now it is time to process the later event (bEvt).
                aEvt = bEvt
                aDir = bDir
                bEvt = null
            } else {
                rEvt = searchRiseSet(body, observer, Direction.Rise, rSearchDate!!, 366.0)
                if (rEvt == null)
                    fail("$filename line $lnum: did not find $body rise event.")

                sEvt = searchRiseSet(body, observer, Direction.Set, sSearchDate!!, 366.0)
                if (sEvt == null)
                    fail("$filename line $lnum: did not find $body set event.")

                // Sort the two events chronologically.
                // We will check the earlier event in this iteration,
                // and check the later event in the next iteration.

                if (rEvt.tt < sEvt.tt) {
                    aEvt = rEvt
                    bEvt = sEvt
                    aDir = Direction.Rise
                    bDir = Direction.Set
                } else {
                    aEvt = sEvt
                    bEvt = rEvt
                    aDir = Direction.Set
                    bDir = Direction.Rise
                }

                // Nudge the event times forward a tiny amount.
                // This prevents us from getting stuck in a loop, finding the same event repeatedly.
                rSearchDate = rEvt.addDays(nudgeDays)
                sSearchDate = sEvt.addDays(nudgeDays)
            }

            // Expect the current search result to match the earlier of the found dates.
            assertTrue(aDir == direction, "$filename line $lnum: expected direction $direction, bound found $aDir")
            val errorMinutes = (24.0 * 60.0) * abs(aEvt.tt - correctDate.tt)
            assertTrue(errorMinutes < 0.57, "$filename line $lnum: excessive prediction time error = $errorMinutes minutes.")
        }
    }

    //----------------------------------------------------------------------------------------

    @Test
    fun `Twilight test`() {
        val filename = dataRootDir + "riseset/twilight.txt"
        val toleranceSeconds = 60.0
        var lnum = 0
        val infile = File(filename)
        for (line in infile.readLines()) {
            ++lnum
            val tokens = tokenize(line, 9, filename, lnum)
            val lat = tokens[0].toDouble()
            val lon = tokens[1].toDouble()
            val observer = Observer(lat, lon, 0.0)
            val searchDate = parseDate(tokens[2])
            val correctTimes = tokens.drop(3).map { it -> parseDate(it) }
            val calcTimes = arrayOf(
                searchAltitude(Body.Sun, observer, Direction.Rise, searchDate, 1.0, -18.0),  // astronomical dawn
                searchAltitude(Body.Sun, observer, Direction.Rise, searchDate, 1.0, -12.0),  // nautical dawn
                searchAltitude(Body.Sun, observer, Direction.Rise, searchDate, 1.0,  -6.0),  // civil dawn
                searchAltitude(Body.Sun, observer, Direction.Set,  searchDate, 1.0,  -6.0),  // civil dawn
                searchAltitude(Body.Sun, observer, Direction.Set,  searchDate, 1.0, -12.0),  // nautical dawn
                searchAltitude(Body.Sun, observer, Direction.Set,  searchDate, 1.0, -18.0),  // astronomical dawn
            )
            assertTrue(correctTimes.size == calcTimes.size, "correctTimes.size = ${correctTimes.size}, but calcTimes.size = ${calcTimes.size}")
            for (i in 0 until correctTimes.size) {
                val correct = correctTimes[i]
                val calc = calcTimes[i] ?: fail("$filename line $lnum: calcTimes[$i] search failed.")
                val diff = 86400.0 * abs(calc.ut - correct.ut)
                assertTrue(diff < toleranceSeconds, "$filename line $lnum: excessive error = $diff seconds.")
            }
        }
    }

    //----------------------------------------------------------------------------------------

    @Test
    fun `Compare geoid state vectors with C`() {
        val filename = dataRootDir + "topostate/geoid.txt"
        val infile = File(filename)
        var lnum = 0
        val speedConvert = KM_PER_AU / 24.0     // converts AU/day into km/hour
        val posTolerance = 1.0e-6   // 1 millimeter
        val velTolerance = 1.0e-6   // 1 mm/hour
        for (line in infile.readLines()) {
            ++lnum
            // 0   -340866.750000000000    1.5     2.7     7.4  2233.390067 -5966.178081   311.806727  1559.005842   591.264250   146.628157
            val tokens = tokenize(line, 11, filename, lnum)
            val equator: EquatorEpoch = when (tokens[0]) {
                "0" -> EquatorEpoch.J2000
                "1" -> EquatorEpoch.OfDate
                else -> fail("$filename line $lnum: unexpected EquatorEpoch code: ${tokens[0]}")
            }
            val ut = tokens[1].toDouble()
            val latitude = tokens[2].toDouble()
            val longitude = tokens[3].toDouble()
            val height = tokens[4].toDouble()
            val x = tokens[5].toDouble()
            val y = tokens[6].toDouble()
            val z = tokens[7].toDouble()
            val vx = tokens[8].toDouble()
            val vy = tokens[9].toDouble()
            val vz = tokens[10].toDouble()

            val time = Time(ut)
            val observer = Observer(latitude, longitude, height)
            val correctState = StateVector(x, y, z, vx, vy, vz, time)
            val calcState = observer.toStateVector(time, equator)
            val posDiff = ((KM_PER_AU * calcState.position()) - correctState.position()).length()
            val velDiff = ((speedConvert * calcState.velocity()) - correctState.velocity()).length()
            assertTrue(posDiff < posTolerance, "$filename line $lnum: excessive position error $posDiff km")
            assertTrue(velDiff < velTolerance, "$filename line $lnum: excessive velocity error $velDiff km/hour")
        }
    }

    //----------------------------------------------------------------------------------------

    private fun VerifyGeoid(observer: Observer, time: Time, equator: EquatorEpoch) {
        val degreeTolerance = 1.0e-12
        val meterTolerance = 1.0e-8
        val vector = observer.toVector(time, equator)
        val check = vector.toObserver(equator)
        val latDiff = abs(check.latitude - observer.latitude)
        val lonDiff = dcos(observer.latitude) * abs(check.longitude - observer.longitude)
        val heightDiff = abs(check.height - observer.height)
        assertTrue(latDiff < degreeTolerance, "excessive latitude error = $latDiff")
        assertTrue(lonDiff < degreeTolerance, "excessive longitude error = $lonDiff")
        assertTrue(heightDiff < meterTolerance, "excessive height error = $heightDiff")
    }

    @Test
    fun `Verify inverse geoid calculations`() {
        // Calculate position vectors for a variety of geographic coordinates.
        // Verify that we can convert each observer back into a matching vector.
        val time = Time(8134.392799058808)     // 2022-04-09T21:25:37.839Z
        var latitude = -85.0
        while (latitude <= +85.0) {
            var longitude = -175.0
            while (longitude <= +175.0) {
                var height = -5000.0
                while (height <= +5000.0) {
                    val observer = Observer(latitude, longitude, height)
                    VerifyGeoid(observer, time, EquatorEpoch.OfDate)
                    VerifyGeoid(observer, time, EquatorEpoch.J2000)
                    height += 1000.0
                }
                longitude += 5.0
            }
            latitude += 5.0
        }
    }

    //----------------------------------------------------------------------------------------

    @Test
    fun `Angle between body and Sun`() {
        val time = Time(8135.494708473634)     // 2022-04-10T23:52:22.812Z
        VerifyAngleFromSun(time, Body.Sun,       0.0)
        VerifyAngleFromSun(time, Body.Moon,    108.91126797125484)
        VerifyAngleFromSun(time, Body.Venus,    45.2998026289299)
        VerifyAngleFromSun(time, Body.Mars,     54.22257479675931)
        VerifyAngleFromSun(time, Body.Jupiter,  27.49401812548837)
    }

    private fun VerifyAngleFromSun(time: Time, body: Body, correctAngle: Double) {
        val angle = angleFromSun(body, time)
        val diff = abs(correctAngle - angle)
        assertTrue(diff < 1.0e-13, "Excessive angle error $diff for $body at $time")
    }

    //----------------------------------------------------------------------------------------

    @Test
    fun `Lunar eclipse search`() {
        val filename = dataRootDir + "eclipse/lunar_eclipse.txt"
        val infile = File(filename)
        var lnum = 0
        val diffLimit = 2.0     // maximum tolerable error in minutes
        var eclipse: LunarEclipseInfo = searchLunarEclipse(Time(1701, 1, 1, 0, 0, 0.0))
        for (line in infile.readLines()) {
            ++lnum
            val tokens = tokenize(line, 3, filename, lnum)
            val peakTime = parseDate(tokens[0])
            val partialMinutes = tokens[1].toDouble()
            val totalMinutes = tokens[2].toDouble()

            // Verify that the calculated semi-durations are consistent with the kind of eclipse.
            val valid: Boolean = when (eclipse.kind) {
                EclipseKind.Penumbral -> (eclipse.sdPenum > 0.0) && (eclipse.sdPartial == 0.0) && (eclipse.sdTotal == 0.0)
                EclipseKind.Partial   -> (eclipse.sdPenum > 0.0) && (eclipse.sdPartial >  0.0) && (eclipse.sdTotal == 0.0)
                EclipseKind.Total     -> (eclipse.sdPenum > 0.0) && (eclipse.sdPartial >  0.0) && (eclipse.sdTotal >  0.0)
                else                  -> fail("Invalid lunar eclipse kind: ${eclipse.kind}")
            }
            assertTrue(valid, "$filename line $lnum: invalid semiduration(s) for kind ${eclipse.kind}")

            // Check eclipse peak time
            val peakDiffDays = eclipse.peak.ut - peakTime.ut

            // Tolerate missing penumbral eclipses.
            // Skip to the next input line without calculating the next eclipse.
            if (partialMinutes == 0.0 && peakDiffDays > 20.0)
                continue

            val peakDiffMinutes = MINUTES_PER_DAY * abs(peakDiffDays)
            assertTrue(peakDiffMinutes < diffLimit, "excessive peak time error = $peakDiffMinutes minutes; expected $peakTime, found ${eclipse.peak}")

            val partialDiffMinutes = abs(partialMinutes - eclipse.sdPartial)
            assertTrue(partialDiffMinutes <= diffLimit, "excessive partial eclipse semiduration error = $partialDiffMinutes minutes")

            val totalDiffMinutes = abs(totalMinutes - eclipse.sdTotal)
            assertTrue(totalDiffMinutes <= diffLimit, "excessive total eclipse semiduration error = $totalDiffMinutes minutes")

            // Find the next consecutive lunar eclipse.
            eclipse = nextLunarEclipse(eclipse.peak)
        }
    }

    //----------------------------------------------------------------------------------------

    @Test
    fun `Global solar eclipse search`() {
        val filename = dataRootDir + "eclipse/solar_eclipse.txt"
        val infile = File(filename)
        var lnum = 0
        var skipCount = 0
        var eclipse: GlobalSolarEclipseInfo = searchGlobalSolarEclipse(Time(1701, 1, 1, 0, 0, 0.0))
        for (line in infile.readLines()) {
            ++lnum
            val token = tokenize(line, 5, filename, lnum)
            val peak = parseDate(token[0])
            val typeChar = token[2]
            val lat = token[3].toDouble()
            val lon = token[4].toDouble()
            val expectedKind = when(typeChar) {
                "P"  -> EclipseKind.Partial
                "A"  -> EclipseKind.Annular
                "T"  -> EclipseKind.Total
                "H"  -> EclipseKind.Total
                else -> fail("$filename line $lnum: invalid eclipse kind '$typeChar'")
            }

            var diffDays = eclipse.peak.ut - peak.ut

            // Sometimes we find marginal eclipses that aren't listed in the test data.
            // Ignore them if the distance between the Sun/Moon shadow axis and the Earth's center is large.
            while (diffDays < -25.0 && eclipse.distance > 9000.0) {
                ++skipCount
                eclipse = nextGlobalSolarEclipse(eclipse.peak)
                diffDays = eclipse.peak.ut - peak.ut
            }

            // Validate the eclipse prediction.
            val diffMinutes = MINUTES_PER_DAY * abs(diffDays)
            assertTrue(diffMinutes < 6.93, "$filename line $lnum: excessive time error = $diffMinutes minutes; expected $peak, found ${eclipse.peak}")

            // Validate the eclipse kind, but only when it is not a "glancing" eclipse.
            if (eclipse.distance < 6360.0)
                assertTrue(eclipse.kind == expectedKind, "$filename line $lnum: wrong eclipse kind; expected $expectedKind, found ${eclipse.kind}")

            if (eclipse.kind == EclipseKind.Total || eclipse.kind == EclipseKind.Annular) {
                // When the distance between the Moon's shadow ray and the Earth's center is beyond 6100 km,
                // it creates a glancing blow whose geographic coordinates are excessively sensitive
                // to slight changes in the ray. Therefore, it is unreasonable to count large errors there.
                if (eclipse.distance < 6100.0) {
                    val diffAngle = angleDiff(lat, lon, eclipse.latitude, eclipse.longitude)
                    assertTrue(diffAngle < 0.247, "$filename line $lnum: excessive geographic location error = $diffAngle degrees.")
                }
            }

            eclipse = nextGlobalSolarEclipse(eclipse.peak)
        }

        val expectedCount = 1180
        assertTrue(lnum == expectedCount, "$filename: expected $expectedCount lines, found $lnum")
        assertTrue(skipCount < 3, "$filename: excessive skip count $skipCount")
    }

    private fun angleDiff(alat: Double, alon: Double, blat: Double, blon: Double): Double {
        // Find the angle between two geocentric locations.
        val avec = Spherical(alat, alon, 1.0).toVector(epoch2000)
        val bvec = Spherical(blat, blon, 1.0).toVector(epoch2000)
        return avec.angleWith(bvec)
    }

    //----------------------------------------------------------------------------------------

    @Test
    fun `Local solar eclipse near peak location`() {
        // Re-use the test data for global solar eclipses, only feed the given coordinates
        // into the local solar eclipse predictor as the observer's location.
        // In each case, start the search 20 days before the expected eclipse.
        // Then verify that the peak time and eclipse type is correct in each case.
        val filename = dataRootDir + "eclipse/solar_eclipse.txt"
        val infile = File(filename)
        var lnum = 0
        var skipCount = 0
        for (line in infile.readLines()) {
            ++lnum
            val token = tokenize(line, 5, filename, lnum)
            assertEquals(5, token.size, "$filename line $lnum: wrong token count")
            val peak = parseDate(token[0])
            //val typeChar = token[2]
            val lat = token[3].toDouble()
            val lon = token[4].toDouble()
            val observer = Observer(lat, lon, 0.0)

            // Start the search 20 days before we know the eclipse should peak.
            val searchStart = peak.addDays(-20.0)
            val eclipse = searchLocalSolarEclipse(searchStart, observer)

            // Validate the predicted peak eclipse time.
            val diffDays = eclipse.peak.time.ut - peak.ut
            if (diffDays > 20.0) {
                // Tolerate missing a small number of marginal eclipses.
                ++skipCount
                continue
            }

            val diffMinutes = MINUTES_PER_DAY * abs(diffDays)
            assertTrue(diffMinutes < 7.14, "$filename line $lnum: excessive time error = $diffMinutes minutes")
        }

        assertTrue(skipCount <= 6, "$filename: excessive skip count = $skipCount")
    }

    @Test
    fun `Local solar eclipse away from peak location`() {
        val filename = dataRootDir + "eclipse/local_solar_eclipse.txt"
        val infile = File(filename)
        var lnum = 0
        var verifyCount = 0
        for (rawLine in infile.readLines()) {
            ++lnum
            val line = stripLine(rawLine)
            if (line.length == 0)
                continue
            val token = tokenize(line, 13, filename, lnum)
            val latitude = token[0].toDouble()
            val longitude = token[1].toDouble()
            val typeCode = token[2]
            val p1 = parseDate(token[3])
            val p1alt = token[4].toDouble()
            val t1 = optionalParseDate(token[5])
            val t1alt = token[6].toDouble()
            val peak = parseDate(token[7])
            val peakalt = token[8].toDouble()
            val t2 = optionalParseDate(token[9])
            val t2alt = token[10].toDouble()
            val p2 = parseDate(token[11])
            val p2alt = token[12].toDouble()

            val expectedKind = when (typeCode) {
                "P"  -> EclipseKind.Partial
                "A"  -> EclipseKind.Annular
                "T"  -> EclipseKind.Total
                else -> fail("$filename line $lnum: invalid eclipse type '$typeCode'")
            }

            val observer = Observer(latitude, longitude, 0.0)
            val searchTime = p1.addDays(-20.0)
            val eclipse = searchLocalSolarEclipse(searchTime, observer)
            assertEquals(expectedKind, eclipse.kind, "$filename line $lnum: wrong eclipse kind")
            checkEvent(filename, lnum, "peak", peak, peakalt, eclipse.peak)
            checkEvent(filename, lnum, "partialBegin", p1, p1alt, eclipse.partialBegin)
            checkEvent(filename, lnum, "partialEnd", p2, p2alt, eclipse.partialEnd)
            if (typeCode != "P") {
                assertTrue(t1 != null, "$filename line $lnum: t1 was not supposed to be null.")
                assertTrue(t2 != null, "$filename line $lnum: t2 was not supposed to be null.")
                checkEvent(filename, lnum, "totalBegin", t1, t1alt, eclipse.totalBegin)
                checkEvent(filename, lnum, "totalEnd", t2, t2alt, eclipse.totalEnd)
            } else {
                assertTrue(eclipse.totalBegin == null, "$filename line $lnum: eclipse.totalBegin was supposed to be null.")
                assertTrue(eclipse.totalEnd == null, "$filename line $lnum: eclipse.totalEnd was supposed to be null.")
            }
            ++verifyCount;
        }
        assertEquals(6, verifyCount)
    }

    private fun checkEvent(
        filename: String,
        lnum: Int,
        name: String,
        expectedTime: Time,
        expectedAltitude: Double,
        evt: EclipseEvent?
    ) {
        assertTrue(evt != null, "$filename line $lnum: eclipse $name was not supposed to be null.")
        val diffMinutes = MINUTES_PER_DAY * abs(expectedTime.ut - evt.time.ut)
        assertTrue(diffMinutes < 1.0, "$filename line $lnum: excessive time error for $name: $diffMinutes minutes.")
        val diffAlt = abs(expectedAltitude - evt.altitude)
        assertTrue(diffAlt < 0.5, "$filename line $lnum: excessive altitude error for $name: $diffAlt degrees.")
    }

    //----------------------------------------------------------------------------------------

    @Test
    fun `Heliocentric ecliptic longitude`() {
        // These test data were generated using the Python version of Astronomy Engine.
        val time = Time(8141.199164526366)   // 2022-04-16T16:46:47.815Z
        verifyEclipticLongitude(time, Body.Mercury,  94.721630860364  )
        verifyEclipticLongitude(time, Body.Venus,   265.7830489493118 )
        verifyEclipticLongitude(time, Body.Earth,   206.37462137776242)
        verifyEclipticLongitude(time, Body.EMB,     206.37458597072734)
        verifyEclipticLongitude(time, Body.Moon,    206.37171443948878)
        verifyEclipticLongitude(time, Body.Mars,    294.9980436292426 )
        verifyEclipticLongitude(time, Body.Jupiter, 348.4960299261051 )
        verifyEclipticLongitude(time, Body.Saturn,  317.84336759149676)
        verifyEclipticLongitude(time, Body.Uranus,   44.28789479945196)
        verifyEclipticLongitude(time, Body.Neptune, 352.78484144085974)
        verifyEclipticLongitude(time, Body.Pluto,   296.58427291483525)
        verifyEclipticLongitude(time, Body.SSB,     344.6690770393668 )
    }

    private fun verifyEclipticLongitude(time: Time, body: Body, expectedLongitude: Double) {
        val calculatedLongitude = eclipticLongitude(body, time)
        val diff = abs(calculatedLongitude - expectedLongitude)
        assertTrue(diff < 1.0e-16, "excessive ecliptic longitude error for $body: $diff degrees")
    }

    //----------------------------------------------------------------------------------------
    // Generic utility function for loading a series of state vectors
    // from a JPL Horizons output file.

    private class JplStateRecord(
        val lnum: Int,              // the line number where the state vector ends in the JPL Horizons text file
        val state: StateVector      // the state vector itself: position, velocity, and time
    )

    private fun jplHorizonsStateVectors(filename: String): List<JplStateRecord> {
        val infile = File(filename)
        var lnum = 0
        var foundBegin = false
        var part = 0
        var time: Time? = null
        val pos = arrayOf(0.0, 0.0, 0.0)
        val vel = arrayOf(0.0, 0.0, 0.0)
        var list = ArrayList<JplStateRecord>()
        var match: MatchResult?
        val regexPos = Regex("""\s*X =\s*(\S+) Y =\s*(\S+) Z =\s*(\S+)""")
        val regexVel = Regex("""\s*VX=\s*(\S+) VY=\s*(\S+) VZ=\s*(\S+)""")
        for (line in infile.readLines()) {
            ++lnum
            if (!foundBegin) {
                if (line == "\$\$SOE")
                    foundBegin = true
            } else {
                // Input comes in triplets of lines:
                //
                // 2444249.500000000 = A.D. 1980-Jan-11 00:00:00.0000 TDB
                // X =-3.314860345089456E-01 Y = 8.463418210972562E-01 Z = 3.667227830514760E-01
                // VX=-1.642704711077836E-02 VY=-5.494770742558920E-03 VZ=-2.383170237527642E-03
                //
                // Track which of these 3 cases we are in using the 'part' variable...
                when (part) {
                    0 -> {
                        if (line == "\$\$EOE")        // end-of-data marker: the list is complete.
                            return list

                        // 2444249.500000000 = A.D. 1980-Jan-11 00:00:00.0000 TDB
                        // Convert JD to J2000 TT.
                        val tt = tokenize(line)[0].toDouble() - 2451545.0
                        time = Time.fromTerrestrialTime(tt)
                    }
                    1 -> {
                        match = regexPos.matchEntire(line)
                        assertTrue(match != null, "$filename line $lnum: cannot parse position vector")
                        pos[0] = groupDouble(match, 1, filename, lnum)
                        pos[1] = groupDouble(match, 2, filename, lnum)
                        pos[2] = groupDouble(match, 3, filename, lnum)
                    }
                    2 -> {
                        match = regexVel.matchEntire(line)
                        assertTrue(match != null, "$filename line $lnum: cannot parse velocity vector")
                        vel[0] = groupDouble(match, 1, filename, lnum)
                        vel[1] = groupDouble(match, 2, filename, lnum)
                        vel[2] = groupDouble(match, 3, filename, lnum)
                        assertTrue(time != null, "$filename line $lnum: time value was never parsed")
                        val state = StateVector(pos[0], pos[1], pos[2], vel[0], vel[1], vel[2], time)
                        list.add(JplStateRecord(lnum, state))
                        time = null
                    }
                    else -> fail("$filename line $lnum: invalid part=$part")
                }
                part = (part + 1) % 3
            }
        }
        fail("$filename: never found end-of-data marker")
    }

    private fun stateVectorDiff(relative: Boolean, expected: Vector, calculated: Vector): Double {
        val diff = (expected - calculated).length()
        return if (relative)
            diff / expected.length()
        else
            diff
    }

    private fun verifyStateBody(
        relativeFileName: String,
        rThresh: Double,
        vThresh: Double,
        func: (Time) -> StateVector
    ) {
        var filename = dataRootDir + relativeFileName
        for (rec in jplHorizonsStateVectors(filename)) {
            val state = func(rec.state.t)
            val rdiff = stateVectorDiff(rThresh > 0.0, rec.state.position(), state.position())
            val vdiff = stateVectorDiff(vThresh > 0.0, rec.state.velocity(), state.velocity())
            assertTrue(rdiff < abs(rThresh), "$filename line ${rec.lnum}: excessive position error = $rdiff")
            assertTrue(vdiff < abs(vThresh), "$filename line ${rec.lnum}: excessive velocity error = $vdiff")
        }
    }

    //----------------------------------------------------------------------------------------

    @Test
    fun `Barycentric state vectors`() {
        verifyStateBody("barystate/Sun.txt",     -1.224e-05, -1.134e-07) { time -> baryState(Body.Sun,     time) }
        verifyStateBody("barystate/Mercury.txt",  1.672e-04,  2.698e-04) { time -> baryState(Body.Mercury, time) }
        verifyStateBody("barystate/Venus.txt",    4.123e-05,  4.308e-05) { time -> baryState(Body.Venus,   time) }
        verifyStateBody("barystate/Earth.txt",    2.296e-05,  6.359e-05) { time -> baryState(Body.Earth,   time) }
        verifyStateBody("barystate/Mars.txt",     3.107e-05,  5.550e-05) { time -> baryState(Body.Mars,    time) }
        verifyStateBody("barystate/Jupiter.txt",  7.389e-05,  2.471e-04) { time -> baryState(Body.Jupiter, time) }
        verifyStateBody("barystate/Saturn.txt",   1.067e-04,  3.220e-04) { time -> baryState(Body.Saturn,  time) }
        verifyStateBody("barystate/Uranus.txt",   9.035e-05,  2.519e-04) { time -> baryState(Body.Uranus,  time) }
        verifyStateBody("barystate/Neptune.txt",  9.838e-05,  4.446e-04) { time -> baryState(Body.Neptune, time) }
        verifyStateBody("barystate/Pluto.txt",    4.259e-05,  7.827e-05) { time -> baryState(Body.Pluto,   time) }
        verifyStateBody("barystate/Moon.txt",     2.354e-05,  6.604e-05) { time -> baryState(Body.Moon,    time) }
        verifyStateBody("barystate/EMB.txt",      2.353e-05,  6.511e-05) { time -> baryState(Body.EMB,     time) }
    }

    @Test
    fun `Geocentric state vectors`() {
        verifyStateBody("barystate/GeoMoon.txt",  4.086e-05,  5.347e-05) { time -> geoMoonState(time) }
        verifyStateBody("barystate/GeoEMB.txt",   4.076e-05,  5.335e-05) { time -> geoEmbState(time)  }
    }

    @Test
    fun `Heliocentric state vectors`() {
        verifyStateBody("heliostate/SSB.txt",     -1.209e-05, -1.125e-07) { time -> helioState(Body.SSB,     time) }
        verifyStateBody("heliostate/Mercury.txt",  1.481e-04,  2.756e-04) { time -> helioState(Body.Mercury, time) }
        verifyStateBody("heliostate/Venus.txt",    3.528e-05,  4.485e-05) { time -> helioState(Body.Venus,   time) }
        verifyStateBody("heliostate/Earth.txt",    1.476e-05,  6.105e-05) { time -> helioState(Body.Earth,   time) }
        verifyStateBody("heliostate/Mars.txt",     3.154e-05,  5.603e-05) { time -> helioState(Body.Mars,    time) }
        verifyStateBody("heliostate/Jupiter.txt",  7.455e-05,  2.562e-04) { time -> helioState(Body.Jupiter, time) }
        verifyStateBody("heliostate/Saturn.txt",   1.066e-04,  3.150e-04) { time -> helioState(Body.Saturn,  time) }
        verifyStateBody("heliostate/Uranus.txt",   9.034e-05,  2.712e-04) { time -> helioState(Body.Uranus,  time) }
        verifyStateBody("heliostate/Neptune.txt",  9.834e-05,  4.534e-04) { time -> helioState(Body.Neptune, time) }
        verifyStateBody("heliostate/Pluto.txt",    4.271e-05,  1.198e-04) { time -> helioState(Body.Pluto,   time) }
        verifyStateBody("heliostate/Moon.txt",     1.477e-05,  6.195e-05) { time -> helioState(Body.Moon,    time) }
        verifyStateBody("heliostate/EMB.txt",      1.476e-05,  6.106e-05) { time -> helioState(Body.EMB,     time) }
    }

    @Test
    fun `Topocentric state vectors`() {
        val observer = Observer(30.0, -80.0, 1000.0)

        verifyStateBody("topostate/Earth_N30_W80_1000m.txt", 2.108e-04, 2.430e-04) {
            time -> -observer.toStateVector(time, EquatorEpoch.J2000)
        }

        verifyStateBody("topostate/EMB_N30_W80_1000m.txt", 7.195e-04, 2.497e-04) {
            time -> geoEmbState(time) - observer.toStateVector(time, EquatorEpoch.J2000)
        }
    }

    //----------------------------------------------------------------------------------------

    @Test
    fun `Search for relative longitude events`() {
        testElongationFile(2018, dataRootDir + "longitude/opposition_2018.txt", 0.0)
        testPlanetLongitudes(Body.Mercury, "inf")
        testPlanetLongitudes(Body.Venus,   "inf")
        testPlanetLongitudes(Body.Mars,    "opp")
        testPlanetLongitudes(Body.Jupiter, "opp")
        testPlanetLongitudes(Body.Saturn,  "opp")
        testPlanetLongitudes(Body.Uranus,  "opp")
        testPlanetLongitudes(Body.Neptune, "opp")
        testPlanetLongitudes(Body.Pluto,   "opp")
    }

    private fun testElongationFile(year: Int, filename: String, targetRelativeLongitude: Double) {
        val infile = File(filename)
        var lnum = 0
        val searchDate = Time(year, 1, 1, 0, 0, 0.0)
        for (line in infile.readLines()) {
            ++lnum
            // 2018-05-09T00:28Z Jupiter
            val token = tokenize(line, 2, filename, lnum)
            val expectedTime = parseDate(token[0])
            val body = Body.valueOf(token[1])
            val searchResult = searchRelativeLongitude(body, targetRelativeLongitude, searchDate)
            val diffMinutes = MINUTES_PER_DAY * abs(searchResult.tt - expectedTime.tt)
            assertTrue(diffMinutes < 6.8, "$filename line $lnum: excessive search time error = $diffMinutes minutes.")
        }
    }

    private fun testPlanetLongitudes(body: Body, zeroLonEventName: String) {
        val startYear = 1700
        val stopYear = 2200
        var rlon = 0.0
        var count = 0
        var minDiff = 1.0e+99
        var maxDiff = 1.0e+99
        var time = Time(startYear, 1, 1, 0, 0, 0.0)
        val stopTime = Time(stopYear, 1, 1, 0, 0, 0.0)
        val filename = dataRootDir + "temp/k_longitude_$body.txt"
        File(filename).printWriter().use { outfile ->
            while (time.tt < stopTime.tt) {
                ++count
                val eventName = if (rlon == 0.0) zeroLonEventName else "sup"
                val searchResult = searchRelativeLongitude(body, rlon, time)
                assertTrue(searchResult.tt > time.tt, "searchRelativeLongitude went backwards from $time to $searchResult")
                if (count >= 2) {
                    // Check for consistent intervals.
                    // Mainly I don't want to skip over an event!
                    val dayDiff = searchResult.tt - time.tt
                    if (count == 2) {
                        minDiff = dayDiff
                        maxDiff = dayDiff
                    } else if (count > 2) {
                        if (dayDiff < minDiff)
                            minDiff = dayDiff
                        if (dayDiff > maxDiff)
                            maxDiff = dayDiff
                    }
                }

                // Write a line of test data to the output file.
                // This output will be checked as a separate test step
                // later by the `unit_test_kotlin` script (Linux/Mac) or `run.bat` (Windows).
                val geo = geoVector(body, searchResult, Aberration.Corrected)
                val dist = geo.length()
                outfile.println("e $body $eventName ${searchResult.tt} $dist")

                // Search for the opposite longitude event next time.
                time = searchResult
                rlon = 180.0 - rlon
            }
        }

        // Adjust the threshold for how much variation in intervals
        // we tolerate, based on how eccentric the planet's orbit is.
        // The two most eccentric orbits are for Mercury and Mars.
        // All the other planets have reasonably circular orbits.
        val thresh = when (body) {
            Body.Mercury -> 1.65
            Body.Mars    -> 1.30
            else         -> 1.07
        }

        val ratio = maxDiff / minDiff
        assertTrue(ratio < thresh, "excessive event interval ratio $ratio")
    }

    //----------------------------------------------------------------------------------------

    @Test
    fun `Transits of Mercury and Venus`() {
        transitFile(Body.Mercury, dataRootDir + "eclipse/mercury.txt", 10.710, 0.2121)
        transitFile(Body.Venus,   dataRootDir + "eclipse/venus.txt",    9.109, 0.6772)
    }

    private fun transitFile(body: Body, filename: String, limitMinutes: Double, limitSep: Double) {
        val infile = File(filename)
        var lnum = 0
        var transit = searchTransit(body, Time(1600, 1, 1, 0, 0, 0.0))
        for (line in infile.readLines()) {
            ++lnum
            // 22:17 1881-11-08T00:57Z 03:38  3.8633
            val token = tokenize(line, 4, filename, lnum)
            val textp = token[1]
            val text1 = textp.substring(0, 11) + token[0] + "Z"
            val text2 = textp.substring(0, 11) + token[2] + "Z"
            val timep = parseDate(textp)
            var time1 = parseDate(text1)
            var time2 = parseDate(text2)
            val separation = token[3].toDouble()

            // If the start time is after the peak time, it really starts on the previous day.
            if (time1.ut > timep.ut)
                time1 = time1.addDays(-1.0)

            // If the finish time is before the peak time, it really starts on the following day.
            if (time2.ut < timep.ut)
                time2 = time2.addDays(+1.0)

            val diffStart  = MINUTES_PER_DAY * abs(time1.ut - transit.start.ut )
            val diffPeak   = MINUTES_PER_DAY * abs(timep.ut - transit.peak.ut  )
            val diffFinish = MINUTES_PER_DAY * abs(time2.ut - transit.finish.ut)
            val diffSep = abs(separation - transit.separation)

            assertTrue(diffStart  < limitMinutes, "$filename line $lnum: excessive error in transit start = $diffStart minutes")
            assertTrue(diffPeak   < limitMinutes, "$filename line $lnum: excessive error in transit peak = $diffPeak minutes")
            assertTrue(diffFinish < limitMinutes, "$filename line $lnum: excessive error in transit finish = $diffFinish minutes")
            assertTrue(diffSep < limitSep, "$filename line $lnum: excessive error in transit separation = $diffSep arcmin")

            transit = nextTransit(body, transit.finish)
        }
    }

    //----------------------------------------------------------------------------------------
}
