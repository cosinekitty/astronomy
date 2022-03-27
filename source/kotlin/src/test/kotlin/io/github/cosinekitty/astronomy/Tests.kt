package io.github.cosinekitty.astronomy

import java.io.File
import org.junit.jupiter.api.Test
import org.junit.jupiter.params.ParameterizedTest
import org.junit.jupiter.params.provider.CsvSource
import kotlin.test.assertEquals
import kotlin.test.assertTrue
import kotlin.math.abs
import kotlin.text.Regex

private const val dataRootDir = "../../generate/"

private fun tokenize(s: String): List<String> {
    return s.split(Regex("""\s+"""))
}

class Tests {

    @Test
    fun `test deltaT calculation`() {
        val time = AstroTime(0.0)
        assertEquals(0.0007389709440951036, time.tt)
    }

    @ParameterizedTest
    @CsvSource(
        value = [
            "2000, 1, 1, 12, 0, 0.0, 0.0, '2000-01-01T12:00:00.000Z'",
            "2022, 1, 1, 12, 0, 0.0, 8036.0, '2022-01-01T12:00:00.000Z'",
            "2022, 1, 1, 18, 0, 0.0, 8036.25, '2022-01-01T18:00:00.000Z'",
            "1970, 12, 13, 23, 45, 12.345, -10610.510273784723, '1970-12-13T23:45:12.345Z'",
            "2022, 1, 1, 18, 59, 59.9999, 8036.291666655093, '2022-01-01T18:59:59.999Z'",
        ]
    )
    fun `universal time calculation should match expectations`(
        year: Int, month: Int, day: Int, hour: Int, minute: Int, second: Double, ut: Double, expectedToString: String
    ) {
        val time = AstroTime(year, month, day, hour, minute, second)
        assertEquals(ut, time.ut)
        assertEquals(time.toString(), expectedToString)
    }

    @Test
    fun `AstroTime should be able to add fractional days`() {
        val time = AstroTime(2000, 1, 1, 12, 0, 0.0)
        assertEquals("2000-01-02T18:00:00.000Z", time.addDays(1.25).toString())
    }

    @Test
    fun `TerseVector methods should work as expected`() {
        val ones = TerseVector(1.0, 1.0, 1.0)
        assertEquals(ones, ones + TerseVector.zero)
        assertEquals(ones, ones + TerseVector.zero)
        assertEquals(TerseVector(6.0, 8.0, 4.0), TerseVector(3.0, 4.0, 2.0) * 2.0)
        assertEquals(TerseVector(-1.5, 2.0, -1.0), TerseVector(-3.0, 4.0, -2.0) / 2.0)
        assertEquals(29.0, TerseVector(-3.0, 4.0, -2.0).quadrature)
        assertEquals(5.744562646538029, TerseVector(-2.0, -2.0, 5.0).magnitude)
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

    private fun compareVectors(testName: String, a: AstroVector, b: AstroVector, tolerance: Double = 1.0e-15) {
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
        var r = RotationMatrix.identity.pivot(2, +90.0)

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
        val time = AstroTime(2000, 1, 1, 0, 0, 0.0)
        val v1 = AstroVector(1.0, 2.0, 3.0, time)

        val v2 = r.rotate(v1)

        compareVectors("Pivot #2", v2, AstroVector(
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
            .identity
            .pivot(0, +37.0)
            .pivot(1, -72.0)
            .pivot(2, +23.0)

        compareMatrices("matrix a", a, RotationMatrix(
            +0.28445164312142457, +0.12074255893448678, +0.9510565162951535,
            -0.8389120034878212,  +0.5115089556077241,  +0.18597106962413537,
            -0.46401930253985274, -0.8507525038429493,  +0.24679194491591758
        ))

        val b = RotationMatrix
            .identity
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
        val time = AstroTime(2022, 3, 15, 21, 50, 0.0)
        val gast = Astronomy.siderealTime(time)
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
                val token = tokenize(line.substring(19))
                assertTrue(token.size == 3, "axisBodyTest($filename line $lnum): expected 3 tokens but found ${token.size}")
                val jd = token[0].toDouble()
                val ra = token[1].toDouble()
                val dec = token[2].toDouble()
                val time = AstroTime(jd - 2451545.0)
                val axis = Astronomy.rotationAxis(body, time)

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
}
