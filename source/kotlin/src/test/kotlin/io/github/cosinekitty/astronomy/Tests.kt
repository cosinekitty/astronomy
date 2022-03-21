package io.github.cosinekitty.astronomy

import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Test
import org.junit.jupiter.params.ParameterizedTest
import org.junit.jupiter.params.provider.CsvSource

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
}
