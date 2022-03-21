package io.github.cosinekitty.astronomy

import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Test
import org.junit.jupiter.params.ParameterizedTest
import org.junit.jupiter.params.provider.CsvSource

class Tests {
    @Test
    fun `test deltaT calculation`() {
        val time = AstroTime(0.0)
        assertEquals(time.tt, 0.0007389709440951036)
    }

    @ParameterizedTest
    @CsvSource(
        value = [
            "2000, 1, 1, 12, 0, 0, 0",
            "2022, 1, 1, 12, 0, 0, 8036",
            "2022, 1, 1, 18, 0, 0, 8036.25",
        ]
    )
    fun `universal time calculation should match expectations`(
        year: Int, month: Int, day: Int, hour: Int, minute: Int, second: Int, deltaT: Double
    ) {
        val time = AstroTime(year, month, day, hour, minute, second)
        assertEquals(time.ut, deltaT)
    }
}
