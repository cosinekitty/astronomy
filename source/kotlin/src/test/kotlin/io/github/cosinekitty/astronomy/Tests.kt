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
            "2000, 1, 1, 12, 0, 0, 0, '2000-01-01 12:00:00.0 +0000'",
            "2022, 1, 1, 12, 0, 0, 8036, '2022-01-01 12:00:00.0 +0000'",
            "2022, 1, 1, 18, 0, 0, 8036.25, '2022-01-01 18:00:00.0 +0000'",
        ]
    )
    fun `universal time calculation should match expectations`(
        year: Int, month: Int, day: Int, hour: Int, minute: Int, second: Int, deltaT: Double, expectedToString: String
    ) {
        val time = AstroTime(year, month, day, hour, minute, second)
        assertEquals(deltaT, time.ut)
        assertEquals(time.toString(), expectedToString)
    }
}
