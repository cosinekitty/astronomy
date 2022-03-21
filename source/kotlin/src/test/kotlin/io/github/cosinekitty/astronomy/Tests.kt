package io.github.cosinekitty.astronomy

import org.junit.jupiter.api.Assertions.assertEquals
import org.junit.jupiter.api.Test
import org.junit.jupiter.params.ParameterizedTest
import org.junit.jupiter.params.provider.CsvSource

class Tests 
{
    @Test
    fun `test DeltaT calculation`()
    {
        val time = AstroTime(0.0);
        assertEquals(time.tt, 0.0007389709440951036);
    }
}
