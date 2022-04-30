package io.github.cosinekitty.astronomy.demo;

import io.github.cosinekitty.astronomy.Time;
import org.junit.jupiter.api.Test;

import java.time.Instant;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class MainTests {
    @Test
    void main() {
        String time = "2022-01-01T12:00:00.000Z";
        assertEquals(time, Time.fromMillisecondsSince1970(Instant.parse(time).toEpochMilli()).toString());
    }
}
