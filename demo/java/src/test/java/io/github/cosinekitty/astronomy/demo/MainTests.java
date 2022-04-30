package io.github.cosinekitty.astronomy.demo;

import io.github.cosinekitty.astronomy.Time;
import org.junit.jupiter.api.Test;

import java.time.Instant;

import static org.junit.jupiter.api.Assertions.assertEquals;

public class MainTests {
    @Test
    void main() {
        timeTest();
    }

    private void timeTest() {
        String text = "2022-04-29T12:34:45.321Z";
        long millis = Instant.parse(text).toEpochMilli();
        String check = Time.fromMillisecondsSince1970(millis).toString();
        assertEquals(text, check);
    }
}
