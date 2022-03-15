package io.github.cosinekitty.astronomy

fun add(a: Int, b: Int): Int = a + b

class Accumulator {
    var state = 0
        private set

    fun accumulate(x: Int) {
        state += x
    }
}
