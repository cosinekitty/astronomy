package io.github.cosinekitty.astronomy

import kotlinx.io.buffered
import kotlinx.io.files.Path
import kotlinx.io.files.SystemFileSystem
import kotlinx.io.readLine
import kotlinx.io.writeString

fun readFileLines(path: String): List<String> {
    val lines = mutableListOf<String>()
    SystemFileSystem.source(Path(path)).buffered().use { source ->
        while (true) {
            val line = source.readLine() ?: break
            lines.add(line)
        }
    }
    return lines
}

fun writeFileText(path: String, text: String) {
    SystemFileSystem.sink(Path(path)).buffered().use { sink ->
        sink.writeString(text)
    }
}
