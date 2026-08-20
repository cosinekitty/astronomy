@file:OptIn(ExperimentalWasmDsl::class, ExperimentalKotlinGradlePluginApi::class)

import org.jetbrains.kotlin.gradle.ExperimentalKotlinGradlePluginApi
import org.jetbrains.kotlin.gradle.ExperimentalWasmDsl


plugins {
    alias(libs.plugins.multiplatform)
}

group = "io.github.cosinekitty.astronomy.demo"
version = "1.0.0"

repositories {
    mavenCentral()
}

kotlin {
    jvm {
        mainRun {
            mainClass.set("MainKt")
        }
    }

    wasmWasi {
        nodejs()
        binaries.executable()
    }

    linuxX64 {
        binaries {
            executable {
                entryPoint = "main"
            }
        }
    }

    mingwX64 {
        binaries {
            executable {
                entryPoint = "main"
            }
        }
    }

    sourceSets {
        commonMain.dependencies {
            implementation("io.github.cosinekitty:astronomy")
            implementation(libs.kotlinx.datetime)
        }
    }
}
