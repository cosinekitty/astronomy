@file:OptIn(org.jetbrains.kotlin.gradle.ExperimentalWasmDsl::class)

import org.jetbrains.kotlin.gradle.dsl.JvmTarget

plugins {
    alias(libs.plugins.multiplatform)
    alias(libs.plugins.android.kmp.library)
    `maven-publish`
}

group = "io.github.cosinekitty"
version = "2.1.19"

repositories {
    mavenCentral()
    google()
}

kotlin {
    android {
        namespace = "io.github.cosinekitty.astronomy"
        compileSdk = 36
        minSdk = 23
        compilerOptions {
            jvmTarget = JvmTarget.JVM_17
        }
    }
    jvm {
        compilerOptions {
            jvmTarget = JvmTarget.JVM_17
        }
    }
    iosArm64()
    iosSimulatorArm64()
    macosArm64()
    linuxX64()
    mingwX64()
    wasmJs { nodejs() }
    wasmWasi { nodejs() }

    sourceSets {
        commonTest.dependencies {
            implementation(kotlin("test"))
            implementation(libs.kotlinx.io.core)
        }
    }
}

tasks.withType<Test>().configureEach {
    workingDir = projectDir
}

// NOTE: do not override artifactId here. The Kotlin Multiplatform plugin creates
// one publication per target and derives its coordinates automatically
// (astronomy, astronomy-jvm, astronomy-android, astronomy-linuxx64, ...).
// Forcing `artifactId = project.name` collapses every target onto the single
// coordinate io.github.cosinekitty:astronomy, so the targets overwrite each
// other and consumers resolve whichever variant happened to publish last.
// groupId and version already come from the project-level `group`/`version`.
