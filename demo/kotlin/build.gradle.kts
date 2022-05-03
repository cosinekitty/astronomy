import org.jetbrains.kotlin.gradle.tasks.KotlinCompile

plugins {
    kotlin("jvm") version "1.6.21"
    application
}

group = "io.github.cosinekitty.astronomy.demo"
version = "1.0.0"

repositories {
    mavenCentral()
}

dependencies {
    implementation(fileTree("../../source/kotlin/build/libs"))
    testImplementation(kotlin("test"))
}

tasks.test {
    useJUnitPlatform()
}

tasks.withType<KotlinCompile> {
    kotlinOptions.jvmTarget = "11"
    kotlinOptions {
        allWarningsAsErrors = true
    }
}

tasks.jar {
    manifest.attributes["Main-Class"] = "MainKt"
    from(configurations.runtimeClasspath.get().map(::zipTree))
    duplicatesStrategy = DuplicatesStrategy.EXCLUDE
}

application {
    mainClass.set("MainKt")
}