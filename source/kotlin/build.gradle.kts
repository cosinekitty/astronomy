plugins {
    java
    kotlin("jvm") version "1.6.10"
    `maven-publish`
    id("org.jetbrains.dokka") version "1.6.10"
}

group = "io.github.cosinekitty"
version = "0.0.1"

repositories {
    mavenCentral()
}

dependencies {
    dokkaHtmlPlugin("org.jetbrains.dokka:kotlin-as-java-plugin:1.6.10")
    val junit5Version = "5.8.2"
    testImplementation("org.junit.jupiter:junit-jupiter-api:$junit5Version")
    testImplementation("org.junit.jupiter:junit-jupiter-params:$junit5Version")
    testRuntimeOnly("org.junit.jupiter:junit-jupiter-engine:$junit5Version")
    testImplementation(kotlin("test"))
}

tasks.test {
    useJUnitPlatform()
}

configure<JavaPluginExtension> {
    sourceCompatibility = JavaVersion.VERSION_1_8
}

publishing {
    publications {
        register("mavenJava", MavenPublication::class) {
            from(components["java"])
        }
    }
}
