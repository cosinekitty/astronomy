plugins {
    java
    application
}

group = "io.github.cosinekitty.astronomy.demo"
version = "1.0.0"

java {
    toolchain {
        languageVersion.set(JavaLanguageVersion.of(21))
    }
}

repositories {
    mavenCentral()
}

dependencies {
    implementation("io.github.cosinekitty:astronomy")
    testImplementation("org.junit.jupiter:junit-jupiter:5.10.2")
    testRuntimeOnly("org.junit.platform:junit-platform-launcher")
}

application {
    mainClass.set("io.github.cosinekitty.astronomy.demo.Main")
}

tasks.jar {
    manifest.attributes["Main-Class"] = "io.github.cosinekitty.astronomy.demo.Main"
    from(configurations.runtimeClasspath.get().map(::zipTree))
    duplicatesStrategy = DuplicatesStrategy.EXCLUDE
}

tasks.getByName<Test>("test") {
    useJUnitPlatform()
}
