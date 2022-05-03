plugins {
    java
    application
}

group = "io.github.cosinekitty.astronomy.demo"
version = "1.0.0"

java {
    toolchain {
        languageVersion.set(JavaLanguageVersion.of(11))
    }
}

repositories {
    mavenCentral()
    // maven("https://jitpack.io")
}

dependencies {
    implementation(fileTree("../../source/kotlin/build/libs"))
    implementation("org.jetbrains.kotlin:kotlin-stdlib:1.6.20") // Not needed if resolved from jitpack
    // In an independent project resolve it from jitpack like,
    //   implementation("com.github.cosinekitty:astronomy:x.y.z")
    testImplementation("org.junit.jupiter:junit-jupiter-api:5.8.2")
    testRuntimeOnly("org.junit.jupiter:junit-jupiter-engine")
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
