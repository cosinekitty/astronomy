plugins {
    java
    application
}

group = "io.github.cosinekitty.astronomy.demo"
version = "0.0.1"

repositories {
    mavenCentral()
    // maven("https://jitpack.io")
}

dependencies {
    implementation(fileTree("../../source/kotlin/build/libs"))
    // Or resolve it from jitpack like,
    //   implementation("com.github.cosinekitty:astronomy:0.0.1")
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
