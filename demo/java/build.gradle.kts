plugins {
    java
    application
}

group = "io.github.cosinekitty.astronomy.demo"
version = "0.0.1"

repositories {
    mavenCentral()
    maven("https://jitpack.io")
}

dependencies {
    implementation("com.github.cosinekitty:astronomy:57a736067a2ef005c05b30d94df08d46d78d5326")
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
