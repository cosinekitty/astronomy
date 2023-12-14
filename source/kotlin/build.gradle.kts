plugins {
    java
    kotlin("jvm") version "1.6.10"
    `maven-publish`
    id("org.jetbrains.dokka") version "1.6.10"
}

group = "io.github.cosinekitty"
version = "2.1.19"

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
    sourceCompatibility = JavaVersion.VERSION_11
}

val sourceJar by tasks.creating(Jar::class) {
    dependsOn(tasks["classes"])
    archiveClassifier.set("sources")
    from(sourceSets["main"].allSource)
}

task("fatJar", type = Jar::class) {
    from(configurations.runtimeClasspath.get().map(::zipTree))
    duplicatesStrategy = DuplicatesStrategy.EXCLUDE
    with(tasks.jar.get())
}

publishing {
    publications {
        register("mavenJava", MavenPublication::class) {
            from(components["kotlin"])
            artifact(sourceJar)
        }
    }
}

tasks.withType<org.jetbrains.kotlin.gradle.tasks.KotlinCompile>().configureEach {
    kotlinOptions {
        allWarningsAsErrors = true
    }
}

tasks.dokkaGfm.configure {
    dokkaSourceSets {
        named("main") {
            includeNonPublic.set(false)
            reportUndocumented.set(true)
        }
    }
}
