#!/usr/bin/env node

/**
 * seoul_solar_time.js - True Solar Time for Seoul, South Korea
 * 
 * This program calculates the true solar time for Seoul, South Korea.
 * Seoul coordinates: 37.5665° N, 126.9780° E
 * 
 * True solar time is based on the actual position of the Sun in the sky,
 * where 12:00 represents when the Sun is at its highest point (solar noon).
 */

const Astronomy = require('./astronomy.js');

// Seoul, South Korea coordinates
const SEOUL_LATITUDE = 37.5665;   // North
const SEOUL_LONGITUDE = 126.9780; // East
const SEOUL_ELEVATION = 38;       // meters above sea level

function formatTime(hours) {
    // Convert decimal hours to HH:MM:SS.mmm format
    const totalMilliseconds = Math.round(hours * 3.6e6);
    const totalSeconds = Math.floor(totalMilliseconds / 1000);
    const milliseconds = totalMilliseconds % 1000;
    const totalMinutes = Math.floor(totalSeconds / 60);
    const seconds = totalSeconds % 60;
    const totalHours = Math.floor(totalMinutes / 60);
    const minutes = totalMinutes % 60;
    const displayHours = totalHours % 24;
    
    return `${displayHours.toString().padStart(2, '0')}:${minutes.toString().padStart(2, '0')}:${seconds.toString().padStart(2, '0')}.${milliseconds.toString().padStart(3, '0')}`;
}

function calculateTrueSolarTime(date) {
    // Create observer for Seoul
    const observer = new Astronomy.Observer(SEOUL_LATITUDE, SEOUL_LONGITUDE, SEOUL_ELEVATION);
    
    // Calculate hour angle of the Sun
    const hourAngle = Astronomy.HourAngle(Astronomy.Body.Sun, date, observer);
    
    // Convert hour angle to true solar time
    // True solar time = (hour angle + 12) % 24
    // This gives us the time since solar midnight
    const solarTimeHours = (hourAngle + 12.0) % 24.0;
    
    return {
        hourAngle: hourAngle,
        solarTimeHours: solarTimeHours,
        formattedTime: formatTime(solarTimeHours)
    };
}

function main() {
    console.log("=== True Solar Time for Seoul, South Korea ===\n");
    console.log(`Location: Seoul, South Korea`);
    console.log(`Coordinates: ${SEOUL_LATITUDE}° N, ${SEOUL_LONGITUDE}° E`);
    console.log(`Elevation: ${SEOUL_ELEVATION} meters\n`);
    
    // Get current time
    const now = new Date();
    console.log(`Current UTC time: ${now.toISOString()}`);
    console.log(`Current local time: ${now.toString()}\n`);
    
    // Calculate true solar time for current moment
    const solarTime = calculateTrueSolarTime(now);

    // Calculate current local time in Seoul (in hours)
    // Seoul is UTC+9
    const localTime = new Date(now.getTime() + 9 * 60 * 60 * 1000);
    const localTimeHours = localTime.getHours() + localTime.getMinutes() / 60 + localTime.getSeconds() / 3600 + localTime.getMilliseconds() / 3600000;

    // Calculate offset (true solar time - local clock time)
    let offsetMinutes = (solarTime.solarTimeHours - localTimeHours) * 60;
    // Normalize offset to [-720, +720) minutes for clarity
    if (offsetMinutes < -720) offsetMinutes += 1440;
    if (offsetMinutes >= 720) offsetMinutes -= 1440;

    // Calculate adjusted local time (what the clock would read if it matched the Sun's position)
    let adjustedLocalTimeHours = (localTimeHours + offsetMinutes / 60) % 24;
    if (adjustedLocalTimeHours < 0) adjustedLocalTimeHours += 24;
    const adjustedLocalTimeStr = formatTime(adjustedLocalTimeHours);

    console.log("=== Current True Solar Time ===");
    console.log(`Hour angle of Sun: ${solarTime.hourAngle.toFixed(4)} hours`);
    console.log(`True solar time: ${solarTime.solarTimeHours.toFixed(4)} hours`);
    console.log(`Formatted: ${solarTime.formattedTime}\n`);
    console.log(`Current local time in Seoul: ${formatTime(localTimeHours)} (hours: ${localTimeHours.toFixed(4)})`);
    console.log(`Offset (true solar time - local time): ${offsetMinutes.toFixed(1)} minutes`);
    console.log(`Adjusted local time (if clock matched Sun): ${adjustedLocalTimeStr}\n`);
    
    // Calculate for different times of day
    console.log("=== True Solar Time Throughout the Day ===");
    
    const times = [
        { name: "Sunrise (approximate)", hour: 6 },
        { name: "Morning", hour: 9 },
        { name: "Solar Noon (should be ~12:00)", hour: 12 },
        { name: "Afternoon", hour: 15 },
        { name: "Sunset (approximate)", hour: 18 },
        { name: "Midnight", hour: 0 }
    ];
    
    times.forEach(({ name, hour }) => {
        const testDate = new Date(now);
        testDate.setUTCHours(hour, 0, 0, 0);
        
        const timeSolarTime = calculateTrueSolarTime(testDate);
        
        console.log(`${name.padEnd(25)} | ${testDate.toISOString().substring(11, 19)} UTC | ${timeSolarTime.formattedTime} solar time`);
    });
    
    console.log("\n=== Explanation ===");
    console.log("True solar time is based on the actual position of the Sun in the sky.");
    console.log("• 12:00 solar time = Sun at highest point (solar noon)");
    console.log("• 00:00 solar time = Sun at lowest point (solar midnight)");
    console.log("• The difference from standard time is due to:");
    console.log("  - Earth's elliptical orbit (equation of time)");
    console.log("  - Observer's longitude relative to time zone meridian");
    console.log("  - Daylight saving time adjustments");
}

if (require.main === module) {
    main();
}

module.exports = {
    calculateTrueSolarTime,
    formatTime,
    SEOUL_LATITUDE,
    SEOUL_LONGITUDE,
    SEOUL_ELEVATION
}; 