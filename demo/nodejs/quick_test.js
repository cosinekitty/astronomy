// Quick test - you can run this with: node -e "const {calculateSolarTerm} = require('./lichun_calculator'); console.log(calculateSolarTerm(315, 2024));"

const { calculateSolarTerm } = require('./lichun_calculator');

// Get command line arguments
const longitude = process.argv[2] || 315;  // Default to Lichun (315°)
const year = process.argv[3] || 2024;      // Default to 2024

console.log(`Calculating solar term for longitude ${longitude}° in year ${year}:`);
const result = calculateSolarTerm(parseInt(longitude), parseInt(year));

if (result) {
    console.log(`Result: ${result.date.toISOString()}`);
    console.log(`Date: ${result.date.toDateString()}`);
} else {
    console.log('No solar term found for the given parameters.');
} 