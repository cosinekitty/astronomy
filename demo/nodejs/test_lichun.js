const { SOLAR_TERMS, calculateSolarTerm } = require('./lichun_calculator');

// Test the function with different solar terms
console.log('=== Solar Terms Calculator Test ===\n');

// Test for Lichun (立春) - Beginning of Spring (longitude 315)
console.log('Testing Lichun (立春) for year 2024:');
const lichun2024 = calculateSolarTerm(315, 2024);
if (lichun2024) {
    console.log(`Lichun 2024: ${lichun2024.date.toISOString()}`);
    console.log(`Date: ${lichun2024.date.toDateString()}`);
} else {
    console.log('Lichun 2024 not found');
}

console.log('\n---');

// Test for Chunfen (春分) - Spring Equinox (longitude 0)
console.log('Testing Chunfen (春分) for year 2024:');
const chunfen2024 = calculateSolarTerm(0, 2024);
if (chunfen2024) {
    console.log(`Chunfen 2024: ${chunfen2024.date.toISOString()}`);
    console.log(`Date: ${chunfen2024.date.toDateString()}`);
} else {
    console.log('Chunfen 2024 not found');
}

console.log('\n---');

// Test for Dongzhi (冬至) - Winter Solstice (longitude 270)
console.log('Testing Dongzhi (冬至) for year 2024:');
const dongzhi2024 = calculateSolarTerm(270, 2024);
if (dongzhi2024) {
    console.log(`Dongzhi 2024: ${dongzhi2024.date.toISOString()}`);
    console.log(`Date: ${dongzhi2024.date.toDateString()}`);
} else {
    console.log('Dongzhi 2024 not found');
}

console.log('\n---');

// Display all solar terms
console.log('All Solar Terms:');
SOLAR_TERMS.forEach((term, index) => {
    console.log(`${index + 1}. ${term.name} (${term.nameEn}) - Longitude: ${term.longitude}° - Season: ${term.season}`);
}); 