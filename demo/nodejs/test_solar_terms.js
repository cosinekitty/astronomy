const { 
    SOLAR_TERMS, 
    calculateSolarTerm, 
    getLichunDate, 
    calculateSearchPeriod, 
    getAllSolarTerms, 
    getSolarTermsForYear 
} = require('./lichun_calculator');

// Test the functions
function testSolarTerms(year = 2024) {
    console.log('=== Solar Terms Calculator Test ===\n');
    
    const testYear = year;
    
    // Test 1: Get Lichun date for a specific year
    console.log(`1. Lichun (Beginning of Spring) for ${testYear}:`);
    const lichunDate = getLichunDate(testYear);
    if (lichunDate) {
        console.log(`   Date: ${lichunDate.date.toISOString()}`);
        console.log(`   Local: ${lichunDate.date.toLocaleString()}`);
    } else {
        console.log('   Not found');
    }
    console.log();
    
    // Test 2: Calculate search period between two Lichun dates
    console.log(`2. Search period between Lichun ${testYear} and Lichun ${testYear + 1}:`);
    try {
        const searchPeriod = calculateSearchPeriod(testYear);
        console.log(`   Start: ${searchPeriod.searchStart.date.toISOString()}`);
        console.log(`   Search Days: ${searchPeriod.searchDays}`);
        
        const endDate = new Date(searchPeriod.searchStart.date.getTime() + searchPeriod.searchDays * 24 * 60 * 60 * 1000);
        console.log(`   End: ${endDate.toISOString()}`);
    } catch (error) {
        console.log(`   Error: ${error.message}`);
    }
    console.log();
    
    // Test 3: Get all 24 solar terms
    console.log(`3. All 24 Solar Terms for ${testYear}:`);
    const allTerms = getAllSolarTerms(testYear);
    console.log(`   Found ${allTerms.length} solar terms:`);
    
    allTerms.forEach((term, index) => {
        const dateStr = term.date.toLocaleDateString('en-US', {
            year: 'numeric',
            month: 'short',
            day: 'numeric',
            hour: '2-digit',
            minute: '2-digit'
        });
        console.log(`   ${(index + 1).toString().padStart(2)}. ${term.nameEn} (${term.name}) - ${dateStr} - ${term.season}`);
    });
    console.log();
    
    // Test 4: Get comprehensive solar terms information
    console.log(`4. Comprehensive Solar Terms Information for ${testYear}:`);
    const comprehensive = getSolarTermsForYear(testYear);
    
    if (comprehensive.error) {
        console.log(`   Error: ${comprehensive.error}`);
    } else {
        console.log(`   Year: ${comprehensive.year}`);
        console.log(`   Search Period:`);
        console.log(`     Start: ${comprehensive.searchPeriod.startDate.toISOString()}`);
        console.log(`     End: ${comprehensive.searchPeriod.endDate.toISOString()}`);
        console.log(`     Duration: ${comprehensive.searchPeriod.searchDays} days`);
        console.log(`   Total Solar Terms Found: ${comprehensive.totalTerms}`);
        
        // Show first few and last few terms
        if (comprehensive.solarTerms.length > 0) {
            console.log(`   First 3 terms:`);
            comprehensive.solarTerms.slice(0, 3).forEach(term => {
                console.log(`     ${term.nameEn} - ${term.date.toLocaleDateString()}`);
            });
            
            if (comprehensive.solarTerms.length > 6) {
                console.log(`   ...`);
                console.log(`   Last 3 terms:`);
                comprehensive.solarTerms.slice(-3).forEach(term => {
                    console.log(`     ${term.nameEn} - ${term.date.toLocaleDateString()}`);
                });
            }
        }
    }
    console.log();
    
    // Test 5: Test individual solar term calculation
    console.log('5. Individual Solar Term Test:');
    const springEquinox = calculateSolarTerm(0, testYear); // 0 degrees = Spring Equinox
    if (springEquinox) {
        console.log(`   Spring Equinox ${testYear}: ${springEquinox.date.toLocaleDateString()}`);
    } else {
        console.log(`   Spring Equinox ${testYear}: Not found`);
    }
}



// Run the test
if (require.main === module) {
    // Get year from command line arguments, default to 2024
    const year = parseInt(process.argv[2]) || 2024;
    
    console.log(`Running tests for year: ${year}\n`);
    
    testSolarTerms(year);
}

module.exports = { testSolarTerms }; 