const Astronomy = require('./astronomy');

const SOLAR_TERMS = [
    { name: 'Lichun (立春)', nameEn: 'Beginning of Spring', longitude: 315, season: 'Spring' },
    { name: 'Yushui (雨水)', nameEn: 'Rain Water', longitude: 330, season: 'Spring' },
    { name: 'Jingzhe (惊蛰)', nameEn: 'Awakening of Insects', longitude: 345, season: 'Spring' },
    { name: 'Chunfen (春分)', nameEn: 'Spring Equinox', longitude: 0, season: 'Spring' },
    { name: 'Qingming (清明)', nameEn: 'Clear and Bright', longitude: 15, season: 'Spring' },
    { name: 'Guyu (谷雨)', nameEn: 'Grain Rain', longitude: 30, season: 'Spring' },
    { name: 'Lixia (立夏)', nameEn: 'Beginning of Summer', longitude: 45, season: 'Summer' },
    { name: 'Xiaoman (小满)', nameEn: 'Grain Buds', longitude: 60, season: 'Summer' },
    { name: 'Mangzhong (芒种)', nameEn: 'Grain in Ear', longitude: 75, season: 'Summer' },
    { name: 'Xiazhi (夏至)', nameEn: 'Summer Solstice', longitude: 90, season: 'Summer' },
    { name: 'Xiaoshu (小暑)', nameEn: 'Slight Heat', longitude: 105, season: 'Summer' },
    { name: 'Dashu (大暑)', nameEn: 'Great Heat', longitude: 120, season: 'Summer' },
    { name: 'Liqiu (立秋)', nameEn: 'Beginning of Autumn', longitude: 135, season: 'Autumn' },
    { name: 'Chushu (处暑)', nameEn: 'Stopping the Heat', longitude: 150, season: 'Autumn' },
    { name: 'Bailu (白露)', nameEn: 'White Dew', longitude: 165, season: 'Autumn' },
    { name: 'Qiufen (秋分)', nameEn: 'Autumn Equinox', longitude: 180, season: 'Autumn' },
    { name: 'Hanlu (寒露)', nameEn: 'Cold Dew', longitude: 195, season: 'Autumn' },
    { name: 'Shuangjiang (霜降)', nameEn: 'Frost\'s Descent', longitude: 210, season: 'Autumn' },
    { name: 'Lidong (立冬)', nameEn: 'Beginning of Winter', longitude: 225, season: 'Winter' },
    { name: 'Xiaoxue (小雪)', nameEn: 'Slight Snow', longitude: 240, season: 'Winter' },
    { name: 'Daxue (大雪)', nameEn: 'Great Snow', longitude: 255, season: 'Winter' },
    { name: 'Dongzhi (冬至)', nameEn: 'Winter Solstice', longitude: 270, season: 'Winter' },
    { name: 'Xiaohan (小寒)', nameEn: 'Slight Cold', longitude: 285, season: 'Winter' },
    { name: 'Dahan (大寒)', nameEn: 'Great Cold', longitude: 300, season: 'Winter' }
];

function calculateSolarTerm(longitude, year) {
    // Create search start date (January 1st of the given year)
    const searchStart = new Date(year, 0, 1); // Month is 0-indexed
    const searchDays = 366; // Search for 366 days
    
    // First attempt: search within the target year
    try {
        const astroTime = Astronomy.SearchSunLongitude(longitude, searchStart, searchDays);
        // Only return if the result is in the target year
        if (astroTime && astroTime.date.getFullYear() === year) {
            return astroTime;
        }
        return astroTime;
    } catch (error) {
        console.error('Error searching for solar term:', error);
    }
    
    return null; // No solar term found in the target year
}

/**
 * Get the Lichun (Beginning of Spring) date for a given year
 * @param {number} year - The year to find Lichun for
 * @returns {AstroTime|null} - The AstroTime object for Lichun, or null if not found
 */
function getLichunDate(year) {
    return calculateSolarTerm(315, year); // 315 degrees is Lichun
}

/**
 * Calculate the search period between two Lichun dates
 * @param {number} year - The starting year
 * @returns {Object} - Object containing searchStart (AstroTime) and searchDays (number)
 */
function calculateSearchPeriod(year) {
    // Get Lichun of the specified year
    const lichunStart = getLichunDate(year);
    if (!lichunStart) {
        throw new Error(`Could not find Lichun date for year ${year}`);
    }
    
    // Get Lichun of the next year
    const lichunEnd = getLichunDate(year + 1);
    if (!lichunEnd) {
        throw new Error(`Could not find Lichun date for year ${year + 1}`);
    }
    
    // Calculate the difference in days between the two Lichun dates
    const searchDays = (lichunEnd.date.getTime() - lichunStart.date.getTime()) / (1000 * 60 * 60 * 24);
    
    // Subtract 1 day from the start date to ensure we capture the first Lichun
    // SearchSunLongitude requires dateStart to be earlier than the desired longitude event
    const adjustedStartDate = new Date(lichunStart.date.getTime() - 24 * 60 * 60 * 1000);
    
    return {
        searchStart: adjustedStartDate,
        searchDays: Math.ceil(searchDays) + 1 // Add 1 day to compensate for the earlier start
    };
}

/**
 * Get all 24 solar terms between two Lichun dates
 * @param {number} year - The starting year (between Lichun of this year and next year)
 * @returns {Array} - Array of objects containing solar term information and dates
 */
function getAllSolarTerms(year) {
    try {
        const { searchStart, searchDays } = calculateSearchPeriod(year);
        const solarTerms = [];
        
        // Search for each of the 24 solar terms
        for (const term of SOLAR_TERMS) {
            try {
                const astroTime = Astronomy.SearchSunLongitude(term.longitude, searchStart, searchDays);
                if (astroTime) {
                    solarTerms.push({
                        ...term,
                        date: astroTime.date,
                        astroTime: astroTime
                    });
                }
            } catch (error) {
                console.error(`Error searching for ${term.name}:`, error);
                // Continue with other terms even if one fails
            }
        }
        
        // Sort by date to ensure chronological order
        solarTerms.sort((a, b) => a.date.getTime() - b.date.getTime());
        
        return solarTerms;
    } catch (error) {
        console.error('Error calculating solar terms:', error);
        return [];
    }
}

/**
 * Get solar terms for a specific year with more detailed information
 * @param {number} year - The year to get solar terms for
 * @returns {Object} - Object containing year info, search period, and all solar terms
 */
function getSolarTermsForYear(year) {
    try {
        const searchPeriod = calculateSearchPeriod(year);
        const solarTerms = getAllSolarTerms(year);
        
        return {
            year: year,
            searchPeriod: {
                startDate: searchPeriod.searchStart,
                endDate: new Date(searchPeriod.searchStart.getTime() + searchPeriod.searchDays * 24 * 60 * 60 * 1000),
                searchDays: searchPeriod.searchDays
            },
            solarTerms: solarTerms,
            totalTerms: solarTerms.length
        };
    } catch (error) {
        console.error(`Error getting solar terms for year ${year}:`, error);
        return {
            year: year,
            error: error.message,
            solarTerms: [],
            totalTerms: 0
        };
    }
}

module.exports = { 
    SOLAR_TERMS, 
    calculateSolarTerm, 
    getLichunDate, 
    calculateSearchPeriod, 
    getAllSolarTerms, 
    getSolarTermsForYear 
};
