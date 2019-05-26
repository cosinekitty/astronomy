# Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`struct `[`astro_angle_result_t`](#structastro__angle__result__t) | 
`struct `[`astro_apsis_t`](#structastro__apsis__t) | 
`struct `[`astro_cheb_coeff_t`](#structastro__cheb__coeff__t) | 
`struct `[`astro_cheb_record_t`](#structastro__cheb__record__t) | 
`struct `[`astro_ecliptic_t`](#structastro__ecliptic__t) | 
`struct `[`astro_elongation_t`](#structastro__elongation__t) | 
`struct `[`astro_equatorial_t`](#structastro__equatorial__t) | 
`struct `[`astro_func_result_t`](#structastro__func__result__t) | 
`struct `[`astro_horizon_t`](#structastro__horizon__t) | 
`struct `[`astro_hour_angle_t`](#structastro__hour__angle__t) | 
`struct `[`astro_illum_t`](#structastro__illum__t) | 
`struct `[`astro_moon_quarter_t`](#structastro__moon__quarter__t) | 
`struct `[`astro_observer_t`](#structastro__observer__t) | Represents a location of an observer on (or near) the surface of the Earth.
`struct `[`astro_search_result_t`](#structastro__search__result__t) | 
`struct `[`astro_seasons_t`](#structastro__seasons__t) | 
`struct `[`astro_time_t`](#structastro__time__t) | A date and time used for astronomical calculations.
`struct `[`astro_utc_t`](#structastro__utc__t) | 
`struct `[`astro_vector_t`](#structastro__vector__t) | 
`struct `[`context_peak_altitude_t`](#structcontext__peak__altitude__t) | 
`struct `[`deltat_entry_t`](#structdeltat__entry__t) | 
`struct `[`earth_tilt_t`](#structearth__tilt__t) | 
`struct `[`MoonContext`](#structMoonContext) | 
`struct `[`vsop_formula_t`](#structvsop__formula__t) | 
`struct `[`vsop_model_t`](#structvsop__model__t) | 
`struct `[`vsop_series_t`](#structvsop__series__t) | 
`struct `[`vsop_term_t`](#structvsop__term__t) | 

# struct `astro_angle_result_t` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public `[`astro_status_t`](#astronomy_8h_1a3605f063314f858dec26db307c4e8214)` `[`status`](#structastro__angle__result__t_1aab07512f9f53b26cc2c863147f774e36) | 
`public double `[`angle`](#structastro__angle__result__t_1af8536a6c9f69f165c2900bef800d7069) | 

## Members

#### `public `[`astro_status_t`](#astronomy_8h_1a3605f063314f858dec26db307c4e8214)` `[`status`](#structastro__angle__result__t_1aab07512f9f53b26cc2c863147f774e36) 

#### `public double `[`angle`](#structastro__angle__result__t_1af8536a6c9f69f165c2900bef800d7069) 

# struct `astro_apsis_t` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public `[`astro_status_t`](#astronomy_8h_1a3605f063314f858dec26db307c4e8214)` `[`status`](#structastro__apsis__t_1ac65536f8f131c7bf91ae2d5744666247) | 
`public `[`astro_time_t`](#structastro__time__t)` `[`time`](#structastro__apsis__t_1a3c22db68e225d218c54f401952b63516) | 
`public `[`astro_apsis_kind_t`](#astronomy_8h_1a637f68d8e61bfb5f1e49d6a115a05330)` `[`kind`](#structastro__apsis__t_1a2bf33f8d08caf0eaf226539c22dceb09) | 
`public double `[`dist_au`](#structastro__apsis__t_1ada2d465ace3e0cfbe1547cd122117a7b) | 
`public double `[`dist_km`](#structastro__apsis__t_1a1ff80b5f6e5d1cb863a5ccbcd362007f) | 

## Members

#### `public `[`astro_status_t`](#astronomy_8h_1a3605f063314f858dec26db307c4e8214)` `[`status`](#structastro__apsis__t_1ac65536f8f131c7bf91ae2d5744666247) 

#### `public `[`astro_time_t`](#structastro__time__t)` `[`time`](#structastro__apsis__t_1a3c22db68e225d218c54f401952b63516) 

#### `public `[`astro_apsis_kind_t`](#astronomy_8h_1a637f68d8e61bfb5f1e49d6a115a05330)` `[`kind`](#structastro__apsis__t_1a2bf33f8d08caf0eaf226539c22dceb09) 

#### `public double `[`dist_au`](#structastro__apsis__t_1ada2d465ace3e0cfbe1547cd122117a7b) 

#### `public double `[`dist_km`](#structastro__apsis__t_1a1ff80b5f6e5d1cb863a5ccbcd362007f) 

# struct `astro_cheb_coeff_t` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public double `[`data`](#structastro__cheb__coeff__t_1a5d2769c57b5d74a488f348b56a7aaada) | 

## Members

#### `public double `[`data`](#structastro__cheb__coeff__t_1a5d2769c57b5d74a488f348b56a7aaada) 

# struct `astro_cheb_record_t` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public double `[`tt`](#structastro__cheb__record__t_1ab5716af5a2fa0b23cc9b7885872634c8) | 
`public double `[`ndays`](#structastro__cheb__record__t_1a18fad328273a7ba45b99db03accb93e8) | 
`public int `[`ncoeff`](#structastro__cheb__record__t_1a72137ec8ccee0c8e5252f99cf4e1364a) | 
`public const `[`astro_cheb_coeff_t`](#structastro__cheb__coeff__t)` * `[`coeff`](#structastro__cheb__record__t_1adbf5052cec3a7a9e0ece662b4712b714) | 

## Members

#### `public double `[`tt`](#structastro__cheb__record__t_1ab5716af5a2fa0b23cc9b7885872634c8) 

#### `public double `[`ndays`](#structastro__cheb__record__t_1a18fad328273a7ba45b99db03accb93e8) 

#### `public int `[`ncoeff`](#structastro__cheb__record__t_1a72137ec8ccee0c8e5252f99cf4e1364a) 

#### `public const `[`astro_cheb_coeff_t`](#structastro__cheb__coeff__t)` * `[`coeff`](#structastro__cheb__record__t_1adbf5052cec3a7a9e0ece662b4712b714) 

# struct `astro_ecliptic_t` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public `[`astro_status_t`](#astronomy_8h_1a3605f063314f858dec26db307c4e8214)` `[`status`](#structastro__ecliptic__t_1ae3e7437e7a364d06845d8e58e65a2599) | 
`public double `[`ex`](#structastro__ecliptic__t_1a53bbb066bc4e10bc82f57e24848e716e) | 
`public double `[`ey`](#structastro__ecliptic__t_1a5747b2459a0bb6be11369d73b77a727d) | 
`public double `[`ez`](#structastro__ecliptic__t_1a234a448ded3d661168fd2a802b433148) | 
`public double `[`elat`](#structastro__ecliptic__t_1a0e8377dd629cf2e71e0b5a5afc4a7e97) | 
`public double `[`elon`](#structastro__ecliptic__t_1ad89a8adf33083aff212ca79a98ae92c3) | 

## Members

#### `public `[`astro_status_t`](#astronomy_8h_1a3605f063314f858dec26db307c4e8214)` `[`status`](#structastro__ecliptic__t_1ae3e7437e7a364d06845d8e58e65a2599) 

#### `public double `[`ex`](#structastro__ecliptic__t_1a53bbb066bc4e10bc82f57e24848e716e) 

#### `public double `[`ey`](#structastro__ecliptic__t_1a5747b2459a0bb6be11369d73b77a727d) 

#### `public double `[`ez`](#structastro__ecliptic__t_1a234a448ded3d661168fd2a802b433148) 

#### `public double `[`elat`](#structastro__ecliptic__t_1a0e8377dd629cf2e71e0b5a5afc4a7e97) 

#### `public double `[`elon`](#structastro__ecliptic__t_1ad89a8adf33083aff212ca79a98ae92c3) 

# struct `astro_elongation_t` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public `[`astro_status_t`](#astronomy_8h_1a3605f063314f858dec26db307c4e8214)` `[`status`](#structastro__elongation__t_1a769739d6a5284136d2e54cb9247acc8a) | 
`public `[`astro_time_t`](#structastro__time__t)` `[`time`](#structastro__elongation__t_1a81a285588e39f9fb3e292ab1305a869e) | 
`public `[`astro_visibility_t`](#astronomy_8h_1a4a3fd16133cdaed41936a10e8520f9ff)` `[`visibility`](#structastro__elongation__t_1a798a24a1449115fbf489083678763cea) | 
`public double `[`elongation`](#structastro__elongation__t_1ad277178bd020945f1e1881e00c323563) | 
`public double `[`relative_longitude`](#structastro__elongation__t_1af14eff240b62a4373e9f130e33730dc5) | 

## Members

#### `public `[`astro_status_t`](#astronomy_8h_1a3605f063314f858dec26db307c4e8214)` `[`status`](#structastro__elongation__t_1a769739d6a5284136d2e54cb9247acc8a) 

#### `public `[`astro_time_t`](#structastro__time__t)` `[`time`](#structastro__elongation__t_1a81a285588e39f9fb3e292ab1305a869e) 

#### `public `[`astro_visibility_t`](#astronomy_8h_1a4a3fd16133cdaed41936a10e8520f9ff)` `[`visibility`](#structastro__elongation__t_1a798a24a1449115fbf489083678763cea) 

#### `public double `[`elongation`](#structastro__elongation__t_1ad277178bd020945f1e1881e00c323563) 

#### `public double `[`relative_longitude`](#structastro__elongation__t_1af14eff240b62a4373e9f130e33730dc5) 

# struct `astro_equatorial_t` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public `[`astro_status_t`](#astronomy_8h_1a3605f063314f858dec26db307c4e8214)` `[`status`](#structastro__equatorial__t_1afab33097965cc037477a3cceefdedf17) | 
`public double `[`ra`](#structastro__equatorial__t_1a3f0d378983c16e181bae1787fc22fc9c) | 
`public double `[`dec`](#structastro__equatorial__t_1af313e351ebecda5b3d4c00dccff26278) | 
`public double `[`dist`](#structastro__equatorial__t_1ae436ef560df6767e82a5e5f87af13a9a) | 

## Members

#### `public `[`astro_status_t`](#astronomy_8h_1a3605f063314f858dec26db307c4e8214)` `[`status`](#structastro__equatorial__t_1afab33097965cc037477a3cceefdedf17) 

#### `public double `[`ra`](#structastro__equatorial__t_1a3f0d378983c16e181bae1787fc22fc9c) 

#### `public double `[`dec`](#structastro__equatorial__t_1af313e351ebecda5b3d4c00dccff26278) 

#### `public double `[`dist`](#structastro__equatorial__t_1ae436ef560df6767e82a5e5f87af13a9a) 

# struct `astro_func_result_t` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public `[`astro_status_t`](#astronomy_8h_1a3605f063314f858dec26db307c4e8214)` `[`status`](#structastro__func__result__t_1a8e86a498accbd6c245ad737b9c775ef4) | 
`public double `[`value`](#structastro__func__result__t_1abb36158e9c74ccc058a05d9e2440b888) | 

## Members

#### `public `[`astro_status_t`](#astronomy_8h_1a3605f063314f858dec26db307c4e8214)` `[`status`](#structastro__func__result__t_1a8e86a498accbd6c245ad737b9c775ef4) 

#### `public double `[`value`](#structastro__func__result__t_1abb36158e9c74ccc058a05d9e2440b888) 

# struct `astro_horizon_t` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public double `[`azimuth`](#structastro__horizon__t_1afd850c814effd476340247474c5a13fc) | 
`public double `[`altitude`](#structastro__horizon__t_1aa92174981aa9e849fd9d148b94e248a2) | 
`public double `[`ra`](#structastro__horizon__t_1a8be5e63cf8fe54fbe7adb7fd34873529) | 
`public double `[`dec`](#structastro__horizon__t_1a447f35e253764bd0b682d932794c9db9) | 

## Members

#### `public double `[`azimuth`](#structastro__horizon__t_1afd850c814effd476340247474c5a13fc) 

#### `public double `[`altitude`](#structastro__horizon__t_1aa92174981aa9e849fd9d148b94e248a2) 

#### `public double `[`ra`](#structastro__horizon__t_1a8be5e63cf8fe54fbe7adb7fd34873529) 

#### `public double `[`dec`](#structastro__horizon__t_1a447f35e253764bd0b682d932794c9db9) 

# struct `astro_hour_angle_t` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public `[`astro_status_t`](#astronomy_8h_1a3605f063314f858dec26db307c4e8214)` `[`status`](#structastro__hour__angle__t_1ac410b728be11fc8bb31c971b275309e4) | 
`public `[`astro_time_t`](#structastro__time__t)` `[`time`](#structastro__hour__angle__t_1aa7a75ac99b4133967ea4e62330535d34) | 
`public `[`astro_horizon_t`](#structastro__horizon__t)` `[`hor`](#structastro__hour__angle__t_1a7e6b30b3c23dc2dacb3bb80b09ef205e) | 
`public int `[`iter`](#structastro__hour__angle__t_1adbddda57a6c50088943d2053e4be56f2) | 

## Members

#### `public `[`astro_status_t`](#astronomy_8h_1a3605f063314f858dec26db307c4e8214)` `[`status`](#structastro__hour__angle__t_1ac410b728be11fc8bb31c971b275309e4) 

#### `public `[`astro_time_t`](#structastro__time__t)` `[`time`](#structastro__hour__angle__t_1aa7a75ac99b4133967ea4e62330535d34) 

#### `public `[`astro_horizon_t`](#structastro__horizon__t)` `[`hor`](#structastro__hour__angle__t_1a7e6b30b3c23dc2dacb3bb80b09ef205e) 

#### `public int `[`iter`](#structastro__hour__angle__t_1adbddda57a6c50088943d2053e4be56f2) 

# struct `astro_illum_t` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public `[`astro_status_t`](#astronomy_8h_1a3605f063314f858dec26db307c4e8214)` `[`status`](#structastro__illum__t_1a2c5171551ffcd52a6f22de82f83f6e11) | 
`public `[`astro_time_t`](#structastro__time__t)` `[`time`](#structastro__illum__t_1af637ea2ce377bbe504dd3f76d71f160b) | 
`public double `[`mag`](#structastro__illum__t_1a4e869d8238236f3b0a6df0be1711abb6) | 
`public double `[`phase_angle`](#structastro__illum__t_1afc4a65f7fe01c4198c95fd2265219046) | 
`public double `[`helio_dist`](#structastro__illum__t_1a2b92ab3fcc89bc6d648433420c19da5d) | 
`public double `[`ring_tilt`](#structastro__illum__t_1a60be7a165fd30aaf279f335f882e39e9) | 

## Members

#### `public `[`astro_status_t`](#astronomy_8h_1a3605f063314f858dec26db307c4e8214)` `[`status`](#structastro__illum__t_1a2c5171551ffcd52a6f22de82f83f6e11) 

#### `public `[`astro_time_t`](#structastro__time__t)` `[`time`](#structastro__illum__t_1af637ea2ce377bbe504dd3f76d71f160b) 

#### `public double `[`mag`](#structastro__illum__t_1a4e869d8238236f3b0a6df0be1711abb6) 

#### `public double `[`phase_angle`](#structastro__illum__t_1afc4a65f7fe01c4198c95fd2265219046) 

#### `public double `[`helio_dist`](#structastro__illum__t_1a2b92ab3fcc89bc6d648433420c19da5d) 

#### `public double `[`ring_tilt`](#structastro__illum__t_1a60be7a165fd30aaf279f335f882e39e9) 

# struct `astro_moon_quarter_t` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public `[`astro_status_t`](#astronomy_8h_1a3605f063314f858dec26db307c4e8214)` `[`status`](#structastro__moon__quarter__t_1acd7a098d48ed7c5ab904744b988b5414) | 
`public int `[`quarter`](#structastro__moon__quarter__t_1acf3ebb2207eb4d23f3be53302a70dcb5) | 
`public `[`astro_time_t`](#structastro__time__t)` `[`time`](#structastro__moon__quarter__t_1ae726bcbd11b18a504d9e619df91d6c64) | 

## Members

#### `public `[`astro_status_t`](#astronomy_8h_1a3605f063314f858dec26db307c4e8214)` `[`status`](#structastro__moon__quarter__t_1acd7a098d48ed7c5ab904744b988b5414) 

#### `public int `[`quarter`](#structastro__moon__quarter__t_1acf3ebb2207eb4d23f3be53302a70dcb5) 

#### `public `[`astro_time_t`](#structastro__time__t)` `[`time`](#structastro__moon__quarter__t_1ae726bcbd11b18a504d9e619df91d6c64) 

# struct `astro_observer_t` 

Represents a location of an observer on (or near) the surface of the Earth.

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public double `[`latitude`](#structastro__observer__t_1a4ece76c5061138651f91cdc0a43fabea) | Geographic latitude in degrees north (positive) or south (negative) of the equator.
`public double `[`longitude`](#structastro__observer__t_1aa54800e341dd1fcdd8ae257f702c3872) | Geographic longitude in degrees east (positive) or west (negative) of the prime meridian at Greenwich, England.
`public double `[`height`](#structastro__observer__t_1aed0ae0473aa7ad733819288803fbff0b) | The height above (positive) or below (negative) sea level, expressed in meters.

## Members

#### `public double `[`latitude`](#structastro__observer__t_1a4ece76c5061138651f91cdc0a43fabea) 

Geographic latitude in degrees north (positive) or south (negative) of the equator.

#### `public double `[`longitude`](#structastro__observer__t_1aa54800e341dd1fcdd8ae257f702c3872) 

Geographic longitude in degrees east (positive) or west (negative) of the prime meridian at Greenwich, England.

#### `public double `[`height`](#structastro__observer__t_1aed0ae0473aa7ad733819288803fbff0b) 

The height above (positive) or below (negative) sea level, expressed in meters.

# struct `astro_search_result_t` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public `[`astro_status_t`](#astronomy_8h_1a3605f063314f858dec26db307c4e8214)` `[`status`](#structastro__search__result__t_1aafdec98ea6e300d1b5b5ce996d6cda0a) | 
`public `[`astro_time_t`](#structastro__time__t)` `[`time`](#structastro__search__result__t_1afa7bbea2c648a30a5d1c3dd575f62b89) | 
`public int `[`iter`](#structastro__search__result__t_1a8f1a16e1e82b38966c64097df2981df1) | 

## Members

#### `public `[`astro_status_t`](#astronomy_8h_1a3605f063314f858dec26db307c4e8214)` `[`status`](#structastro__search__result__t_1aafdec98ea6e300d1b5b5ce996d6cda0a) 

#### `public `[`astro_time_t`](#structastro__time__t)` `[`time`](#structastro__search__result__t_1afa7bbea2c648a30a5d1c3dd575f62b89) 

#### `public int `[`iter`](#structastro__search__result__t_1a8f1a16e1e82b38966c64097df2981df1) 

# struct `astro_seasons_t` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public `[`astro_status_t`](#astronomy_8h_1a3605f063314f858dec26db307c4e8214)` `[`status`](#structastro__seasons__t_1a1e597ee698867cde851e6c8b0d43695c) | 
`public `[`astro_time_t`](#structastro__time__t)` `[`mar_equinox`](#structastro__seasons__t_1ab986e420149ee6217090ec28d8d44721) | 
`public `[`astro_time_t`](#structastro__time__t)` `[`jun_solstice`](#structastro__seasons__t_1aa70d91b741af06f465c1ce697c85259f) | 
`public `[`astro_time_t`](#structastro__time__t)` `[`sep_equinox`](#structastro__seasons__t_1acb83b7271927a79e4eccce7d86d4a322) | 
`public `[`astro_time_t`](#structastro__time__t)` `[`dec_solstice`](#structastro__seasons__t_1a05f1ac67ae8180267ac6140dc0643d3b) | 

## Members

#### `public `[`astro_status_t`](#astronomy_8h_1a3605f063314f858dec26db307c4e8214)` `[`status`](#structastro__seasons__t_1a1e597ee698867cde851e6c8b0d43695c) 

#### `public `[`astro_time_t`](#structastro__time__t)` `[`mar_equinox`](#structastro__seasons__t_1ab986e420149ee6217090ec28d8d44721) 

#### `public `[`astro_time_t`](#structastro__time__t)` `[`jun_solstice`](#structastro__seasons__t_1aa70d91b741af06f465c1ce697c85259f) 

#### `public `[`astro_time_t`](#structastro__time__t)` `[`sep_equinox`](#structastro__seasons__t_1acb83b7271927a79e4eccce7d86d4a322) 

#### `public `[`astro_time_t`](#structastro__time__t)` `[`dec_solstice`](#structastro__seasons__t_1a05f1ac67ae8180267ac6140dc0643d3b) 

# struct `astro_time_t` 

A date and time used for astronomical calculations.

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public double `[`ut`](#structastro__time__t_1a0874edc886a134a3a09ed3c1e4c41d24) | UT1/UTC number of days since noon on January 1, 2000
`public double `[`tt`](#structastro__time__t_1a7efe42a250b382ca1052fd4ae5da4d8d) | Terrestrial Time days since noon on January 1, 2000

## Members

#### `public double `[`ut`](#structastro__time__t_1a0874edc886a134a3a09ed3c1e4c41d24) 

UT1/UTC number of days since noon on January 1, 2000

#### `public double `[`tt`](#structastro__time__t_1a7efe42a250b382ca1052fd4ae5da4d8d) 

Terrestrial Time days since noon on January 1, 2000

# struct `astro_utc_t` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public int `[`year`](#structastro__utc__t_1aa636614dbe7bd040fdbdcac42f8b3069) | 
`public int `[`month`](#structastro__utc__t_1a16d083a87186eb7a22a527127fba0574) | 
`public int `[`day`](#structastro__utc__t_1a0c92aac85896c7b1f70630035dbcdaa8) | 
`public int `[`hour`](#structastro__utc__t_1a126fe88a150a70a7c1644a2771dd6a31) | 
`public int `[`minute`](#structastro__utc__t_1a1c503dad9736c86aed0630957e1f2646) | 
`public double `[`second`](#structastro__utc__t_1a0cc18479f7a558050021acf9b32caadc) | 

## Members

#### `public int `[`year`](#structastro__utc__t_1aa636614dbe7bd040fdbdcac42f8b3069) 

#### `public int `[`month`](#structastro__utc__t_1a16d083a87186eb7a22a527127fba0574) 

#### `public int `[`day`](#structastro__utc__t_1a0c92aac85896c7b1f70630035dbcdaa8) 

#### `public int `[`hour`](#structastro__utc__t_1a126fe88a150a70a7c1644a2771dd6a31) 

#### `public int `[`minute`](#structastro__utc__t_1a1c503dad9736c86aed0630957e1f2646) 

#### `public double `[`second`](#structastro__utc__t_1a0cc18479f7a558050021acf9b32caadc) 

# struct `astro_vector_t` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public `[`astro_status_t`](#astronomy_8h_1a3605f063314f858dec26db307c4e8214)` `[`status`](#structastro__vector__t_1a56282eb2bcf0dce0170952f9d1aa8b3c) | 
`public double `[`x`](#structastro__vector__t_1a9519f6730bb6e466b810c9b5b053d7a5) | 
`public double `[`y`](#structastro__vector__t_1a9f78a5f7f0320b6f5194e1006488baa5) | 
`public double `[`z`](#structastro__vector__t_1aa96463d4940475e13e1df0fdaabbb955) | 
`public `[`astro_time_t`](#structastro__time__t)` `[`t`](#structastro__vector__t_1ad90a34927631a1940ae0378574d42103) | 

## Members

#### `public `[`astro_status_t`](#astronomy_8h_1a3605f063314f858dec26db307c4e8214)` `[`status`](#structastro__vector__t_1a56282eb2bcf0dce0170952f9d1aa8b3c) 

#### `public double `[`x`](#structastro__vector__t_1a9519f6730bb6e466b810c9b5b053d7a5) 

#### `public double `[`y`](#structastro__vector__t_1a9f78a5f7f0320b6f5194e1006488baa5) 

#### `public double `[`z`](#structastro__vector__t_1aa96463d4940475e13e1df0fdaabbb955) 

#### `public `[`astro_time_t`](#structastro__time__t)` `[`t`](#structastro__vector__t_1ad90a34927631a1940ae0378574d42103) 

# struct `context_peak_altitude_t` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public `[`astro_body_t`](#astronomy_8h_1abc11a4185b6fd5086da5ee9d3c7ce928)` `[`body`](#structcontext__peak__altitude__t_1aa47c0f204bd7230b55dd7d42a5a4d39d) | 
`public int `[`direction`](#structcontext__peak__altitude__t_1a500c08071e6eb2916623c1a235e26412) | 
`public `[`astro_observer_t`](#structastro__observer__t)` `[`observer`](#structcontext__peak__altitude__t_1a2482ce7ba3aac1be3636a2f07a1c8b53) | 
`public double `[`body_radius_au`](#structcontext__peak__altitude__t_1aea0fe57205509ba0753730e84c9e8b97) | 

## Members

#### `public `[`astro_body_t`](#astronomy_8h_1abc11a4185b6fd5086da5ee9d3c7ce928)` `[`body`](#structcontext__peak__altitude__t_1aa47c0f204bd7230b55dd7d42a5a4d39d) 

#### `public int `[`direction`](#structcontext__peak__altitude__t_1a500c08071e6eb2916623c1a235e26412) 

#### `public `[`astro_observer_t`](#structastro__observer__t)` `[`observer`](#structcontext__peak__altitude__t_1a2482ce7ba3aac1be3636a2f07a1c8b53) 

#### `public double `[`body_radius_au`](#structcontext__peak__altitude__t_1aea0fe57205509ba0753730e84c9e8b97) 

# struct `deltat_entry_t` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public double `[`mjd`](#structdeltat__entry__t_1aaea1f67366018adb1df1b8ba59d551b8) | 
`public double `[`dt`](#structdeltat__entry__t_1a2e9ba25a4ea80205f1eabe7878f7c9f8) | 

## Members

#### `public double `[`mjd`](#structdeltat__entry__t_1aaea1f67366018adb1df1b8ba59d551b8) 

#### `public double `[`dt`](#structdeltat__entry__t_1a2e9ba25a4ea80205f1eabe7878f7c9f8) 

# struct `earth_tilt_t` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public double `[`tt`](#structearth__tilt__t_1a672d10b8a865c97d10f940f850e23218) | 
`public double `[`dpsi`](#structearth__tilt__t_1aecf0599dd0aaf72c6da8dee770e22154) | 
`public double `[`deps`](#structearth__tilt__t_1a0d522907b75ba96ec4878e4d3896e881) | 
`public double `[`ee`](#structearth__tilt__t_1a174d39ec59199e212ce12eed15770e15) | 
`public double `[`mobl`](#structearth__tilt__t_1aa87141a38e1a48ed0cea091d662cec61) | 
`public double `[`tobl`](#structearth__tilt__t_1a1220c4af0e00156f1d2d527b4c9176f9) | 

## Members

#### `public double `[`tt`](#structearth__tilt__t_1a672d10b8a865c97d10f940f850e23218) 

#### `public double `[`dpsi`](#structearth__tilt__t_1aecf0599dd0aaf72c6da8dee770e22154) 

#### `public double `[`deps`](#structearth__tilt__t_1a0d522907b75ba96ec4878e4d3896e881) 

#### `public double `[`ee`](#structearth__tilt__t_1a174d39ec59199e212ce12eed15770e15) 

#### `public double `[`mobl`](#structearth__tilt__t_1aa87141a38e1a48ed0cea091d662cec61) 

#### `public double `[`tobl`](#structearth__tilt__t_1a1220c4af0e00156f1d2d527b4c9176f9) 

# struct `MoonContext` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public double `[`t`](#structMoonContext_1aaaf561e8c74cf5b38efeabf945fbe7c7) | 
`public double `[`dgam`](#structMoonContext_1a72984bcb2447cb83383d440afde31f7a) | 
`public double `[`dlam`](#structMoonContext_1aa83859066b6b745f6a865ed61ae9fb0b) | 
`public double `[`n`](#structMoonContext_1a6c440a84e40d26cea2fc0cc2d2b074e4) | 
`public double `[`gam1c`](#structMoonContext_1a591cb67157db09b81e73c208fafb0f30) | 
`public double `[`sinpi`](#structMoonContext_1aaf3e242c96281538ecbae825d158a095) | 
`public double `[`l0`](#structMoonContext_1a124c9586686ae0506029ab98ce4993fe) | 
`public double `[`l`](#structMoonContext_1ac6b65ca94db24548276265c7f75f7ccb) | 
`public double `[`ls`](#structMoonContext_1aa4d999955bb2bca04d25eb00291a5308) | 
`public double `[`f`](#structMoonContext_1a9ef4eefb9f624a98f2823f90c1928417) | 
`public double `[`d`](#structMoonContext_1ab3cc1fec2e00544413ccfcf313bd1d85) | 
`public double `[`s`](#structMoonContext_1a7a2ebea092b91340db5d0146e763f17b) | 
`public double `[`dl0`](#structMoonContext_1af46590adaa21bf6b42fd3783f7d8db9c) | 
`public double `[`dl`](#structMoonContext_1a86bec90c26afe7c88475bd9131a5b392) | 
`public double `[`dls`](#structMoonContext_1a2b17842eef72ebc0e8084afcb18805a7) | 
`public double `[`df`](#structMoonContext_1a8726d33dec64f105c8bfd47ae8845227) | 
`public double `[`dd`](#structMoonContext_1aaee8afdea6e3275079da08934d6b4497) | 
`public double `[`ds`](#structMoonContext_1ae4f02d1ae1f84d3b53a79a2b11f24095) | 
`public  `[`DECLARE_PASCAL_ARRAY_2`](#structMoonContext_1a56219dd4d574ba8f969a346be8e8a94b)`(double,co,- 6,6,1,4)` | 
`public  `[`DECLARE_PASCAL_ARRAY_2`](#structMoonContext_1a6bbd892350493626fccf00c5b68cbd74)`(double,si,- 6,6,1,4)` | 

## Members

#### `public double `[`t`](#structMoonContext_1aaaf561e8c74cf5b38efeabf945fbe7c7) 

#### `public double `[`dgam`](#structMoonContext_1a72984bcb2447cb83383d440afde31f7a) 

#### `public double `[`dlam`](#structMoonContext_1aa83859066b6b745f6a865ed61ae9fb0b) 

#### `public double `[`n`](#structMoonContext_1a6c440a84e40d26cea2fc0cc2d2b074e4) 

#### `public double `[`gam1c`](#structMoonContext_1a591cb67157db09b81e73c208fafb0f30) 

#### `public double `[`sinpi`](#structMoonContext_1aaf3e242c96281538ecbae825d158a095) 

#### `public double `[`l0`](#structMoonContext_1a124c9586686ae0506029ab98ce4993fe) 

#### `public double `[`l`](#structMoonContext_1ac6b65ca94db24548276265c7f75f7ccb) 

#### `public double `[`ls`](#structMoonContext_1aa4d999955bb2bca04d25eb00291a5308) 

#### `public double `[`f`](#structMoonContext_1a9ef4eefb9f624a98f2823f90c1928417) 

#### `public double `[`d`](#structMoonContext_1ab3cc1fec2e00544413ccfcf313bd1d85) 

#### `public double `[`s`](#structMoonContext_1a7a2ebea092b91340db5d0146e763f17b) 

#### `public double `[`dl0`](#structMoonContext_1af46590adaa21bf6b42fd3783f7d8db9c) 

#### `public double `[`dl`](#structMoonContext_1a86bec90c26afe7c88475bd9131a5b392) 

#### `public double `[`dls`](#structMoonContext_1a2b17842eef72ebc0e8084afcb18805a7) 

#### `public double `[`df`](#structMoonContext_1a8726d33dec64f105c8bfd47ae8845227) 

#### `public double `[`dd`](#structMoonContext_1aaee8afdea6e3275079da08934d6b4497) 

#### `public double `[`ds`](#structMoonContext_1ae4f02d1ae1f84d3b53a79a2b11f24095) 

#### `public  `[`DECLARE_PASCAL_ARRAY_2`](#structMoonContext_1a56219dd4d574ba8f969a346be8e8a94b)`(double,co,- 6,6,1,4)` 

#### `public  `[`DECLARE_PASCAL_ARRAY_2`](#structMoonContext_1a6bbd892350493626fccf00c5b68cbd74)`(double,si,- 6,6,1,4)` 

# struct `vsop_formula_t` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public int `[`nseries`](#structvsop__formula__t_1a56afa21306a89bf18ddfac5ca85ee677) | 
`public const `[`vsop_series_t`](#structvsop__series__t)` * `[`series`](#structvsop__formula__t_1acf5e628306f23d4583199dbe2048e5ce) | 

## Members

#### `public int `[`nseries`](#structvsop__formula__t_1a56afa21306a89bf18ddfac5ca85ee677) 

#### `public const `[`vsop_series_t`](#structvsop__series__t)` * `[`series`](#structvsop__formula__t_1acf5e628306f23d4583199dbe2048e5ce) 

# struct `vsop_model_t` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public const `[`vsop_formula_t`](#structvsop__formula__t)` `[`formula`](#structvsop__model__t_1afc3a605659bf015b0439ee087419af94) | 

## Members

#### `public const `[`vsop_formula_t`](#structvsop__formula__t)` `[`formula`](#structvsop__model__t_1afc3a605659bf015b0439ee087419af94) 

# struct `vsop_series_t` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public int `[`nterms`](#structvsop__series__t_1a7aeb051802ffd6ef6ae02a9403031b58) | 
`public const `[`vsop_term_t`](#structvsop__term__t)` * `[`term`](#structvsop__series__t_1a4adb390ffedd92b0726b125c04c9a76f) | 

## Members

#### `public int `[`nterms`](#structvsop__series__t_1a7aeb051802ffd6ef6ae02a9403031b58) 

#### `public const `[`vsop_term_t`](#structvsop__term__t)` * `[`term`](#structvsop__series__t_1a4adb390ffedd92b0726b125c04c9a76f) 

# struct `vsop_term_t` 

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
`public double `[`amplitude`](#structvsop__term__t_1aecd5a26eb7821a3c526b439790d99aab) | 
`public double `[`phase`](#structvsop__term__t_1a67104f91069ca921319a4035cb836ec8) | 
`public double `[`frequency`](#structvsop__term__t_1a82d3dd5cc631234a6185c01c33930532) | 

## Members

#### `public double `[`amplitude`](#structvsop__term__t_1aecd5a26eb7821a3c526b439790d99aab) 

#### `public double `[`phase`](#structvsop__term__t_1a67104f91069ca921319a4035cb836ec8) 

#### `public double `[`frequency`](#structvsop__term__t_1a82d3dd5cc631234a6185c01c33930532) 

Generated by [Moxygen](https://github.com/sourcey/moxygen)