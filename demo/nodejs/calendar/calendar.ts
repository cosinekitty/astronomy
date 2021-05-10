/*
    calendar.ts   -   Don Cross   -   2021-05-09

    A demo of using Astronomy Engine to find a series of
    interesting events for a calendar.
*/

import { AstroTime } from "./astronomy";


function RunTest(): void {
    const now = new AstroTime(new Date());
    console.log(now);
}


RunTest();
