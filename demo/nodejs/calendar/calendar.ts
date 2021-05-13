/*
    calendar.ts   -   Don Cross   -   2021-05-09

    A demo of using Astronomy Engine to find a series of
    interesting events for a calendar.
*/

import {
    AstroTime,
    Body,
    Observer,
    SearchMoonQuarter, NextMoonQuarter, MoonQuarter,
    SearchRiseSet,
    Seasons,
} from "./astronomy";


class AstroEvent {
    constructor(
        public time: AstroTime,
        public title: string,
        public creator: AstroEventEnumerator)
        {}
}


interface AstroEventEnumerator {
    FindFirst(startTime: AstroTime): AstroEvent;
    FindNext(): AstroEvent;
}


class EventCollator implements AstroEventEnumerator {
    private eventQueue: AstroEvent[];

    constructor(private enumeratorList: AstroEventEnumerator[]) {
    }

    FindFirst(startTime: AstroTime): AstroEvent {
        this.eventQueue = [];
        for (let enumerator of this.enumeratorList)
            this.InsertEvent(enumerator.FindFirst(startTime));
        return this.FindNext();
    }

    FindNext(): AstroEvent {
        if (this.eventQueue.length === 0)
            return null;

        const evt = this.eventQueue.shift();
        const another = evt.creator.FindNext();
        this.InsertEvent(another);
        return evt;
    }

    InsertEvent(evt: AstroEvent): void {
        if (evt !== null) {
            // Insert the event in time order -- after anything that happens before it.

            let i = 0;
            while (i < this.eventQueue.length && this.eventQueue[i].time.tt < evt.time.tt)
                ++i;

            this.eventQueue.splice(i, 0, evt);
        }
    }
}


class RiseSetEnumerator implements AstroEventEnumerator {
    private nextSearchTime: AstroTime;

    constructor(private observer: Observer, private body: Body, private direction: number, private title: string) {
    }

    FindFirst(startTime: AstroTime): AstroEvent {
        this.nextSearchTime = SearchRiseSet(this.body, this.observer, this.direction, startTime, 366.0);
        if (this.nextSearchTime)
            return new AstroEvent(this.nextSearchTime, this.title, this);
        return null;
    }

    FindNext(): AstroEvent {
        if (this.nextSearchTime) {
            const startTime = this.nextSearchTime.AddDays(0.01);
            return this.FindFirst(startTime);
        }
        return null;
    }
}


class SeasonEnumerator implements AstroEventEnumerator {
    private slist: AstroEvent[];
    private year: number;
    private index: number;

    FindFirst(startTime: AstroTime): AstroEvent {
        this.year = startTime.date.getUTCFullYear();
        this.LoadYear(this.year);
        while (this.index < this.slist.length && this.slist[this.index].time.tt < startTime.tt)
            ++this.index;
        return this.FindNext();
    }

    FindNext(): AstroEvent {
        if (this.index === this.slist.length)
            this.LoadYear(++this.year);
        return this.slist[this.index++];
    }

    private LoadYear(year: number): void {
        const seasons = Seasons(year);
        this.slist = [
            new AstroEvent(seasons.mar_equinox,  'March equinox', this),
            new AstroEvent(seasons.jun_solstice, 'June solstice', this),
            new AstroEvent(seasons.sep_equinox,  'September equinox', this),
            new AstroEvent(seasons.dec_solstice, 'December solstice', this)
        ];
        this.index = 0;
    }
}


class MoonQuarterEnumerator implements AstroEventEnumerator {
    private mq: MoonQuarter;

    FindFirst(startTime: AstroTime): AstroEvent {
        this.mq = SearchMoonQuarter(startTime);
        return this.MakeEvent();
    }

    FindNext(): AstroEvent {
        this.mq = NextMoonQuarter(this.mq);
        return this.MakeEvent();
    }

    private MakeEvent(): AstroEvent {
        return new AstroEvent(
            this.mq.time,
            ['new moon', 'first quarter', 'full moon', 'third quarter'][this.mq.quarter],
            this
        );
    }
}


function RunTest(): void {
    const startTime = new AstroTime(new Date('2021-05-12T00:00:00Z'));
    const observer = new Observer(28.6, -81.2, 10.0);
    const collator = new EventCollator([
        new RiseSetEnumerator(observer, Body.Sun, +1, 'sunrise'),
        new RiseSetEnumerator(observer, Body.Sun, -1, 'sunset'),
        new RiseSetEnumerator(observer, Body.Moon, +1, 'moonrise'),
        new RiseSetEnumerator(observer, Body.Moon, -1, 'moonset'),
        new SeasonEnumerator(),
        new MoonQuarterEnumerator()
    ]);

    const stopYear = startTime.date.getUTCFullYear() + 11;
    let evt:AstroEvent = collator.FindFirst(startTime);
    while (evt !== null && evt.time.date.getUTCFullYear() < stopYear) {
        console.log(`${evt.time} ${evt.title}`);
        evt = collator.FindNext();
    }
}

RunTest();
