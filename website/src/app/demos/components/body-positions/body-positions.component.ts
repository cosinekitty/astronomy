import { ChangeDetectionStrategy, ChangeDetectorRef, Component, OnInit } from '@angular/core';
import { FormBuilder, FormGroup } from '@angular/forms';
import { UntilDestroy, untilDestroyed } from '@ngneat/until-destroy';
import { Equator, Horizon, MakeObserver } from 'astronomy-engine';
import { combineLatest, interval } from 'rxjs';

interface BodyPositionsArgs {
  datetime: string;
  livetime: boolean;
  lat: number;
  lng: number;
  alt: number;
}

interface BodyPositionsData {
  body: string;
  ra: string;
  dec: string;
  az: string;
  alt: string;
}

@UntilDestroy()
@Component({
  selector: 'ae-demos-body-positions',
  templateUrl: './body-positions.component.html',
  styleUrls: ['./body-positions.component.scss'],
  changeDetection: ChangeDetectionStrategy.OnPush,
})
export class BodyPositionsComponent implements OnInit {
  form!: FormGroup;

  bodies = [
    'Sun', 'Moon', 'Mercury', 'Venus', 'Mars',
    'Jupiter', 'Saturn', 'Uranus', 'Neptune', 'Pluto',
  ];

  columns: string[] = ['body', 'ra', 'dec', 'az', 'alt'];

  data: BodyPositionsData[] = [];

  constructor(private cdr: ChangeDetectorRef, private fb: FormBuilder) {}

  ngOnInit(): void {
    this.form  = this.fb.group({
      datetime: this.formatDate(new Date()),
      livetime: true,
      lat: 0,
      lng: 0,
      alt: 0,
    });

    this.getCurrentLocation().then(pos => this.form.patchValue(pos));

    // listen changes
    combineLatest([
      interval(1000),
      this.form.valueChanges,
    ])
    .pipe(untilDestroyed(this))
    .subscribe(([, values]: [number, BodyPositionsArgs]) => {
      if (values.livetime) {
        this.updateTable(values);
      }
    });
  }

  getCurrentLocation(): Promise<any> {
    return new Promise((resolve, reject) => {
      navigator.geolocation.getCurrentPosition(
        resp => resolve({
          lat: resp.coords.latitude,
          lng: resp.coords.longitude,
        }),
        err => reject(err)
      );
    });
  }

  updateTable({ datetime, livetime, lat, lng, alt }: BodyPositionsArgs): void {
    let date: Date;
    if (livetime) {
      date = new Date();
      this.form.patchValue({ datetime: this.formatDate(date) }, { emitEvent: false });
    } else {
      date = new Date(datetime);
    }

    this.data = [];

    try {
      const observer = MakeObserver(Number(lat), Number(lng), Number(alt));

      this.bodies.forEach(body => {
        const equ2000 = Equator(body, date, observer, false, true);
        const equofdate = Equator(body, date, observer, true, true);
        const hor = Horizon(date, observer, equofdate.ra, equofdate.dec, 'normal');

        this.data.push({
          body,
          ra: equ2000.ra.toFixed(2),
          dec: equ2000.dec.toFixed(2),
          az: hor.azimuth.toFixed(2),
          alt: hor.altitude.toFixed(2),
        });
      });

      this.cdr.detectChanges();

    } catch (e) {
      console.error(e);
    }
  }

  private pad(val: number, len: number): string {
    return val.toFixed(0).padStart(len, '0');
  }

  private formatDate(date: Date): string {
    const year = this.pad(date.getFullYear(), 4);
    const month = this.pad(1 + date.getMonth(), 2);
    const day = this.pad(date.getDate(), 2);
    const hour = this.pad(date.getHours(), 2);
    const minute = this.pad(date.getMinutes(), 2);
    const second = this.pad(date.getSeconds(), 2);
    return `${year}-${month}-${day} ${hour}:${minute}:${second}`;
  }
}
