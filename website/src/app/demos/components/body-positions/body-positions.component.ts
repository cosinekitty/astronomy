import { ChangeDetectionStrategy, Component, OnInit } from '@angular/core';
import { FormBuilder } from '@angular/forms';
import { Bodies, Equator, Horizon, MakeObserver } from 'astronomy-engine';

interface BodyPositionsArgs {
  date: Date;
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

@Component({
  selector: 'ae-demos-body-positions',
  templateUrl: './body-positions.component.html',
  styleUrls: ['./body-positions.component.scss'],
  changeDetection: ChangeDetectionStrategy.OnPush,
})
export class BodyPositionsComponent implements OnInit {
  form = this.fb.group({
    date: new Date(),
    lat: 0,
    lng: 0,
    alt: 0,
  });

  columns: string[] = ['body', 'ra', 'dec', 'az', 'alt'];

  data: BodyPositionsData[] = [];

  constructor(private fb: FormBuilder) {}

  ngOnInit(): void {
    // TODO format and parse string date
    this.getCurrentLocation().then(pos => this.form.patchValue(pos));

    this.form.valueChanges.subscribe(values => this.updateTable(values));
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

  updateTable({ date, lat, lng, alt }: BodyPositionsArgs): void {
    // TODO validate inputs
    this.data = [];
    const observer = MakeObserver(lat, lng, alt);

    Bodies.forEach(body => {
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
  }
}
