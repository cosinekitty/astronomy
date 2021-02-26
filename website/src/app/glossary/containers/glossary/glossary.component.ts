import { ChangeDetectionStrategy, Component } from '@angular/core';

@Component({
  selector: 'ae-glossary',
  templateUrl: './glossary.component.html',
  styleUrls: ['./glossary.component.scss'],
  changeDetection: ChangeDetectionStrategy.OnPush
})
export class GlossaryComponent {
  menu = [
    {
      text: 'Observers',
      link: './observers',
    },
    {
      text: 'Coordinate Systems',
      link: './coordinate-systems',
    },
    {
      text: 'Time Scales',
      link: './time-scales',
    },
    {
      text: 'Barycenters',
      link: './barycenters',
    },
    {
      text: 'Apsides',
      link: './apsides',
    },
    {
      text: 'Equinox',
      link: './equinox',
    },
  ];
}
