import { ChangeDetectionStrategy, Component } from '@angular/core';

@Component({
  selector: 'ae-demos',
  templateUrl: './demos.component.html',
  styleUrls: ['./demos.component.scss'],
  changeDetection: ChangeDetectionStrategy.OnPush
})
export class DemosComponent {
  menu = [
    {
      text: 'Body Positions',
      link: './body-positions',
    },
  ];
}
