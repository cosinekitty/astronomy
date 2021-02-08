import { ChangeDetectionStrategy, Component } from '@angular/core';

@Component({
  selector: 'ae-demos-layout',
  templateUrl: './layout.component.html',
  styleUrls: ['./layout.component.scss'],
  changeDetection: ChangeDetectionStrategy.OnPush
})
export class DemosLayoutComponent {
  menu = [
    {
      text: 'Body Positions',
      link: './body-positions',
    },
  ];
}
