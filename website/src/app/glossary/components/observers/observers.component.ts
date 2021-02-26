import { ChangeDetectionStrategy, Component, OnInit } from '@angular/core';
import { SidenavItem } from '../../../shared/interfaces/sidenav-item.interface';

@Component({
  selector: 'ae-observers',
  templateUrl: './observers.component.html',
  styleUrls: ['./observers.component.scss'],
  changeDetection: ChangeDetectionStrategy.OnPush
})
export class ObserversComponent implements OnInit {
  children: SidenavItem[] = [
    {
      text: 'Heliocentric',
    },
    {
      text: 'Geocentric',
    },
    {
      text: 'Topocentric',
    },
  ];

  ngOnInit(): void {
    this.onChange(0);
  }

  onChange(index: number): void {
    // console.log(this.children[index]);
  }
}
