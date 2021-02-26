import { ChangeDetectionStrategy, Component, OnInit } from '@angular/core';
import { SidenavItem } from '../../../shared/interfaces/sidenav-item.interface';

@Component({
  selector: 'ae-coordinate-systems',
  templateUrl: './coordinate-systems.component.html',
  styleUrls: ['./coordinate-systems.component.scss'],
  changeDetection: ChangeDetectionStrategy.OnPush
})
export class CoordinateSystemsComponent implements OnInit {
  children: SidenavItem[] = [
    {
      text: 'Equatorial',
    },
    {
      text: 'Ecliptic',
    },
    {
      text: 'Horizontal',
    },
  ];

  ngOnInit(): void {
    this.onChange(0);
  }

  onChange(index: number): void {
    // console.log(this.children[index]);
  }
}

