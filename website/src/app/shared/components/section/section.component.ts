import { ChangeDetectionStrategy, Component, Input } from '@angular/core';
import { SidenavItem } from '../../interfaces/sidenav-item.interface';

@Component({
  selector: 'ae-section',
  templateUrl: './section.component.html',
  styleUrls: ['./section.component.scss'],
  changeDetection: ChangeDetectionStrategy.OnPush
})
export class SectionComponent {
  @Input() title!: string;
  @Input() children?: SidenavItem[];
  activeLink?: string;
}
