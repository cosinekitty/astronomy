import { ChangeDetectionStrategy, Component, Input } from '@angular/core';
import { SidenavItem } from '../../interfaces/sidenav-item.interface';

@Component({
  selector: 'ae-sidenav',
  templateUrl: './sidenav.component.html',
  styleUrls: ['./sidenav.component.scss'],
  changeDetection: ChangeDetectionStrategy.OnPush
})
export class SidenavComponent {
  @Input() menu!: SidenavItem[];
}
