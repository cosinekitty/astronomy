import { ChangeDetectionStrategy, Component, Input } from '@angular/core';
import { SidenavItem } from '../../interfaces/sidenav-item.interface';

@Component({
  selector: 'ae-layout',
  templateUrl: './layout.component.html',
  styleUrls: ['./layout.component.scss'],
  changeDetection: ChangeDetectionStrategy.OnPush
})
export class LayoutComponent {
  @Input() menu!: SidenavItem[];
}
