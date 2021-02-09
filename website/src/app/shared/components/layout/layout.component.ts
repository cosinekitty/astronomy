import { ChangeDetectionStrategy, Component, HostBinding, Input } from '@angular/core';
import { SidenavItem } from '../../interfaces/sidenav-item.interface';

@Component({
  selector: 'ae-layout',
  templateUrl: './layout.component.html',
  styleUrls: ['./layout.component.scss'],
  changeDetection: ChangeDetectionStrategy.OnPush
})
export class LayoutComponent {
  @HostBinding('class')
  @Input() color!: string;
  @Input() section!: string;
  @Input() icon!: string;
  @Input() menu!: SidenavItem[];
}
