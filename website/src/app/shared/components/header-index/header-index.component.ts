import { ChangeDetectionStrategy, Component, HostBinding } from '@angular/core';

@Component({
  selector: 'ae-header-index',
  templateUrl: './header-index.component.html',
  styleUrls: ['./header-index.component.scss'],
  changeDetection: ChangeDetectionStrategy.OnPush
})
export class HeaderIndexComponent {
  @HostBinding('class.header') isHeader = true;
}
