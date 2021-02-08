import { ChangeDetectionStrategy, Component } from '@angular/core';

@Component({
  selector: 'ae-theme',
  templateUrl: './theme.component.html',
  styleUrls: ['./theme.component.scss'],
  changeDetection: ChangeDetectionStrategy.OnPush
})
export class ThemeComponent {}
