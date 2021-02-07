import { CommonModule } from '@angular/common';
import { NgModule } from '@angular/core';
import { MatIconModule } from '@angular/material/icon';
import { RouterModule } from '@angular/router';
import { HeaderIndexComponent } from './components/header-index/header-index.component';
import { HeaderComponent } from './components/header/header.component';
import { ThemeComponent } from './components/theme/theme.component';

@NgModule({
  imports: [
    CommonModule,
    MatIconModule,
    RouterModule,
  ],
  declarations: [
    ThemeComponent,
    HeaderIndexComponent,
    HeaderComponent,
  ],
  exports: [
    ThemeComponent,
    HeaderIndexComponent,
    HeaderComponent,
  ]
})
export class SharedModule { }
