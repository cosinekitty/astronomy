import { CommonModule } from '@angular/common';
import { NgModule } from '@angular/core';
import { FlexLayoutModule } from '@angular/flex-layout';
import { MatCheckboxModule } from '@angular/material/checkbox';
import { MatFormFieldModule } from '@angular/material/form-field';
import { MatIconModule } from '@angular/material/icon';
import { MatInputModule } from '@angular/material/input';
import { MatTableModule } from '@angular/material/table';
import { RouterModule } from '@angular/router';
import { HeaderIndexComponent } from './components/header-index/header-index.component';
import { HeaderComponent } from './components/header/header.component';
import { LayoutComponent } from './components/layout/layout.component';
import { SidenavComponent } from './components/sidenav/sidenav.component';
import { ThemeComponent } from './components/theme/theme.component';

@NgModule({
  imports: [
    CommonModule,
    FlexLayoutModule,
    MatIconModule,
    RouterModule,
  ],
  declarations: [
    ThemeComponent,
    HeaderIndexComponent,
    HeaderComponent,
    LayoutComponent,
    SidenavComponent,
  ],
  exports: [
    ThemeComponent,
    HeaderIndexComponent,
    HeaderComponent,
    LayoutComponent,
    SidenavComponent,
    // material modules
    MatCheckboxModule,
    MatFormFieldModule,
    MatIconModule,
    MatInputModule,
    MatTableModule,
  ]
})
export class SharedModule { }
