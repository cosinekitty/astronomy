import { CommonModule } from '@angular/common';
import { NgModule } from '@angular/core';
import { FlexLayoutModule } from '@angular/flex-layout';
import { MatCheckboxModule } from '@angular/material/checkbox';
import { MatFormFieldModule } from '@angular/material/form-field';
import { MatIconModule } from '@angular/material/icon';
import { MatInputModule } from '@angular/material/input';
import { MatTableModule } from '@angular/material/table';
import { MatTabsModule } from '@angular/material/tabs';
import { RouterModule } from '@angular/router';
import { HeaderIndexComponent } from './components/header-index/header-index.component';
import { HeaderComponent } from './components/header/header.component';
import { LayoutComponent } from './components/layout/layout.component';
import { SectionComponent } from './components/section/section.component';
import { SidenavComponent } from './components/sidenav/sidenav.component';
import { ThemeComponent } from './components/theme/theme.component';

@NgModule({
  imports: [
    CommonModule,
    FlexLayoutModule,
    MatIconModule,
    MatTabsModule,
    RouterModule,
  ],
  declarations: [
    ThemeComponent,
    HeaderIndexComponent,
    HeaderComponent,
    LayoutComponent,
    SidenavComponent,
    SectionComponent,
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
    MatTabsModule,
    SectionComponent,
  ]
})
export class SharedModule { }
