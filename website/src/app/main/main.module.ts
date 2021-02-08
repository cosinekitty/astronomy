import { CommonModule } from '@angular/common';
import { NgModule } from '@angular/core';
import { FlexLayoutModule } from '@angular/flex-layout';
import { MatIconModule } from '@angular/material/icon';
import { MatTooltipModule } from '@angular/material/tooltip';
import { RouterModule, Routes } from '@angular/router';
import { ThemeComponent } from '../shared/components/theme/theme.component';
import { SharedModule } from '../shared/shared.module';
import { IndexComponent } from './components/index/index.component';

const routes: Routes = [
  {
    path: '',
    component: ThemeComponent,
    children: [
      {
        path: '',
        component: IndexComponent,
      }
    ]
  }
];

@NgModule({
  imports: [
    CommonModule,
    FlexLayoutModule,
    MatIconModule,
    MatTooltipModule,
    RouterModule.forChild(routes),
    SharedModule,
  ],
  declarations: [
    IndexComponent,
  ],
  exports: [
    RouterModule,
  ],
})
export class MainModule { }
