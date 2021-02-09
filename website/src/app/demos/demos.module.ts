import { CommonModule } from '@angular/common';
import { NgModule } from '@angular/core';
import { FlexLayoutModule } from '@angular/flex-layout';
import { ReactiveFormsModule } from '@angular/forms';
import { RouterModule, Routes } from '@angular/router';
import { SharedModule } from '../shared/shared.module';
import { BodyPositionsComponent } from './components/body-positions/body-positions.component';
import { DemosComponent } from './containers/demos/demos.component';

const routes: Routes = [
  {
    path: '',
    component: DemosComponent,
    children: [
      {
        path: '',
        pathMatch: 'full',
        redirectTo: 'body-positions',
      },
      {
        path: 'body-positions',
        component: BodyPositionsComponent,
      },
    ]
  }
];

@NgModule({
  imports: [
    CommonModule,
    FlexLayoutModule,
    ReactiveFormsModule,
    RouterModule.forChild(routes),
    SharedModule,
  ],
  declarations: [
    DemosComponent,
    BodyPositionsComponent,
  ]
})
export class DemosModule {}
