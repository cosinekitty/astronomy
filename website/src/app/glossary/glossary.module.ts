import { CommonModule } from '@angular/common';
import { NgModule } from '@angular/core';
import { RouterModule, Routes } from '@angular/router';
import { MarkdownModule } from 'ngx-markdown';
import { SharedModule } from '../shared/shared.module';
import { ApsidesComponent } from './components/apsides/apsides.component';
import { BarycentersComponent } from './components/barycenters/barycenters.component';
import { CoordinateSystemsComponent } from './components/coordinate-systems/coordinate-systems.component';
import { EquinoxComponent } from './components/equinox/equinox.component';
import { ObserversComponent } from './components/observers/observers.component';
import { TimeScalesComponent } from './components/time-scales/time-scales.component';
import { GlossaryComponent } from './containers/glossary/glossary.component';

const routes: Routes = [
  {
    path: '',
    component: GlossaryComponent,
    children: [
      {
        path: '',
        pathMatch: 'full',
        redirectTo: 'observers',
      },
      {
        path: 'observers',
        component: ObserversComponent,
      },
      {
        path: 'coordinate-systems',
        component: CoordinateSystemsComponent,
      },
      {
        path: 'time-scales',
        component: TimeScalesComponent,
      },
      {
        path: 'barycenters',
        component: BarycentersComponent,
      },
      {
        path: 'apsides',
        component: ApsidesComponent,
      },
      {
        path: 'equinox',
        component: EquinoxComponent,
      },
    ]
  }
];


@NgModule({
  imports: [
    CommonModule,
    RouterModule.forChild(routes),
    MarkdownModule.forChild(),
    SharedModule,
  ],
  declarations: [
    GlossaryComponent,
    ObserversComponent,
    CoordinateSystemsComponent,
    TimeScalesComponent,
    BarycentersComponent,
    ApsidesComponent,
    EquinoxComponent,
  ],
})
export class GlossaryModule {}
