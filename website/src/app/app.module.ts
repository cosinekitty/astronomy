import { NgModule } from '@angular/core';
import { BrowserAnimationsModule } from '@angular/platform-browser/animations';
import { RouterModule, Routes } from '@angular/router';
import { AppComponent } from './app.component';
import { MainModule } from './main/main.module';

const routes: Routes = [];

@NgModule({
  imports: [
    BrowserAnimationsModule,
    RouterModule.forRoot(routes),
    MainModule,
  ],
  declarations: [AppComponent],
  bootstrap: [AppComponent]
})
export class AppModule { }
