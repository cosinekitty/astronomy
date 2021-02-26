import { HttpClient, HttpClientModule } from '@angular/common/http';
import { NgModule } from '@angular/core';
import { MAT_FORM_FIELD_DEFAULT_OPTIONS } from '@angular/material/form-field';
import { BrowserAnimationsModule } from '@angular/platform-browser/animations';
import { RouterModule, Routes } from '@angular/router';
import { MarkdownModule, MarkedOptions, MarkedRenderer } from 'ngx-markdown';
import { AppComponent } from './app.component';
import { MainModule } from './main/main.module';

const routes: Routes = [
  {
    path: 'demos',
    loadChildren: () => import('./demos/demos.module').then(m => m.DemosModule)
  },
  { path: 'glossary', loadChildren: () => import('./glossary/glossary.module').then(m => m.GlossaryModule) },
];

export function markedOptionsFactory(): MarkedOptions {
  const renderer = new MarkedRenderer();

  renderer.heading = (text: string, level: number) => {
    const escapedText = text.toLowerCase().replace(/[^\w]+/g, '-').substr(0, 20);
    return `<h${ level + 1 }>` +
      `<a name="${ escapedText }" class="anchor">` +
        `<span class="header-link"></span>` +
      `</a>${ text }` +
    `</h${ level + 1 }>`;
  };

  return {
    renderer,
    gfm: true,
    headerIds: true,
    breaks: false,
    pedantic: false,
    smartLists: true,
    smartypants: false,
  };
}

@NgModule({
  imports: [
    HttpClientModule,
    BrowserAnimationsModule,
    RouterModule.forRoot(routes, { useHash: true }),
    MarkdownModule.forRoot({
      loader: HttpClient,
      markedOptions: {
        provide: MarkedOptions,
        useFactory: markedOptionsFactory,
      },
      // sanitize: SecurityContext.NONE,
    }),
    MainModule,
  ],
  declarations: [AppComponent],
  bootstrap: [AppComponent],
  providers: [
    {
      provide: MAT_FORM_FIELD_DEFAULT_OPTIONS,
      useValue: {
        appearance: 'standard',
        floatLabel: 'always'
      }
    }
  ]
})
export class AppModule { }
