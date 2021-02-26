import { ChangeDetectionStrategy, Component } from '@angular/core';

@Component({
  selector: 'ae-main-index',
  templateUrl: './index.component.html',
  styleUrls: ['./index.component.scss'],
  changeDetection: ChangeDetectionStrategy.OnPush
})
export class IndexComponent {
  items = [
    {
      link: '/glossary',
      icon: 'description',
      color: 'blue',
      title: 'Glossary',
      description: 'Concepts and terminology',
      disabled: false,
    },
    {
      link: '/demos',
      icon: 'travel_explore',
      color: 'green',
      title: 'Demos',
      description: 'See AstronomyEngine in action',
      disabled: false,
    },
    {
      link: '/showcase',
      icon: 'campaign',
      color: 'pink',
      title: 'Showcase',
      description: 'Applications using AstronomyEngine',
      disabled: true,
    },
    {
      link: '/start',
      icon: 'assistant',
      color: 'purple',
      title: 'Quick Start',
      description: 'Requirements and Installation',
      disabled: true,
    },
    {
      link: '/documentation',
      icon: 'extension',
      color: 'aqua',
      title: 'API',
      description: 'Documentation',
      disabled: true,
    },
    {
      link: '/credits',
      icon: 'volunteer_activism',
      color: 'orange',
      title: 'Credits',
      description: 'Author and contributors',
      disabled: true,
    },
  ];
}
