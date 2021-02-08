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
      link: '/demos',
      icon: 'travel_explore',
      color: 'blue',
      title: 'Demos',
      description: 'See AstronomyEngine in action',
      disabled: false,
    },
    {
      link: '/start',
      icon: 'assistant',
      color: 'green',
      title: 'Quick Start',
      description: 'Requirements and Installation',
      disabled: true,
    },
    {
      link: '/documentation',
      icon: 'extension',
      color: 'pink',
      title: 'API',
      description: 'Documentation',
      disabled: true,
    },
    {
      link: '/faq',
      icon: 'support',
      color: 'purple',
      title: 'FAQs',
      description: 'Frequently asked questions',
      disabled: true,
    },
    {
      link: '/showcase',
      icon: 'campaign',
      color: 'aqua',
      title: 'Showcase',
      description: 'Some software using AstronomyEngine',
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
