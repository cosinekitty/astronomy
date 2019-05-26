# {{kind}} `{{name}}` {{anchor refid}}

{{briefdescription}}

{{detaileddescription}}

{{#if filtered.members}}

## Summary

 Members                        | Descriptions                                
--------------------------------|---------------------------------------------
{{#each filtered.members}}{{cell proto}}            | {{cell summary}}
{{/each}}{{#each filtered.compounds}}{{cell proto}} | {{cell summary}}
{{/each}}

## Members

{{#each filtered.members}}
#### {{title proto}} {{anchor refid}}

{{#if enumvalue}}
 Values                         | Descriptions                                
--------------------------------|---------------------------------------------
{{#each enumvalue}}{{cell name}}            | {{cell summary}}
{{/each}}
{{/if}}

{{briefdescription}}

{{detaileddescription}}

{{/each}}
{{/if}}
