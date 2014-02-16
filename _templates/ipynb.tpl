---
layout: lesson
root: ../..
---
{% extends 'markdown.tpl' %}

{% block in_prompt %}
{% endblock in_prompt %}

{% block input %}
~~~
{{ cell.input }}
~~~
{% endblock input %}

{% block markdowncell scoped %}
{% if 'cell_tags' in cell.metadata %}
<div class="{{ cell.metadata['cell_tags'][0] }}">
{{ cell.source  | markdown2html | strip_files_prefix }}
</div>
{% else %}
<div>
{{ cell.source  | markdown2html | strip_files_prefix }}
</div>
{% endif %}
{%- endblock markdowncell %}
