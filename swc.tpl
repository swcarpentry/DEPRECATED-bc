---
layout: lesson
root: ../..
---
{% extends 'markdown.tpl' %}

{% block input %}
<div class="input">
<span>In [{{ cell.prompt_number }}]:</span>
<pre>{{ cell.input }}</pre>
</div>
{% endblock input %}

{% block output_group %}
<div class="output">
<span>Out [{{ cell.prompt_number }}]:</span>
<pre>{{ super() }}</pre>
</div>
{% endblock output_group %}

{% block stream %}
{{ output.text }}
{% endblock stream %}

{% block pyout %}
{{ output.text }}
{% endblock pyout %}

{% block pyerr %}
{{ output.traceback | join('\n') | strip_ansi | escape }}
{% endblock pyerr %}

{% block markdowncell %}
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
