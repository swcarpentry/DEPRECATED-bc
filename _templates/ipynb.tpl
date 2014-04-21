---
layout: lesson
root: ../..
---
{% extends 'markdown.tpl' %}

{% block input %}
<div class="in">
<pre>{{ cell.input | escape }}</pre>
</div>
{% endblock input %}

{% block output_group %}
<div class="out">
<pre>{{- super() -}}</pre>
</div>
{% endblock output_group %}

{%- block stream -%}{{- output.text | escape -}}{%- endblock stream -%}

{%- block pyout -%}{{- output.text | escape -}}{%- endblock pyout -%}

{%- block pyerr -%}{{- output.traceback | join('\n') | strip_ansi | escape -}}{%- endblock pyerr -%}

{%- block data_svg -%}<img src="../../{{ output.svg_filename | path2url }}">{%- endblock data_svg -%}

{%- block data_png -%}<img src="../../{{ output.png_filename | path2url }}">{%- endblock data_png -%}

{%- block data_jpg -%}<img src="../../{{ output.jpeg_filename | path2url }}">{%- endblock data_jpg -%}

{%- block data_text -%}{{ output.html }}{%- endblock data_text -%}

{%- block display_data -%}
{%- if output.html -%}
{{- output.html -}}
{%- else -%}
{{- super() -}}
{%- endif -%}
{%- endblock display_data -%}

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
