---
layout: lesson
root: ../..
---
{% extends 'markdown.tpl' %}

{% block input %}
<pre class="in"><code>{{ cell.input | escape }}</code></pre>
{% endblock input %}

{% block output_group %}
<div class="out">{{- super() -}}</div>
{% endblock output_group %}

{%- block stream -%}<pre class='out'><code>{{- output.text | escape -}}</code></pre>{%- endblock stream -%}

{%- block pyout -%}<pre class='out'><code>{{- output.text | escape -}}</code></pre>{%- endblock pyout -%}

{%- block pyerr -%}<pre class='err'><code>{{- output.traceback | join('\n') | strip_ansi | escape -}}</code></pre>{%- endblock pyerr -%}

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
{% if 'cell_tags' in cell.metadata and cell.metadata['cell_tags'] %}
<div class="{{ cell.metadata['cell_tags'][0] }}" markdown="1">
{{ cell.source | strip_files_prefix }}
</div>
{% else %}
{{ cell.source | strip_files_prefix }}
{% endif %}
{%- endblock markdowncell %}
