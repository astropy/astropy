{% if render_title %}
{% if versiondata.name %}
{{ versiondata.name }} {{ versiondata.version }} ({{ versiondata.date }})
{{ top_underline * ((versiondata.name + versiondata.version + versiondata.date)|length + 4)}}
{% else %}
{{ versiondata.version }} ({{ versiondata.date }})
{{ top_underline * ((versiondata.version + versiondata.date)|length + 3)}}
{% endif %}
{% endif %}

{% for category, val in definitions.items() %}

{% set underline = underlines[0] %}
{{ definitions[category]['name'] }}
{{ underline * definitions[category]['name']|length }}
{% set underline = underlines[1] %}

{% for section, _ in sections.items() %}
{% if section and category in sections[section] %}
{{section}}
{{ underline * section|length }}

{% endif %}
{% if sections[section] and category in sections[section] %}
{% if definitions[category]['showcontent'] %}
{% for text, values in sections[section][category].items() %}
- {{ text }} [{{ values|join(', ') }}]

{% endfor %}
{% else %}
- {{ sections[section][category]['']|join(', ') }}

{% endif %}
{% if sections[section][category]|length == 0 %}
No significant changes.

{% else %}
{% endif %}
{% else %}
{# No significant changes. #}
{% endif %}
{% endfor %}
{% endfor %}
