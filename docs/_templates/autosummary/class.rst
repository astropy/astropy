{% if referencefile %}
.. include:: {{ referencefile }}
{% endif %}

{{ objname }}
{{ underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}

   {% block attributes %}
   {% if attributes %}

   .. rubric:: Attributes

   .. autosummary::
   {% for item in attributes %}
      ~{{ name }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block methods %}

   {% if methods %}

   .. rubric:: Methods Summary

   .. autosummary::
   {% for item in methods %}
      {% if item != "__init__" %}
      ~{{ name }}.{{ item }}
      {% endif %}
   {%- endfor %}

   .. rubric:: Methods Documentation

   {% for item in methods %}
      {% if item != "__init__" %}
   .. automethod:: {{ item }}
      {% endif %}
   {%- endfor %}
   {% endif %}
   {% endblock %}



