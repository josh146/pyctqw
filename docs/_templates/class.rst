{{ name }}
{{ underline }}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
	:show-inheritance:

{% block methods %}

Method Summary
~~~~~~~~~~~~~~

.. autosummary::
	{% for item in methods %}
	{% if item != '__init__' %}
	{{ name }}.{{ item }}
	{% endif %}
	{%- endfor %}

Method Details
~~~~~~~~~~~~~~

{% for item in methods %}
{% if item != '__init__' %}
.. automethod::	{{ name }}.{{ item }}

{% endif %}
{%- endfor %}

{% endblock %}
