======================
Standard configuration
======================

There are various settings in ``dolfen`` which can be configured. These can be grouped in to general settings and proposal settings. The former controls general aspects of the sampler such as the model being sampler or how many live points are used. The latter affect the proposal process and how new points are drawn.

All of the settings are controlled when creating an instance of :py:class:`~dolfen.likelihood`. The most important settings to consider when using ``dolfen`` are .....

General configuration
=====================

These are general settings which apply to the whole algorithm and are parsed to :py:class:`~dolfen.likelihood`. However some of these settings, such as :code:`prsrv_FIM` which defines what method is used, will affect the ability to find solutions given the number of samples chosen.


.. autoapiclass:: dolfen.likelihood
    :members: None

