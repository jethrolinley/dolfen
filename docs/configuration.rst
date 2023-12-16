==============================
Standard sampler configuration
==============================

There are various settings in ``nessai`` which can be configured. These can be grouped in to general settings and proposal settings. The former controls general aspects of the sampler such as the model being sampler or how many live points are used. The latter affect the proposal process and how new points are drawn.

All of the settings are controlled when creating an instance of :py:class:`~nessai.flowsampler.FlowSampler`. The most important settings to consider when using ``nessai`` are the :doc:`reparameterisations<reparameterisations>` used for the proposals.

General configuration
=====================

These are general settings which apply to the whole algorithm and are parsed to :py:class:`~dolfen.likelihood`. However some of these settings, such as :code:`training_frequency` which defines how often the proposal method is retrained, will affect the normalising flow used in the proposal class.

.. autoapiclass:: dolfen.likelihood
    :members: None

