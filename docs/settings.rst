========
Settings
========

Various settings in ``dolfen`` can be configured, although the default settings should be acceptable for use with systems that roughly comply with the description of the ideal sort of system dolfen is designed for (i.e., slow evolution signals consisting of a large number of datapoints). Some understanding of the background theory is required to know how changing the settings could affect the accuracy of dolfen's representation of the likelihood, for example, uniform downsampling can be expected to lead to aliasing.

All settings are passed when creating an instance of :py:class:`~dolfen.likelihood`. The most important settings are the number of samples, the downsampling scheme (how to select remaining datapoints) and method (approximate FIM preservation, or using FIM Jeffreys' divergence minimisation), and the MCS override - see more details in the API below.

API and standard configuration
==============================

These are general settings which apply to the whole algorithm and are parsed to :py:class:`~dolfen.likelihood`. However some of these settings, such as :code:`prsrv_FIM` which defines what method is used, will affect the ability to find solutions given the number of samples chosen.


.. autoapiclass:: dolfen.likelihood
    :members: None

