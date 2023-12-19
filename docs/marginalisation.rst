===============
Marginalisation
===============


Gravitational Wave Phase Marginalisation
----------------------------------------

In the :doc:`bilby example<bilby-example>`, the time-domain waveform from `this paper <https://arxiv.org/abs/2004.08302>`_ is coded up (using aligned-spins) and used as the waveform model, and posteriors are produced with 4 free parameters. However, the phase is marginalised out of the likelihood in the 'pre-sampling stage', giving the sampler one less parameter to search over and thus reducing the cost of computing the phase-marginalised posterior by a factor of a few (compared to marginalising post-sampling). ``dolfen`` has an inbuilt (log-)likelihood function that is designed to receive individually the plus and cross polarisations of the the waveform from the ``GWmodel_margphase`` callable, and is activated when a ``GWmodel_margphase`` callable is passed as a parameter to ``dolfen``. This allows for a quick and easy pre-sampling numerical marginalisation/integration over the phase. To use a different model with a (pre-sampler) numerical marginalisation, see :ref:`marginalisations for arbitrary models<arbmarg>`.

.. _arbmarg:
Marginalisations For Arbitrary Models
-------------------------------------


.. note::
    This section has not been completed yet 
