================
Marginalisations
================


Gravitational Wave Phase Marginalisation
========================================

In the :doc:`bilby example<bilby-example>`, the time-domain waveform from `this paper<https://arxiv.org/abs/2004.08302>_` is coded up (using aligned-spins) and used as the waveform model, with posteriors produce with 4 free parameters. However, the phase is marginalised out of the likelihood in the 'pre-sampling' phase, giving the sampler one less parameter to search over and thus reducing the cost of produing the posterior by a factor of a few. ``dolfen`` has an inbuilt 
