.. dolfen documentation master file, created by
   sphinx-quickstart on Sat Dec 16 10:42:56 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Welcome to dolfen's documentation!
==================================

======
Dolfen
======

``dolfen``: Downsampling Likelihood Function Estimation

``dolfen`` is an algorithm that infers the likelihood function used in Bayesian parameter estimation (PE) by downsampling, i.e., discarding most of the datapoints, then redefining the likelihood's inner product operator on the remaining datapoints to compensate for lost information. It is thus for use in simulated environments only where one has control of the signal and noise models and is particularly useful for rapid likelihood calculation of long duration, slowly evolving signals consisting of large numbers of datapoints, where the Bayesian likelihood is expensive.



Time-domain Gravitational Wave Inference
----------------------------------------

Although ``dolfen`` may be applied generally, it was created for the specific problem of simulated parameter estimation of low to intermediate mass compact binary coalescences (such as `black hole binaries <https://en.wikipedia.org/wiki/Binary_black_hole>`_), which will be present in the data recorded by the `LISA gravitational wave observatory <https://en.wikipedia.org/wiki/Laser_Interferometer_Space_Antenna>`_. In this case, working in the time-, rather than frequency-domain, can be particularly helpful for detailed modelling of the waveform, as some realistic and expected time-dependent general relativistic and astrophysical processes will have significant impact on GW production and propagation, and thus the morphology of the waveforms. Depending on the details of the signal one is analysing, one can expect ``dolfen`` to give an accurate representation of realistic likelihood functions for these kinds of LISA signals around 200 to 5000 times faster than the corresponding frequency-domain analysis using all the datapoints.


The code is available at: https://github.com/jethrolinley/dolfen.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   installation
   settings
   usage
   marginalisation

.. toctree::
   :maxdepth: 1
   :caption: Examples:

   linear-model
   bilby-example
