.. dolfen documentation master file, created by
   sphinx-quickstart on Sat Dec 16 10:42:56 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


Welcome to dolfen's documentation!
========================================

======
Dolfen
======

``dolfen``: Downsampling Likelihood Function Estimation

``dolfen`` is an algorithm that infers the likelihood function used in Bayesian parameter estimation (PE) by downsampling. It is for use in simulated environments only and is particularly useful for rapid likelihood calculation of long duration, slowly evolving signals consisting of large numbers of datapoints, where the Bayesian likelihood is expensive.

The code is available at: https://github.com/jethrolinley/dolfen.

.. toctree::
      :maxdepth: 0
   :caption: Contents:

   installation
   settings

.. toctree::
      :maxdepth: 0
   :caption: Examples:

   bilby-example
