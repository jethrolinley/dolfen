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

``dolfen`` is an algorithm for Bayesian analysis that infers the likelihood function by downsampling. It is designed in particular for simulated environments only and is particularly useful for rapid likelihood calculation of, long duration, slowly evolving signals consisting of very many datapoints, where the Bayesian likelihood is expensive.

The code is available at: https://github.com/jethrolinley/dolfen.

.. toctree::
      :maxdepth: 1
   :caption: Contents:

   installation

.. toctree::
      :maxdepth: 1
   :caption: Examples:

   bilby-example
