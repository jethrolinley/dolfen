===============
Suggested Usage
===============

There are various simplifying settings that can be chosen which can drastically speed up the likelihood evaluation, at the cost of assured accuracy. This provides a rapid feedback arena for working out some useful properties of the likelihood, for example, giving a good idea of the parameter errors/standard deviations, so that appropriate priors can be chosen (note also that a rough first idea of the priors can be derived from the Fisher information matrix). Once the PE setup is arranged, one can then revert to the default ``dolfen`` settings which may take a while longer to evaluate (but still likely very much faster than not using ``dolfen``!), however, the fine details of the returned likelihood function will then be likely to be highly accurate. 

The process of obtaining a high accuracy likelihood/posterior with ``dolfen`` may go as follows:

.. note::

    One may obtain mostly indistinguishable likelihoods/posteriors from different choice of ``MCS``, the number of maximum correlated samples, especially when using the FIM preservation method. However, whilst lower ``MCS`` means faster likelihood evaluation, unfortunately as a result of the theoretical background, one can only expect high accuracy from dolfen when ``MCS`` is high (or left to default settings).

#. Choose a low number of samples, around 100, say, and choose zero maximum correlated samples when initialising ``dolfen``; :code:`numsamps=100, MCS=0`. Then 100 data points only will be used. Note that without choosing a downsampling method, ``dolfen`` by default first tries to find a solution via approximate FIM preservation, then reverts to single factor noise reduction to find a solution. If your model is very expensive and ``dolfen`` is taking too long to find a FIM preserving solution with the low :code:`numsamps`, you could try again with single factor method selected on initialisation of ``dolfen``, i.e. :code:`prsrv_FIM=False`. 

#. When ``dolfen`` is initialised, start computing posteriors (with bilby, for example) and work out the appropriate priors and any other PE settings.

#. With appropriate PE settings found, we can produce another likelihood/posterior using default ``dolfen`` settings to return a result which is likely highly accurate. Unfortunately, the accuracy can never be guaranteed (we would have to know the fully sampled, originally defined likelihood to compare them in order to give this guarantee, but if we had this, we wouldn't need ``dolfen``!). So, we can produce multiple posteriors to ensure consistency, implying convergence and thus an accurate representation of the original likelihood. 
