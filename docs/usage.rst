=================
Recommended Usage
=================

There are various simplifying settings that can be chosen which can drastically speed up the likelihood evaluation, at the cost of assured accuracy. This provides a rapid feedback arena for working out some useful properties of the likelihood, for example, giving a good idea of the parameter errors/standard deviations, so that appropriate priors can be chosen (note also that a rough first idea of the priors can be derived from the inverse of the Fisher information matrix). Once the PE setup is arranged, one can then revert to the default ``dolfen`` settings which may take a little longer to evaluate (but still much faster than not using ``dolfen``!), however, the fine details of the returned likelihood function will then be likely to be highly accurate. 

The process of obtaining a high accuracy likelihood/posterior with ``dolfen`` could proceed as follows:

.. note::

    Whilst one may obtain essentially indistinguishable likelihoods/posteriors from different choice of MCS (the number of *maximum correlated samples*), especially when using the FIM preservation method, one can only theoretically *expect* high accuracy from dolfen when MCS is larger than some PSD defined constant (i.e., left to default settings).

#. Choose a low number of samples, around 100, say, and choose zero maximum correlated samples when initialising ``dolfen``; :code:`numsamps=100, MCS_override=0`. Then 100 data points only will be used. Note that without choosing a downsampling method, ``dolfen`` by default first tries to find a solution via approximate FIM preservation, then reverts to single factor noise reduction to find a solution. If your model is very expensive and ``dolfen`` is taking too long to find a FIM preserving solution with the low :code:`numsamps`, you could try again with single factor method selected on initialisation of ``dolfen``, i.e. :code:`prsrv_FIM=False`. 

#. After successful initialisation of ``dolfen``, start computing posteriors (with bilby, for example) to work out and fine tune the appropriate priors and any other PE environment and/or sampler settings.

#. With appropriate PE settings found, we can now produce another likelihood/posterior using default ``dolfen`` settings, to return a result which is likely highly accurate. Unfortunately, the accuracy can never be guaranteed (we would of course have to know the fully sampled, originally defined likelihood, in order to compare them so that we can give this guarantee; however if we had this, we wouldn't need ``dolfen``!). What we can now do however is to produce multiple posteriors to ensure consistency, implying convergence and thus an accurate representation of the original likelihood. 
