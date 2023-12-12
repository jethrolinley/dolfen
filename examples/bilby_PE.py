#!/usr/bin/env python

# Example of using dolfen with bilby (needs a separate install)


import bilby
from bilby.core.prior import PriorDict
import numpy as np
import dolfen
import LISAPSD
import TaylorT3hpwf_margphase as TaylorT3hpwf
import nessai




##-----------------------------------
## Signal model injection parameters
##-----------------------------------
t_c = 0.0
chirp_mass = 10.0
mass_ratio = 0.82
spin1 = 0.32
injection_parameters = dict(
    chirp_mass=chirp_mass, mass_ratio=mass_ratio, a_1=spin1, a_2=spin1, tilt_1=0., tilt_2=0.,
    phi_12=0., phi_jl=0., luminosity_distance=20.0, theta_jn=0.68, psi=0.659,
    phase=0.5, geocent_time=t_c, ra=1.375, dec=-0.7108)




##-----------------
## Bayesian priors
##-----------------
priors = PriorDict()
for kparam in injection_parameters:
    exec("priors[kparam] = injection_parameters[\""+kparam+"\"]")

priors["chirp_mass"] = bilby.prior.Uniform(name='chirp_mass', latex_label='$\\mathcal{M}_\\mathrm{c}$', minimum=chirp_mass-0.0001, maximum=chirp_mass+0.0001)
priors["mass_ratio"] = bilby.prior.Uniform(name='mass_ratio', latex_label='$q$', minimum=0.5, maximum=1.0)
priors["a_1"] = bilby.prior.Uniform(name='a_1', latex_label='$\chi_\mathrm{eff}$', minimum=0.3, maximum=0.34)




##------------------------------
## Here are the waveform models
##------------------------------

# wrapper for TaylorT3 for reparameterisationsi (spin aligned!)
def wf_model(times, geocent_time, chirp_mass, mass_ratio, a_1, tilt_1, tilt_2, phi_12, phi_jl, theta_jn,
              luminosity_distance, phase, psi, ra, dec, **kwargs):
    a_2 = a_1
    mass_1, mass_2 = chirpmass_massratio_to_componentmasses(chirp_mass, mass_ratio)
    return TaylorT3hpwf.TaylorT3hpwaveform(times,geocent_time, mass_1,mass_2,a_1, a_2, tilt_1, tilt_2, phi_12,
                                           phi_jl, theta_jn, luminosity_distance, phase, psi, ra, dec)

# wrapper for TaylorT3 for reparameterisations for phase-marginalisation
def wf_model_margphase(times, geocent_time, chirp_mass, mass_ratio, a_1, tilt_1, tilt_2, phi_12, phi_jl, theta_jn,
              luminosity_distance, psi, ra, dec, **kwargs):
    a_2 = a_1
    mass_1, mass_2 = chirpmass_massratio_to_componentmasses(chirp_mass, mass_ratio)
    return TaylorT3hpwf.TaylorT3hpwaveform_margphase(times,geocent_time, mass_1,mass_2,a_1, a_2, tilt_1, tilt_2, phi_12,
                                           phi_jl, theta_jn, luminosity_distance, psi, ra, dec, **kwargs)

def chirpmass_massratio_to_componentmasses(chirpmass,q):
    m1 = chirpmass*(1.+q)**0.2 / q**0.6
    m2 = q*m1
    return m1,m2





##--------------------------------------
## Signal time and duration information
##--------------------------------------

year = 3.15581497632e7
signal_duration = 4.*year
signal_end_to_t_c = 0.1*year
signal_start_time = t_c - signal_duration - signal_end_to_t_c
signal_end_time = signal_start_time + signal_duration

## Here we just want to get the Nyquist frequency of the injected signal
m1, m2 = chirpmass_massratio_to_componentmasses(chirp_mass, mass_ratio)
injection_parameters['mass_1'] = m1
injection_parameters['mass_2'] = m2
LISANyqFreq = TaylorT3hpwf.T3NyquistFrequency(signal_end_time,**injection_parameters)
del injection_parameters['mass_1']
del injection_parameters['mass_2']

## Now define the vector of time stamps at which the signal is sampled
LISANyqsampletimes = np.arange(signal_start_time, signal_end_time, 1./LISANyqFreq)





##----------------------------------------
## Initialise the dolfen/bilby likelihood
##----------------------------------------

## Bilby makes a couple of specific calls to the likelihood function.
## To use dolfen as the likelihood for bilby, you can use a subclass inheriting from bilby.Likelihood
class dolfilby_likelihood(dolfen.likelihood, bilby.Likelihood):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

outdir = "./dolfilby_PE_example"
likelihood = dolfilby_likelihood(LISANyqsampletimes, wf_model, injection_parameters, LISAPSD.PSDfunc, priors=priors,
                                    save_load_FIM=True, overwritesavedFIM=False, numprocs=0, numsamps=120,
                                    GWmodel_margphase=wf_model_margphase, MCS_override=-1,
                                    resume_dir=outdir)

likelihood.parameters = injection_parameters.copy()
print("Likelihood value at injection parameters:",likelihood.log_likelihood())




##-----------------
## Run the sampler
##-----------------

bilby.utils.check_directory_exists_and_if_not_mkdir(outdir)
label = "bilby_example"

result = bilby.run_sampler(
    likelihood=likelihood, priors=likelihood.priors, sampler='nessai', seed=int(10000*np.random.rand()),
    injection_parameters=injection_parameters, resume=True,
    outdir=outdir, label=label,
    flow_class='GWFlowProposal',
    checkpointing=True, plot=False,
    analytic_priors=True)

## Plot the result
result.plot_corner()

