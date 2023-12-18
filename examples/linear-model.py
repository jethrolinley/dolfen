#!/usr/bin/env python

# Example of using dolfen with simple linear model

import numpy as np
import dolfen

def model(t,d):
    return t*d

inj_param = {'d':0.0}

def PSD(f):
    return np.ones(len(f))

class prior():
    def __init__(self, minimum, maximum):
        self.minimum = minimum
        self.maximum = maximum

priors = {'d' : prior(-0.3, 0.3)}

t = np.linspace(0., 1., 2000)

llhood = dolfen.likelihood(t, model, inj_param, PSD, priors)
llhood.parameters = {'d': 0.1}
llhood.log_likelihood()

