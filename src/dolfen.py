# -*- coding: utf-8 -*-

import copy
import inspect
import os
import string
import multiprocessing
import random
import scipy
import scipy.signal
from scipy.special import logsumexp
import hashlib
import numpy as np
import scipy.integrate
import pickle
import warnings
from scipy.linalg import LinAlgWarning


class likelihood:
    """
    The 'downsampled' approximation of the Bayesian likelihood function defined by the signal and noise models, and injection parameters.

    Parameters
    ----------
    t : 1D numpy array
        The set of time-stamps of the 'original' dataset, sampled at least
        at the Nyquist frequency of all waveforms with non-negligible probability
        (just above the Nyquist frequency of the injected signal should do).
    model : callable
        The user-defined waveform model, the first argument must be the times-stamps, t.
    inj_params : dict
        The signal's injection parameters for each parameter of the
        signal model: key=param name, value=param value.
    PSD : callable
        The 'Power Spectral Density', a function of frequency, f. This encodes the
        noise model. Must return 1D array when passed a 1D array.
    priors : dict
        The parameters (dict keys) and the priors (values). If the parameter is known
        the prior should be a float, otherwise it should be an object with a ".minimum"
        and a ".maximum" property.
    numsamps : int
        The number of datapoints to use for downsampling.
    scheme : string
        The 'downsampling scheme' that should be employed by dolfen. Currently accepted
        schemes are 
            "rand": random downsampling;

            "unif": uniform downsampling (a.k.a. decimation);

            "hybrid": a 50/50 split of uniform and random downsampling;

            "cluster,n": random samples are taken from 'n' uniformly spaced blocks of
            size 'blocksize'. if n is not specified, n is set to the number of model parameters;

            "prand": 'pseudo-random' downsampling. Mainly used for testing dolfen, for
            any given system, a seed is fixed for random sample selection.
    prsrv_FIM : Bool
        True to use the (approximate) FIM preservation method of downsampling. This is the most
        accurate method. False to use the minimisation of Jeffreys' divergence of the FIMs of fully-
        and down-sampled datasets. If not supplied, dolfen will try to find FIM preservation
        solutions a number of times, then revert to FIM Jeffreys' divergence minimisation if no
        solution can be found.
    MCS_override : int
        The number of maximum correlated samples of any given sample as determined by the inverse
        autocorrelation function. Any negative value means dolfen will automatically compute this.
    addnoise : Bool
        True to add a random noise realisation to the data from the PSD.
    f_low_cut, f_high_cut : floats
        Lowest and highest frequencies at which the PSD is irrelevant to the likelihood; the PSD
        below f_low_cut is set equal to the PSD at f_low_cut, the PSD above f_high_cut is set equal
        to the PSD at f_high_cut.
    GWmodel_margphase : callable
        Dolfen was created for studying likelihoods of gravitational waves of low-mass CBCs in LISA.
        Numerically marginalising the phase of the time-domain CBC GW waveform is possible as per
        the example included with dolfen; the model is required to separately return the plus and
        cross GW polarisations so that the phase marginalised likelihood can be computed quickly.
    blocksize : int
        For use with "cluster" downsampling; size of blocks from which to sample
    param_diffs : list
        The step-sizes for numerical derivatives of signal model with respect to its parameters,
        in case dolfen has difficulty in automatically computing them.
    save_load_FIM : Bool
        True if one wishes to write to disk the computed Fisher matrix and derivative step-sizes.
    overwritesavedFIM : Bool
        If FIM data has already been saved to disk, overwrite it.
    phase_int_points : int
        The number of discrete integration points to use for marginalising GW phase.
    resume_dir : string
        Directory for writing/reading downsampling solution information for given dolfen inputs.
    numprocs : int
        Number of processors to use to compute the Fisher information matrix. Uses all
        available if not set.
    forcerun : Bool
        True to force run even if the MCS is large and downsampling would not be expected to yield
        a significant reduction in likelihood evaluation time.
    """

    def makeSubSet(self, numsamps, numblocks=0, blocksize=1, scheme="rand", resume_dir='', prsrv_FIM=True):
        if numsamps < 1:
            print("Number of datapoints entered less than 1!!")
            exit()
        if not resume_dir == '':
            if not os.path.exists(resume_dir):
                os.mkdir(resume_dir)
                print("Created resume directory:",resume_dir)
            elif os.path.isfile(resume_dir.rstrip('/')+'/'+self.hashfilename+'.pickle'):
                with open(resume_dir.rstrip('/')+'/'+self.hashfilename+'.pickle', 'rb') as f:
                    resume_data = pickle.load(f)
                    f.close()
                print("Resuming downsampling subset data from:",resume_dir.rstrip('/')+'/'+self.hashfilename+'.pickle')
                self.required_samps =        resume_data['srs']
                self.selsamps_idx_reqsamps = resume_data['sir']
                self.subsampletimes =        resume_data['subsampletimes']
                self.subdata =               resume_data['subdata']
                self.sample_weights =        resume_data['sample_weights']
                self.factor =                resume_data['factor']
                self.Fconst =                resume_data['Fconst']
                self.leftpad = [0. for p in self.required_samps if p < 0]
                self.rightpad = [0. for p in self.required_samps if p >= self.N]
                self.invACFfunc *= self.factor
                self.max2sidedInvACF *= self.factor
                self.invACFfuncsqrt *= np.sqrt(self.factor)
                self.max2sidedInvACFsqrt *= np.sqrt(self.factor)
                return True

        self.subN = numsamps
        print("\nDatapoint subset proposal",self.nt+1,"of",self.maxtries,"... ",end="")

        # choose the samples wanted in the subset (uniformly spaced chunks here)
        self.selected_samps = np.zeros(self.subN)
        srs = list()
        sir = list()
        if scheme == "unif":
            if blocksize == 0:
                blocksize = 1
            dbbs = int((self.N-blocksize)/numsamps) # "distance between block starts"
            for i in range(numsamps-1):
                self.selected_samps[i*blocksize:(i+1)*blocksize]=np.arange(i*dbbs,i*dbbs+blocksize,dtype='int')
                if len(srs)>0:
                    add = np.arange(max(i*dbbs-self.max_corr_samps,srs[-1]+1),i*dbbs+blocksize+self.max_corr_samps+1,dtype='int')
                else:
                    add = np.arange(i*dbbs-self.max_corr_samps,i*dbbs+blocksize+self.max_corr_samps+1,dtype='int')
            
            self.selected_samps[-int(blocksize):] = np.arange(self.N-blocksize,self.N,dtype='int')
        elif scheme[:7] == "cluster":
            try:
                numclusters = int(scheme.split(",")[1])
            except:
                numclusters = self.FN
            if blocksize*numclusters < numsamps and blocksize > 0:
                print("Please ensure that blocksize*numofclusters > samps")
                exit()
            if blocksize < 2: # i.e., not passed parameter! Or badly chosen
                blocksize = int(np.floor(self.N/numclusters))
            startpoints = np.linspace(0,self.N-blocksize,numclusters, dtype=int).tolist()
            sampspercluster = numsamps//numclusters*np.ones(numclusters, dtype=int)
            sampspercluster[:numsamps % numclusters] += 1
            random.shuffle(list(sampspercluster))
            self.selected_samps = []
            for i in range(numclusters):
                self.selected_samps += list(np.asarray(sorted(list(set(random.sample(range(blocksize),sampspercluster[i]))))) + startpoints[i])
        elif scheme == "prand" or scheme == "rand" or scheme == "hybrid":
            if scheme == "prand":
                random.seed(self.nt)
            if blocksize > 1:
                print("Warning: random/hybrid downsampling with blocksize > 1 not yet supported, changing blocksize to 1 and numblocks =",self.subN)
                blocksize=1
                numblocks=self.subN
            hybrid_even_samps = []
            if scheme == "hybrid":
                hybrid_even_samps = np.linspace(0, self.N, numsamps//2, endpoint=False, dtype=int).tolist()
                hybrid_even_samps.pop(0)
            self.selected_samps = sorted(list(set([0,self.N-1]+random.sample(range(1,self.N-1),self.subN-2-len(hybrid_even_samps))+hybrid_even_samps))) # select first and last sample then random from the middle
            self.subN = len(self.selected_samps)
        else:
            print("Please choose a downsampling scheme! From: unif, rand, hybrid, cluster.")
            exit()
    
        for i in self.selected_samps:
            if len(srs)>0:
                add = np.arange(max(i-self.max_corr_samps,srs[-1]+1),i+self.max_corr_samps+1,dtype='int')
            else:
                add = np.arange(i-self.max_corr_samps,i+self.max_corr_samps+1,dtype='int')
            srs += list(add)
            sir.append(len(srs)-self.max_corr_samps-1)
 
        self.required_samps = srs
        self.selsamps_idx_reqsamps = sir
        self.subsampletimes = self.t[[ti for ti in srs if ti >= 0 and ti < self.N]]

        self.leftpad = [0. for p in self.required_samps if p < 0]
        self.rightpad = [0. for p in self.required_samps if p >= self.N]

        self.subFisher = np.array(self.getFisherInfoSub(),dtype='float64')

        np.set_printoptions( linewidth=186)

        singular_params = [self.FIMparams.index(spar) for spar in self.singular_params]
        singular_params += list(np.where(self.subFisher.diagonal() <= 0.)[0])
        singular_params = sorted(set(singular_params))
        for i in range(len(singular_params)):
            j = singular_params[len(singular_params) - i - 1]
            self.subFisher = np.delete(self.subFisher,j,0)
            self.subFisher = np.delete(self.subFisher,j,1)
            self.fullFisher = np.delete(self.fullFisher,j,0)
            self.fullFisher = np.delete(self.fullFisher,j,1)


        # Now compute the noise reduction factor, or the time sample weight
        wf, Pf = scipy.linalg.eig(self.fullFisher)
        Pf_inv = scipy.linalg.inv(Pf)
        subFishDiag_f = np.linalg.multi_dot((Pf_inv,self.subFisher,Pf))
        fullFishDiag_f = np.linalg.multi_dot((Pf_inv,self.fullFisher,Pf))

        ws, Ps = scipy.linalg.eig(self.subFisher)
        try:
            Ps_inv = scipy.linalg.inv(Ps)
        except:
            return False
        subFishDiag_s = np.linalg.multi_dot((Ps_inv,self.subFisher,Ps))
        fullFishDiag_s = np.linalg.multi_dot((Ps_inv,self.fullFisher,Ps))

        significance_ratio = 1.e-10
        significant_eigvals = np.where(wf > np.max(wf)*significance_ratio, 1., 0.) # An extra check, essentially PCA, to remove v. low info (diagonalised) params, to prevent numerical breakdown.

        if np.sum(significant_eigvals) < 1.:
            return False

        self.factor = 0.0
        try:
            filtered_fullFishDiag = [i for indx,i in enumerate(np.diag(fullFishDiag_f)) if significant_eigvals[indx] == 1.]
            filtered_subFishDiag  = [i for indx,i in enumerate(np.diag(subFishDiag_s)) if significant_eigvals[indx] == 1.]
            self.factor_det = (np.prod(filtered_fullFishDiag)/np.prod(filtered_subFishDiag))**(1./np.sum(significant_eigvals))

            numerator_diag = significant_eigvals*(fullFishDiag_s.diagonal()/subFishDiag_s.diagonal())
            denominator_diag = significant_eigvals*(subFishDiag_f.diagonal()/fullFishDiag_f.diagonal())
            if any(t < 0 for t in numerator_diag):
                return False
            if any(t < 0 for t in denominator_diag):
                return False
            self.factor_jeffreys = np.sum(numerator_diag)/np.sum(denominator_diag)
            if self.factor_jeffreys <= 0.0:
                return False
            self.factor_jeffreys = np.sqrt(self.factor_jeffreys)
            self.factor = self.factor_jeffreys
        except:# Exception as excpt:
            return False

        if not (np.isfinite(self.factor) and self.factor > 0.0):
            return False
        if not prsrv_FIM:# or np.sum(significant_eigvals) < 2:
            print("Noise reduction factor computed:",self.factor)
            self.sample_weights = 1.#np.ones(self.subN)
            self.invACFfunc *= self.factor
            self.max2sidedInvACF *= self.factor
            self.invACFfuncsqrt *= np.sqrt(self.factor)
            self.max2sidedInvACFsqrt *= np.sqrt(self.factor)
            self.subFisher *= self.factor
        else:
            diag_singular_params = list(np.where(subFishDiag_f.diagonal() <= 0.)[0])
            fullFishColMaxOffDiagonal = np.max(fullFishDiag_f-np.diag(fullFishDiag_f.diagonal()), axis=0)
            diag_singular_params += list(np.where(fullFishDiag_f.diagonal() <= fullFishColMaxOffDiagonal)[0])
            diag_singular_params = sorted(set(diag_singular_params))
            for i in range(len(diag_singular_params)):
                j = diag_singular_params[i]#len(diag_singular_params) - i - 1]
                subFishDiag_f[j,:] = 0.   # = np.delete(subFishDiag_f,j,0)
                subFishDiag_f[:,j] = 0.   # = np.delete(subFishDiag_f,j,1)
                fullFishDiag_f[j,:] = 0.  # = np.delete(fullFishDiag_f,j,0)
                fullFishDiag_f[:,j] = 0.  # = np.delete(fullFishDiag_f,j,1)

            if not self.getSampleWeightsApprox(fullFishDiag_f, Pf_inv):#significant_Fisher):
                print("Get sample weights failed ("+str(self.subN),"samples)")
                return False

            self.factor = 1.

        print("success!\n\nTotal number of samples required to be computed:",len(srs))

        self.subdata = np.asarray(self.whitenVec(self.model(self.subsampletimes,*self.orig_ordered_inj_params)))

        if self.addnoise:
            noisetimeseries = np.random.randn(self.subN)/(2.*np.sqrt(self.dt))
        else:
            noisetimeseries = np.zeros(self.subN)

        self.Fconst = np.sum(self.subdata*self.sample_weights*self.subdata) - np.log(2.*np.pi/self.margphaseintpoints)

        self.parameters = self.inj_params
        if self.log_likelihood() > 0.0:
            print("\nWarning: the maximum log likelihood value is greater than zero. This is benign and can happen as a result of the downsampling process, especially when the SNR of the injected signal is very low. Is the signal injected correctly?\n")

        self.subdata = self.subdata + noisetimeseries


        if not resume_dir == '':
            resume_data = dict()
            resume_data['srs'] = self.required_samps
            resume_data['sir'] = self.selsamps_idx_reqsamps
            resume_data['subsampletimes'] = self.subsampletimes
            resume_data['factor'] = self.factor
            resume_data['Fconst'] = self.Fconst
            resume_data['subdata'] = self.subdata
            resume_data['sample_weights'] = self.sample_weights
            with open(resume_dir.rstrip('/')+'/'+self.hashfilename+'.pickle', 'wb') as f:
               pickle.dump(resume_data, f, pickle.HIGHEST_PROTOCOL)
               f.close()

        print("Finished initialising likelihood for data subset                   ")
        return True

    def getSampleWeightsApprox(self, Fisher, Pf_inv):
        FIMlen = len(Fisher)
        whitened_sub_dh_dx_i = np.zeros((FIMlen, self.subN))
        for i in range(FIMlen):
            whitened_sub_dh_dx_i[i] = self.whitenVec(self.sub_dh_dx_i[i])

        sing_FIMparams = np.where(Fisher.diagonal() <= 0.0)[0]
        sigFIMlen = FIMlen - len(sing_FIMparams)

        # whitened AND DIAGONALISED derivatives!! (do the linear transf. with Pf_inv)
        whitened_sub_dh_dx_i = np.matmul(Pf_inv, whitened_sub_dh_dx_i)
        nonsing_Fish = Fisher.copy()
        for i in range(len(sing_FIMparams)):
            j = sing_FIMparams[len(sing_FIMparams) - i - 1]
            whitened_sub_dh_dx_i = np.delete(whitened_sub_dh_dx_i,j,0)
            nonsing_Fish = np.delete(nonsing_Fish,j,0)
            nonsing_Fish = np.delete(nonsing_Fish,j,1)
            

        coeff_mat = np.zeros((sigFIMlen,sigFIMlen))#, dtype='float128')
        timelist = np.linspace(10.0, 11.0, self.subN)
        for i in range(sigFIMlen):
            for j in range(sigFIMlen):
                if nonsing_Fish[i][i] > 0.0 or j == i:
                    coeff_mat[i][j] = np.dot((timelist**j)*whitened_sub_dh_dx_i[i], whitened_sub_dh_dx_i[i])

        try:
            warnings.filterwarnings(action='ignore', category=LinAlgWarning)#, module='sklearn')
            time_coeffs = scipy.linalg.solve(coeff_mat, np.diag(nonsing_Fish))
        except Exception as e:
            print(e)
            return False
        diffs = np.abs(np.ones(sigFIMlen) - np.dot(coeff_mat,time_coeffs)/np.diag(nonsing_Fish))
        if np.max(diffs) > 0.0001:
            print("Failed accuracy check.")
            return False
        else:
            print("Passed accuracy check.")
        self.sample_weights = np.zeros(self.subN)
        for i in range(sigFIMlen):
            self.sample_weights += time_coeffs[i]*(timelist**i)

        ######################    Now do some checks on the solution!!
        if np.min(self.sample_weights) <0.:
            print("Some sample weights are negative")
            return False

        subFish_00 = np.zeros((sigFIMlen,sigFIMlen))
        for i in range(sigFIMlen):
            for j in range(sigFIMlen):
                subFish_00[i][j] = np.dot(self.sample_weights*whitened_sub_dh_dx_i[i],whitened_sub_dh_dx_i[j])
            
        if np.min(subFish_00.diagonal()) < 0.:
            return False
        normaliseFish = subFish_00.copy()
        for i in range(sigFIMlen):
            normaliseFish[i,:] /= np.sqrt(subFish_00[i][i])
            normaliseFish[:,i] /= np.sqrt(subFish_00[i][i])

        if np.max(np.abs(normaliseFish)) > 1.00001:# or np.sum(self.sample_weights) < 0.:
            print("Bad subsample selection for FIM preservation!")
            return False
        return True


    def whitenVec(self, vec):
        vec = np.concatenate((self.leftpad, vec, self.rightpad))
        decorr_vec = [np.dot(self.max2sidedInvACFsqrt, vec[i-self.max_corr_samps:i+self.max_corr_samps+1]) for i in self.selsamps_idx_reqsamps]
        return decorr_vec

    def subFIMInnerProd(self,vec1,vec2):
        subvec1 = [np.dot(self.max2sidedInvACFsqrt, vec1[i-self.max_corr_samps:i+self.max_corr_samps+1]) for i in self.selsamps_idx_reqsamps]
        subvec2 = [np.dot(self.max2sidedInvACFsqrt, vec2[i-self.max_corr_samps:i+self.max_corr_samps+1]) for i in self.selsamps_idx_reqsamps]
        return np.dot(subvec1,subvec2)

    def getFisherInfoSub(self):
        print("\nComputing Fisher matrix for data subset....")
        subFN = self.FN
        Fisher = np.zeros((subFN,subFN),dtype='float128')
        self.sub_dh_dx_i = np.zeros((subFN,len(self.required_samps)),dtype='float128')
        for i in range(subFN):
            self.sub_dh_dx_i[i] = np.concatenate(([0. for p in self.required_samps if p < 0],self.getDeriv(self.t[[p for p in self.required_samps if p >= 0 and p < self.N]],i),[0. for p in self.required_samps if p >= self.N]))
        for i in range(subFN):
            for j in range(subFN):
                if i >= j:
                    Fisher[i][j] = self.subFIMInnerProd(self.sub_dh_dx_i[i],self.sub_dh_dx_i[j])
                    Fisher[j][i] = Fisher[i][j]
        return Fisher


    def getInvACFfunc(self, t, PSDfunc):
        # this is in order to split up (shorten) the vector so that numerical breakdown doesn't occur in the FFTs
        self.samplefreqs = np.linspace(0.,1./self.dt,self.max_corr_samps,endpoint=False)[:self.max_corr_samps//2+self.max_corr_samps%2]
        print("Calculating PSD....                                               \r",end="")
        PSD = np.array(PSDfunc(self.samplefreqs),dtype="float128")
        print("Calculating \"Inverse\" noise Autocorrelation Function....  (FFT size",len(PSD),")                  \r",end="")
        invACFfunc=np.fft.irfft(1./PSD)[:len(PSD)]
        invACFfuncsqrt=np.fft.irfft(1./np.sqrt(PSD))[:len(PSD)]
        return np.array(invACFfunc[:len(invACFfunc)//2],dtype='float128'), np.array(invACFfuncsqrt[:len(invACFfunc)//2],dtype='float128')
    
    def resetsampletimesanddata(self, t, PSD):
        self.N = len(t)
        self.t = t
        self.dt = np.abs(t[1]-t[0])

        corr_cut_off_ok = False
        while not corr_cut_off_ok:
            self.invACFfunc, self.invACFfuncsqrt  = self.getInvACFfunc(t, PSD)
            invACFsqrt_int = np.cumsum(np.abs(self.invACFfuncsqrt))
            if self.MCS_override >= 0:
                corrcutoffsamps = self.MCS_override
                corr_cut_off_ok = True
            try:
                corrcutoffsamps = next(x[0] for x in enumerate(invACFsqrt_int) if x[1] > 0.97*invACFsqrt_int[-1])
            except:
                corrcutoffsamps = len(invACFfuncsqrt)
            if corrcutoffsamps < self.max_corr_samps/10.: # if 97% of invACFsqrt_int is NOT within first 10% of "test max_corr_samps", then probably NOT useful to downsample
                corr_cut_off_ok = True
            elif self.max_corr_samps == self.N:
                if self.forcerun:
                    corr_cut_off_ok = True
                    print("Continuing using probably bad ACF function! Since (forcerun=True)")
                else:
                    print("Cannot find MCS < orig num of samps! (I.e., need to compute all datapoints and downsampling is of no benefit) Aborting! If you want to go ahead anyway, rerun using \"forcerun=True\"")
                    return False
            else:
                self.max_corr_samps += 100
                if self.max_corr_samps >= self.N:
                    self.max_corr_samps = self.N

        self.corrcutoff = corrcutoffsamps*self.dt # set the max correlation length to when correlations drop below about 1/10000 of the max correlation
        self.max_corr_samps=min(corrcutoffsamps,len(t))
        if self.MCS_override >= 0:
            print("Maximum correlation time lag computed as ~",round(self.corrcutoff,2),"seconds. However, the MCS_override argument has been set, therefore...:               ")
            self.max_corr_samps = self.MCS_override
        else:
            print("Maximum correlation time lag computed and set to ~",round(self.corrcutoff,2),"seconds.                                ")
        print("    - maximum number of correlated neighbouring samples:", self.max_corr_samps)

        self.invACFfunc = self.invACFfunc[:self.max_corr_samps+1]
        self.invACFfuncsqrt = self.invACFfuncsqrt[:self.max_corr_samps+1]
        self.max2sidedInvACF = np.concatenate((np.flip(self.invACFfunc),self.invACFfunc[1:]))
        self.max2sidedInvACFsqrt = np.concatenate((np.flip(self.invACFfuncsqrt),self.invACFfuncsqrt[1:]))
        self.lenSIACF = len(self.max2sidedInvACF)

        print("Calculating Fisher matrices (full dataset)....                        \r",end="")
        return self.getFisherInfo(self.t)


    def getFIMchunk(self,times, ids, q):
        Fisher = np.zeros((self.FN,self.FN),dtype='float128')
        N = len(times)
        blocksize = min(100000,N)
        remainingsamples = N
        blockaddsize = 0
        blockremovesize = 0
        blockdh_dx_i = np.zeros((self.FN,blocksize))
        while remainingsamples > 0:
            if remainingsamples <= blocksize:
                blockaddsize = remainingsamples
                blockremovesize = 0
                blockdh_dx_i = np.zeros((self.FN,blockaddsize),dtype='float128')
            else:
                blockaddsize = blocksize
                blockremovesize = self.max_corr_samps

            blocktimes = times[N-remainingsamples:N-remainingsamples+blockaddsize]
            for i in range(self.FN):
                blockdh_dx_i[i] = self.getDeriv(blocktimes,i)
            blocksig = self.model(blocktimes,*self.orig_ordered_inj_params)

            for i in range(self.FN):
                for j in range(self.FN):
                    if i >= j:
                        Fisher[i][j] += self.approxInvAcfInnerProd(blockdh_dx_i[i],blockdh_dx_i[j])
                        if blockremovesize > 0:
                            Fisher[i][j] -= self.approxInvAcfInnerProd(blockdh_dx_i[i][-blockremovesize:],blockdh_dx_i[j][-blockremovesize:])
            remainingsamples -= blockaddsize - blockremovesize
        print("  "+u'\u2713',end="")
        q.put(Fisher)

    def opt_stepsize(self, model, boundsfactor, diffparam, times, params, numtests=250):
        lentimes = len(times)
        derivs = np.zeros((numtests,lentimes))
        stepsizes=10.**np.linspace(-20.0,-4.0,numtests)*boundsfactor
        ## check between 1e-20*boundsfactor and 1e-4*boundsfactor
        failed_attempts = 0
        for i in range(numtests):
            try:
                leftparams = params.copy()
                leftparams[diffparam] -= stepsizes[i]
                leftsig = model(times,**leftparams)
                rightparams = params.copy()
                rightparams[diffparam] += stepsizes[i]
                rightsig = model(times,**rightparams)
                derivs[i] = (rightsig-leftsig)/(2.*stepsizes[i])
                if np.sum(np.abs(derivs[i])) == 0.:
                    derivs[i] = 10.*i*np.ones(lentimes)
                    failed_attempts += 1
            except:
                derivs[i] = 10.*i*np.ones(lentimes)
                failed_attempts += 1
        if failed_attempts >= numtests:
            return False
        derivdiffs = np.diff(derivs,axis=0)
        flatness = np.sum(derivdiffs**2,axis=1)

        flatness_smoothed = np.zeros(len(flatness))
        for i in range(-2,3):
            flatness_smoothed += np.roll(flatness,i)
        flatness = np.ma.masked_equal(flatness_smoothed, 0.)
        flatness = np.ma.masked_equal(flatness,np.min(flatness))
        flatness = flatness[np.logical_not(np.isnan(flatness))]
        flatness = flatness[2:-2]
        return stepsizes[2+np.where(flatness == flatness.min())[0]][0]



    def getFisherInfo(self, times):
        self.FN = self.numparams
        self.FIMparams = list()
        knownparams=""
        for i in range(self.numparams):
            if not isinstance(self.origpriors[str(self.ordered_inj_keys[i])], (int, float, complex)):
                self.FIMparams.append(self.ordered_inj_keys[i])
            else:
                self.FN -= 1
                knownparams += ", "+str(self.ordered_inj_keys[i])
        knownparams = knownparams[2:]
        self.initFIMparams = self.FIMparams.copy()

        PSDstr = str(self.PSD(self.samplefreqs[1:])/self.samplefreqs[1:])
        PSDnum = hashlib.md5()
        PSDnum.update(PSDstr.encode('utf-8'))
        PSDencode = PSDnum.hexdigest()

        # this should contain all of the unique signal/PSD/FIM param info for any given system...
        signal_det_info = knownparams + "st" + str(times[0]) + "t0" + str(times[1]-times[0]) + "PSD" + str(PSDencode) + str(self.orig_ordered_inj_params)
        hash_md5 = hashlib.md5()
        hash_md5.update(signal_det_info.encode('utf-8'))
        hashfilename = hash_md5.hexdigest()
        self.hashfilename = hashfilename

        if not knownparams == "":
            print("Parameter(s)",knownparams,"are known and will not contribute to FIM                       ")
        Fisher = np.zeros((self.FN,self.FN),dtype='float128')

        derivchecklen = 500

        # only need 2-sided/central diff, since need to calculate 2 signals anyway! (very large signal not passed to likelihood anymore)
        self.derivtype=-np.ones(self.FN) # deriv calc type:-1=not set,0=centraldiff,1=symbolic/auto


        fileloc = str(os.path.dirname(os.path.abspath(__file__))) + "/dsFIMinf/"
        dthetasloaded = False
        if self.sv_ld_FIM and not self.overwriteFIM:
            if not os.path.exists(fileloc):
                os.mkdir(fileloc)
            elif os.path.isfile(fileloc+hashfilename+'.dt.npy'):
                self.dthetas = np.load(fileloc+hashfilename+'.dt.npy')
                self.derivtype = np.zeros(self.FN)
                print("Stepsizes already evaluated in, and loaded from \"./dsFIMinf/"+"".join(str.split(hashfilename))+".dt.npy\"")
                dthetasloaded = True
        if not dthetasloaded:
            if self.sv_ld_FIM:
                print("\nEvaluating stepsizes...")
            self.dthetas=np.empty(self.FN)
            for i in range(self.FN):
                modelparamnum=self.ordered_inj_keys.index(self.FIMparams[i])
                derivativesuccess = False
                if self.param_diffs_supplied:
                    self.derivtype[i] = 0
                    self.dthetas[i] = self.param_diffs[self.FIMparams[i]]
                    derivativesuccess=True
                else:#if self.diff_evol_opt_step_size:
                    self.found_opt_stepsize = -1.
                    print("Optimising step-size for numerical derivative of model w.r.t param",str(self.FIMparams[i]),"....", end="")
                    derivchecktimes = times[::len(times)//derivchecklen]
                    prior = self.origpriors[str(self.FIMparams[i])]
                    boundsfactor = np.abs(prior.minimum-prior.maximum)
                    res = self.opt_stepsize(self.model, boundsfactor, self.FIMparams[i], derivchecktimes, self.inj_params)
                    if res:
                        self.derivtype[i] = 0
                        derivativesuccess = True
                        print("done")
                        self.dthetas[i]=res#.x[0]
                    else:
                        print("failed! Removing parameter from FIM.")
                        self.derivtype[i] = -1

            if self.sv_ld_FIM:
                np.save(fileloc+hashfilename+'.dt', self.dthetas)

        # Check if FIM has already been computed and saved for this system - if so, load it!
        if self.sv_ld_FIM and not self.overwriteFIM:
            if not os.path.exists(fileloc):
                os.mkdir(fileloc)
            elif os.path.isfile(fileloc+hashfilename+'.npy'):
                self.fullFisher = np.load(fileloc+hashfilename+'.npy')
                print("FIM already evaluated in, and loaded from \"./dsFIMinf/"+"".join(str.split(hashfilename))+".npy\"")
                return True

        Fisher = np.zeros((self.FN,self.FN),dtype='float128')
       
        print("\nCalculating Fisher Information Matrix:                                                       ")
        timesets = list()
        timesetoverlaps = list()
        numprocs = self.numprocs
        sampsperset = self.N//numprocs
        while sampsperset <= self.max_corr_samps*10 and numprocs > 1:
            numprocs -= 1
            sampsperset = self.N//numprocs
        for i in range(numprocs-1):
            timesets.append(times[i*sampsperset:(i+1)*sampsperset+self.max_corr_samps])
            timesetoverlaps.append(times[(i+1)*sampsperset:(i+1)*sampsperset+self.max_corr_samps])
        timesets.append(times[(numprocs-1)*sampsperset:])

        # remove overlap regions
        print("Preparing tidy-up of segment overlap regions...                         ")
        blockdh_dx_i = np.zeros((self.FN,self.max_corr_samps))
        for k in range(len(timesetoverlaps)):
            for i in range(self.FN):
                blockdh_dx_i[i] = self.getDeriv(timesetoverlaps[k],i)
            for i in range(self.FN):
                for j in range(self.FN):
                    if i >= j:
                        Fisher[i][j] -= self.approxInvAcfInnerProd(blockdh_dx_i[i],blockdh_dx_i[j])

        # add all regions with multiprocessing
        print("Continuing with bulk FIM calculation (using",numprocs,"threads)...                                      ")
        
        q = multiprocessing.Queue()
        processes = []
        for i in range(numprocs):
            p = multiprocessing.Process(target=self.getFIMchunk, args=(timesets[i],i,q,))
            processes.append(p)
        
        for x in processes:
            x.start()

        print("\nChunk(s)                         ",end="")
        for i in range(numprocs):
            print(" "*(3-len(str(i+1)))+str(i+1),end="")
        print("\nCompleted the following chunk(s):",end="")

        for x in processes:
            Fisher += q.get()

        for x in processes:
            x.join()

        print("\033[F\033[F"," "*400,end="\r")
        print(" "*400,end="\r")
        print("\033[F\033[F\033[F\033[F\033[F ... done computing Fisher Information Matrix                                             ")

        for i in range(self.FN):
            for j in range(self.FN):
                 if i < j:
                    Fisher[i][j]=Fisher[j][i]
        self.fullFisher = np.array(Fisher,dtype='float64')
        if self.sv_ld_FIM:
            np.save(fileloc+hashfilename+'.npy', self.fullFisher)

        return True
        
    def approxInvAcfInnerProd(self,sig1,sig2):
        invACFmat_x_sig2 = scipy.signal.convolve(sig2,self.max2sidedInvACF,mode="same")
        return np.dot(sig1,invACFmat_x_sig2)
    
    def flatPSD_osb(self,f):        
        flatten_len = len(np.where(f < self.f_low_cut))
        f[:flatten_len] = np.ones(flatten_len)*self.f_low_cut
        flatten_len = len(np.where(f > self.f_high_cut))
        if flatten_len > 0:
            f[-flatten_len:] = np.ones(flatten_len)*self.f_high_cut
        return self.origPSD(f)


    def __init__(self,
                 t,
                 model,
                 inj_params,
                 PSD,
                 priors,
                 numsamps=400,
                 scheme="rand",
                 prsrv_FIM=None,
                 MCS_override=-1,
                 addnoise=False,
                 f_low_cut=0.0,
                 f_high_cut=0.0,
                 GWmodel_margphase=None,
                 blocksize=1,
                 param_diffs=None,
                 save_load_FIM=True,
                 overwritesavedFIM=False,
                 phase_int_points=2357,
                 resume_dir='',
                 numprocs=0,
                 forcerun=False,
    ):

        print("\nInitialising dolfen....                     \r",end="")
        super(likelihood, self).__init__()
        self.model = model
        self.inj_params = inj_params
        self.origPSD = PSD
        self.scheme = scheme
        self.addnoise = addnoise
        self.sv_ld_FIM = save_load_FIM
        self.param_diffs = param_diffs
        self.f_low_cut = f_low_cut
        self.f_high_cut = f_high_cut
        self.overwriteFIM = overwritesavedFIM
        self.GWmodel_margphase = GWmodel_margphase
        self.margphaseintpoints = phase_int_points
        self.forcerun = forcerun # force run even if invACF func has inadequate padding
        self.prsrv_FIM = copy.deepcopy(prsrv_FIM)
        self.MCS_override = MCS_override
        self.origpriors = copy.deepcopy(priors)
        self.max_corr_samps = 1000 # first guess, will be incremented if no good solution (approximation) to invACFfunc is found

        if numprocs == 0:
            self.numprocs = max(1,(multiprocessing.cpu_count()-1)//2)
        else:
            self.numprocs = numprocs

        if f_low_cut > 0.0 or f_high_cut > 0.0:
            self.PSD = self.flatPSD_osb # outside_sig_band
        else:
            self.PSD = self.origPSD

        # These lines of code infer the parameters from the provided function
        parameters = inspect.getfullargspec(model).args
        parameters.pop(0)
        self.numparams = len(parameters) # number of param space/FIM dimensions
        self.parameters = dict.fromkeys(parameters)

        self.orig_ordered_inj_params = list()
        for i in range(self.numparams):
            self.orig_ordered_inj_params = self.orig_ordered_inj_params + list([inj_params.get(str(parameters[i]))])

        self.GWmodel_hasmargphase = False
        if GWmodel_margphase != None:
            self.parameters.pop('phase')
            parameters.remove('phase')
            self.numparams -= 1
            self.GWmodel_hasmargphase = True

        self.function_keys = self.parameters.keys()
        self.ordered_inj_params = list()
        self.ordered_inj_keys = list()
        for i in range(self.numparams):
            self.ordered_inj_params = self.ordered_inj_params + list([inj_params.get(str(parameters[i]))])
            self.ordered_inj_keys = self.ordered_inj_keys + list([str(parameters[i])])

        print("\nOrdered model parameters/keys:",self.ordered_inj_keys)
        if param_diffs == None:
            self.param_diffs_supplied = False
        else:
            self.param_diffs_supplied = True
        diagFIMresult = -1
        print("Original number of samples (before downsampling):",len(t))

        if not self.resetsampletimesanddata(t, self.PSD):
            self.fullFisher = None
            print("\nCould not initialise!\n")
            return

        self.coslist = np.cos(np.linspace(0., 2.*np.pi, self.margphaseintpoints, endpoint=False)) # use endpoint=True if using np.trapz though!
        self.sinlist = np.sin(np.linspace(0., 2.*np.pi, self.margphaseintpoints, endpoint=False))
        self.cos2list = self.coslist*self.coslist
        self.sin2list = self.sinlist*self.sinlist
        self.cossinlist = self.sinlist*self.coslist

        singular_params = np.sort(np.where(np.sum(np.abs(self.fullFisher), axis=1)==0.)[0])
        self.singular_params = []
        for i in range(len(singular_params)):
            j = singular_params[len(singular_params) - i - 1]
            print("Zero information on parameter",self.FIMparams[j],"!! Removing from FIM (makes FIM singular or poorly representative)")
            self.singular_params.append(self.FIMparams[j])
            self.fullFisher = np.delete(self.fullFisher,j,0)
            self.fullFisher = np.delete(self.fullFisher,j,1)
        self.initFIM = self.fullFisher.copy()
        print("FIM parameters:", self.initFIMparams)

        if self.GWmodel_hasmargphase:
            self.log_likelihood = self.log_likelihood_marg
        else:
            self.log_likelihood = self.log_likelihood_nomarg

        self.nt = 0 # nt = current Number of Tries
        self.maxtries = 200
        print("\nTrying to find downsampling solution...")#, end=" ")
        xtrasamps_even = 0
        while self.nt < self.maxtries:
            if self.prsrv_FIM == None: ##### if no downsampling method was selected
                if self.nt < self.maxtries//2:
                    prsrv_FIM = True
                else:
                    prsrv_FIM = False
        
            if self.scheme == "unif":
                xtrasamps_even = int((-1)**(self.nt-1)*np.ceil(self.nt/2))
            if not self.makeSubSet(numsamps=numsamps+xtrasamps_even,blocksize=blocksize, scheme=self.scheme, resume_dir=resume_dir,prsrv_FIM=prsrv_FIM):
                self.nt += 1
            else:
                break

        if self.nt >= self.maxtries:
            print("\nUh-oh!! Critical Error: Unable to find random sampling subset passing minimal accuracy criteria",
              "or producing non-singular Fisher information matrix after",self.nt,"tries. Aborting.")
            print("Could not find downsampling solution. Solutions are more easily found using:\n   1. more samples,\n   2. random sampling,\n   3. single factor rather than FIM preservation (less accurate)")
            return
        elif self.scheme == "unif":
            numsamps = numsamps + int((-1)**(self.nt-1)*np.ceil(self.nt/2))
            print("Found solution using",numsamps,"uniformly spaced samples.")

        for i in range(len(singular_params)):
            j = singular_params[len(singular_params) - i - 1]
            self.FIMparams = np.delete(self.FIMparams,j)
        self.FN -= len(singular_params)

        self.priors = copy.deepcopy(priors)
 
        if self.prsrv_FIM == None:
            if prsrv_FIM == True:
                print("\nSuccessful initialisation using approximate FIM preservation!\n")
            else:
                print("\nCould not find FIM preservation solution, successfully initialised using single factor noise reduction instead!\n")
        else:
            print("\nSuccessful initialisation!\n")

    def getDeriv(self, times, i):
        if self.derivtype[i] == 0:
            k = self.ordered_inj_keys.index(self.FIMparams[i])
            diff_inj_params = list(self.orig_ordered_inj_params).copy()
            diff_inj_params[k] += self.dthetas[i]
            cdiff = self.model(times,*diff_inj_params)
            diff_inj_params[k]-=2*self.dthetas[i]
            cdiff -= self.model(times,*diff_inj_params)
            cdiff /= (2.*self.dthetas[i])
            return cdiff
        elif self.derivtype[i] == 1:
            return self.getFIMAGDeriv(times,i)
        else:
            return np.zeros(len(times))
    
    def log_likelihood_marg(self):
        model_parameters = {k: self.parameters[k] for k in self.function_keys}# if not k=='phase'}
        h_cos, h_sin = self.GWmodel_margphase(self.subsampletimes,**model_parameters)

        white_cos = np.array(self.whitenVec(h_cos))#*self.sample_weights
        white_sin = np.array(self.whitenVec(h_sin))#*self.sample_weights

        A = np.sum(white_cos*self.sample_weights*white_cos)
        B = -2.*np.sum(white_cos*self.sample_weights*white_sin)
        C = -2.*np.sum(white_cos*self.sample_weights*self.subdata)
        D = np.sum(white_sin*self.sample_weights*white_sin)
        E = 2.*np.sum(white_sin*self.sample_weights*self.subdata)

        exponent = -self.dt*(A*self.cos2list + B*self.cossinlist + C*self.coslist + D*self.sin2list + E*self.sinlist + self.Fconst)
        return logsumexp(exponent)


    def log_likelihood_nomarg(self):
        model_parameters = {k: self.parameters[k] for k in self.function_keys}
        diff = self.subdata - self.whitenVec(self.model(self.subsampletimes,**model_parameters))
        return -np.sum(diff*diff)*self.dt

