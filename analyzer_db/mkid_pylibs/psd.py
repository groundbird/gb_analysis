import numpy as np
import scipy
from scipy import signal
from . import misc
from .psddata import PsdData
from . import fit

def calc_psd_welch(amp, phs, rate, fold=None, is_norm=False):
    '''
    Calculate power-spectrum-density.
    :param amp: list of amplitude
    :param phs: list of phase
    :param rate: sampling rate
    :param fold: Folded-number for the PSD calculation. Lower number will provide more precise spectrum.
                 If not specified, the optimum number will be used.
    '''
    size = 2**int(np.floor(np.log2(len(amp))))
    if fold is None: fold = 2**(int(len(str(len(amp))))-2)
    myamp=amp[:size]
    myphs=np.unwrap(phs[:size])
    if is_norm:
        myamp = myamp/np.mean(myamp)
    frq_, amp_ = scipy.signal.welch(myamp, fs=rate, nperseg=size/(fold+1))
    frq_, phs_ = scipy.signal.welch(myphs, fs=rate, nperseg=size/(fold+1))
    return frq_,amp_,phs_

def calc_psd_fromdata(toddata, kind=None,fold=None, is_norm=False):
    data,kind = misc.get_data_type(toddata, kind)
    toddata.info['type_psd'] = kind

    ### PSD calc
    #return calc_psd_welch(data.amplitude, data.phase, toddata.rate, fold)
    if is_norm:
        return calc_psd_welch(data.amplitude, data.phase, toddata.rate, fold, is_norm)
    else:
        return calc_psd_welch(data.coramplitude, data.corphase, toddata.rate, fold, is_norm)

def fit_psd(psd, swpfitresult=None, Qr=None, fr=None, fullfit=True, varname='phase', varrange=None, nsigma=5., silent=False):
    '''
    Fit PSD to obtain the quasi-particle lifetime
    To calculate tauqp, fr and Qr (or swpfitresult including them) should be assigned
    You can select one of two fitting functions:
       |--- full fitting: fitting together TLS, GR, and white noises.
       |--- simple fitting: fitting only GR and white noises. No 1/f components in the lower frequency region.
    Fitting will be performed with following steps:
     - Pre-fitting
     - Highly spiked data points (>=`nsigma` [sigma] deviation) will be removed
     - Skimmed data will be fitted again
     - If `redchi2`>1 or any fitted param < 0 in full-fitting, simple fitting will be performed. (The first fitting will be kept as the result if the simple one has no benefits).
    :param psd: PSDData variable or dictionary including 'frq' and varname.
    :param swpfitresult: kidfit.KidFitResult variable. `Qr` and `fr` will be used.
    :param Qr and fr: Can be directly specified for `Qr` and `fr` instead of specifying swpfitresult.
    :param fullfit: True or False. If True(False), full(simple) fitting with(without) TLS==1/f noise curve will be used.
    :param varname: variable name to be fitted. one of ['amp','phs','amplitude','phase']
    :param varrange: fittig range
    :param nsigma: If the data points >= `nsigma` deviation from the fitted function, they will be removed.
    :param silent: If False, detail fitting information will be dumped.
    '''

    ## Fixed parameters from sweep fitting
    if swpfitresult is not None:
        if 'Qr' not in swpfitresult.params.keys() or 'fr' not in swpfitresult.params.keys():
            raise RuntimeError(f'ERROR:: Qr and fr are required for psd fit. Check swpfitresult parameters: {swpfitresult.params.keys()}')
        Qr = swpfitresult.params['Qr'].value
        fr = swpfitresult.params['fr'].value #[Hz]
    elif Qr is None or fr is None:
        raise RuntimeError(f'Invalid input for Qr and fr. Specify swpfitresult OR Qr and fr')

    tres = Qr/(np.pi*fr) #[sec]
    if not silent: print(f'Used fixed parameters: Qr={Qr}, fr={fr}, tres={tres}')

    ## Fitting data
    if type(psd) is PsdData:
        xs = psd.f[1:]
        if varname == 'phase' or varname == 'phs':
            ys_orig = psd.phase[1:]
        elif varname == 'amplitude' or varname == 'amp':
            ys_orig = psd.amplitude[1:]
        else:
            raise RuntimeError("ERROR:: Assign \'amplitude\'(\'amp\') or \'phase\'(\'phs\') for PSD fitting.")
    elif type(psd) is dict:
        xs = np.asarray(psd['frq'])
        ys_orig = np.asarray(psd[varname])

    ys = 10*np.log10(ys_orig)

    if varrange is not None:
        xs = xs[varrange]
        ys = ys[varrange]
        ys_orig = ys_orig[varrange]

    ## Initial values (may need to tune for your KIDs...) and fitting function
    i_WhiteNoise = np.array(ys_orig[-10:]).mean()*0.8
    i_tqp = 2e-6
    if fullfit:
        i_GRconst = 10**((np.median(ys))/10)
        i_TLSexp = 1.
        tmpy = ys_orig[:20].mean()
        tmpx = xs[:20].mean()
        i_TLSconst = (tmpy - i_GRconst/((1+(2*np.pi*tmpx*i_tqp)**2)*(1+(2*np.pi*tmpx*tres)**2)) - i_WhiteNoise)*(tmpx**i_TLSexp)
        ### Fitting functions
        def fitfunc(f,TLSconst=i_TLSconst,TLSexp=i_TLSexp, GRconst=i_GRconst,tqp=i_tqp,WhiteNoise=i_WhiteNoise):
            fys = TLSconst/(f**TLSexp) + GRconst/((1+(2*np.pi*f*tqp)**2)*(1+(2*np.pi*f*tres)**2)) + WhiteNoise
            logfys = 10*np.log10(np.abs(fys))
            return logfys
        if not silent: print(f'Full fitting initial params: WhiteNoise={i_WhiteNoise}, tqp={i_tqp}, GRconst={i_GRconst}, TLSexp={i_TLSexp}, TLSconst={i_TLSconst}')
    else:
        i_GRconst = 10**((np.array(ys[:10]).mean())/10)
        ### Fitting functions
        def fitfunc(f,GRconst=i_GRconst,tqp=i_tqp,WhiteNoise=i_WhiteNoise):
            fys = GRconst/((1+(2*np.pi*f*tqp)**2)*(1+(2*np.pi*f*tres)**2)) + WhiteNoise
            logfys = 10*np.log10(np.abs(fys))
            return logfys
        if not silent: print(f'Simple fitting initial params: WhiteNoise={i_WhiteNoise}, tqp={i_tqp}, GRconst={i_GRconst}')

    ## First fitting
    if not silent: print("FIRST FITTING::")
    preresult = fit.fit(fitfunc,xs,ys, silent=silent)
    result = preresult

    ## Remove too high spikes (>= nsigma [sigma] from fitted function)
    diff = ys - preresult.eval(np.array(xs))
    mean = diff.mean()
    sigm = diff.std()
    xs2 = xs[abs(diff-mean)<nsigma*sigm]
    ys2 = ys[abs(diff-mean)<nsigma*sigm]
    #if not silent: print(f"Found spike (>={nsigma} deviation): x = {xs[i]}, y = {ys[i]}")

    ## If found high spikes, fitting is re-run
    if len(xs2) != len(xs):
        if not silent: print("FITTING WItH SKIMMED DATA::")
        result = fit.fit(fitfunc,np.array(xs2),np.array(ys2), silent=silent)
        if preresult.info['redchi']<result.info['redchi']:
            result = preresult
            if not silent: print(f'The first fitting is better: redchi2 in 1st fit = {preresult.info["redchi"]} <--> redchi2 in skimmed data fit = {result.info["redchi"]}')

    ## subtract TLS component and re-fit
    #ys2 = 10*np.log10(np.abs(10**(np.array(ys2)/10) - result.params['TLSconst'].value/(xs2**(result.params['TLSexp'].value))))
    #def tqpfunc(f,GRconst=result.params['GRconst'].value,tqp=result.params['tqp'].value,WhiteNoise=result.params['WhiteNoise'].value):
    #    fys = GRconst/((1+(2*np.pi*f*tqp)**2)*(1+(2*np.pi*f*tres)**2)) + WhiteNoise
    #    if any(fys<0):
    #        if not silent:
    #            print('WARNING:: psd fit function returns minus value. (unexpected) ')
    #            print(f' -------- GRconst={GRconst:.3e}, tqp={tqp:.3e}, WhiteNoise={WhiteNoise:.3e}')
    #    logfys = 10*np.log10(np.abs(fys))
    #    return logfys
    #tqp_preresult = fit.fit(tqpfunc,np.array(xs2),np.array(ys2), silent=silent)
    # remove too high spikes and re-fit
    #xs3 = []
    #ys3 = []
    #for i,b_ind in enumerate(misc.is_difflarge(ys2,tqp_preresult.eval(np.array(xs2)))):
    #    if b_ind:
    #        xs3.append(xs2[i])
    #        ys3.append(ys2[i])
    #    elif not silent:
    #        print(f"Found spike (>={nsigma} deviation): x = {xs2[i]}, y = {ys2[i]}")
    #tqp_result = tqp_preresult
    #if len(xs3) != len(xs2):
    #    tqp_result = fit.fit(tqpfunc,np.array(xs3),np.array(ys3), silent=silent)
    #    if tqp_preresult.info['redchi']<tqp_result.info['redchi']: tqp_result = tqp_preresult

    ## if full-fit is not good (reducedChi2>1 or fitparam<0), try simple-fit
    minusval = any(np.logical_and([prm.value < 0 for prm in result.params.values()], result.params.keys() != 'tqp'))
    if fullfit and (result.info['redchi']>1 or minusval):
        if not silent: print(f"WARNING:: full-fit looks not good. try simple fit.")
        ### pick up only high-frequency region (>1kHz)
        xs3 = xs2[xs2>1e3]
        ys3 = ys2[xs2>1e3]

        ### initial values and fitting function
        i_GRconst = 10**((np.array(ys3[:10]).mean())/10)
        def fitfunc(f,GRconst=i_GRconst,tqp=i_tqp,WhiteNoise=i_WhiteNoise):
            fys = GRconst/((1+(2*np.pi*f*tqp)**2)*(1+(2*np.pi*f*tres)**2)) + WhiteNoise
            logfys = 10*np.log10(np.abs(fys))
            return logfys
        if not silent: print(f'Simple fitting initial params: WhiteNoise={i_WhiteNoise}, tqp={i_tqp}, GRconst={i_GRconst}')

        ### Fitting
        if not silent: print("SIMPLE FITTING::")
        simpleresult = fit.fit(fitfunc,np.array(xs3),np.array(ys3), silent=silent)

        ### Check fitting quality
        minusval = any(np.logical_and([prm.value < 0 for prm in result.params.values()], result.params.keys() != 'tqp'))
        if simpleresult.info['redchi']<1 and not minusval:
            result = simpleresult
            print(f"fit-psd():: Full-fit looks not good. Adopt simple fit.")
        elif not silent: print(f'The first fitting is better: redchi2 in 1st fit = {result.info["redchi"]} <--> redchi2 in simple fit = {simpleresult.info["redchi"]}')

    return result
