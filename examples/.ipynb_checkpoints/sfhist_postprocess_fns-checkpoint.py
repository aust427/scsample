import numpy as np
import matplotlib.pyplot as plt
import dense_basis as db
import hickle
from scipy.interpolate import PchipInterpolator, interp1d
from scipy.signal import savgol_filter
from george import kernels
import george
from scipy.optimize import minimize

from astropy.cosmology import FlatLambdaCDM

cosmo = db.cosmo

convert_to_microjansky = db.convert_to_microjansky
make_filvalkit_simple = db.make_filvalkit_simple
calc_fnu_sed_fast = db.calc_fnu_sed_fast
convert_to_splined_spec = db.convert_to_splined_spec

def tofloatarr(line):
    temp = line.split()
    return np.array([float(x) for x in temp])

def ingest_sfhist_file(work_dir, sfhfile):

    id1 = []
    id2 = []
    galz = []
    sfharrs = []
    issfh = False

    file_path = work_dir + sfhfile

    ctr = 0
    try:
        with open(file_path, 'r') as file:
            for line in file:

                if ctr == 0:
                    temp = line.split()
                    H_0 = float(temp[0])
                    Omega_m = float(temp[1])

                if ctr == 1:
                    ntbins = int(line)
                    print(ntbins)
                    # Zbins = tofloatarr(line)

                if ctr == 2:
                    tbins = tofloatarr(line)

                if ctr > 2:
                    temp = tofloatarr(line)
                    if len(temp) == 3:
                        id1.append(temp[0])
                        id2.append(temp[1])
                        galz.append(temp[2])
                        if issfh == True:
                            sfharrs.append(np.array(sfharr))
                        sfharr = []
                        issfh = True
                    else:
                        sfharr.append(temp)

                ctr = ctr+1
            sfharrs.append(np.array(sfharr))
    except FileNotFoundError:
        print(f"File '{file_path}' not found.")
    except IOError:
        print(f"Error reading file '{file_path}'.")
        
    return id1, id2, galz, sfharrs, tbins



def calctimes(
    timeax, sfh, nparams, 
    sfr_Myr = 10, 
    scale = 'log', log_stretch = -2
):
    """
    Function to convert an input time-series into a tuple.
    In case of an SFH, this is (M*, <SFR>|_Myr, tx/tuniv)
    
    Inputs:
        timeax: (np array) times
        sfh: (np array) SFR(t) at those times
        nparams: (int) number of tx params 
        sfr_Myr: (optional, float [10]) time in Myr over which to average recent SFR
        scale: (optional, 'lin' or ['log']) scale for dividing the tx params
        log_stretch: (optional [-2]) how finely to divide the last few percentiles. only for 'log' scale
        
    Outputs:
        mass: (float) 
        sfr: (float) SFR averaged over sfr_Myr
        tx: (array of floats) tx parameters normalized to 1
    
    """
    

    sfr_index = np.argmin(np.abs(timeax-sfr_Myr/1e3))
    mass = np.log10(np.trapz(sfh,timeax*1e9))
    sfr = np.log10(np.nanmedian(sfh[-sfr_index:]))
    
    massint = np.cumsum(sfh)
    massint_normed = massint/np.amax(massint)
    tx = np.zeros((nparams,))
    
    if sfr != -np.inf:
        mass_formed_sfr_Myr = sfr_Myr * 1e6 * (10**sfr)
        mass_frac_rem = 1 - (mass_formed_sfr_Myr / (10**mass))
    else:
        mass_frac_rem = 1.0
    
    if scale == 'log':
        txints = (1-np.flip(np.logspace(log_stretch,0,nparams+2),0))/(1-10**log_stretch)*mass_frac_rem
        txints = txints[1:-1]
    elif scale == 'lin':
        txints = np.linspace(0,mass_frac_rem,nparams+2)[1:-1]
        
    # calculating the tXs from the mass quantiles
    for i in range(nparams):
        tx[i] = timeax[np.argmin(np.abs(massint_normed - txints[i]))]
    return mass, sfr, tx/np.amax(timeax)

def calctuple_sfhZt(sfhZt, tvals, nparam = 3, sfr_Myr = 100, scale='log', log_stretch = -2):
    """
    Function to calculate sfh tuple at each metallicity value. 
    Expected input is a log SFR numpy array with dimensions Ntime x Nmet
    """

    nmet = sfhZt.shape[1]
    sfhZtuple = np.zeros((nparam+3, nmet))
    for i in range(nmet):
    
        m, sfr, tx = calctimes(tvals, 10**sfhZt[0:,i], nparam, sfr_Myr = sfr_Myr, scale=scale, log_stretch=log_stretch)
        Z_tuple = db.calctimes_to_tuple([m, sfr, tx])
        sfhZtuple[0:,i] = Z_tuple
        
    return sfhZtuple


def nll(p, gp, y):
    gp.set_parameter_vector(p)
    ll = gp.log_likelihood(y, quiet=True)
    return -ll if np.isfinite(ll) else 1e25

def grad_nll(p, gp, y):
    gp.set_parameter_vector(p)
    return -gp.grad_log_likelihood(y, quiet=True)

def linear_interpolator(x,y,res = 1000):

    interpolator = interp1d(x,y)
    x_pred = np.linspace(np.amin(x), np.amax(x), res)
    y_pred = interpolator(x_pred)

    return x_pred, y_pred

def Pchip_interpolator(x,y,res = 1000):

    interpolator = PchipInterpolator(x,y)
    x_pred = np.linspace(np.amin(x), np.amax(x), res)
    y_pred = interpolator(x_pred)

    return x_pred, y_pred

def gp_george_interpolator(x,y,res = 1000, Nparam = 3, decouple_sfr = False):

    yerr = np.zeros_like(y)
    yerr[2:(2+Nparam)] = 0.001/np.sqrt(Nparam)
    if len(yerr) > 26:
        yerr[2:(2+Nparam)] = 0.1/np.sqrt(Nparam)
    if decouple_sfr == True:
        yerr[(2+Nparam):] = 0.1

#     kernel = np.var(y) * (kernels.Matern32Kernel(np.median(y)) *( kernels.LinearKernel(np.median(y), order=2)))
#     kernel = np.var(y) * (kernels.Matern32Kernel(np.median(y)) *( kernels.DotProductKernel(np.median(y))))
    kernel = np.var(y) * (kernels.Matern32Kernel(np.median(y)) *( kernels.PolynomialKernel(np.median(y), order=2)))
    gp = george.GP(kernel, solver=george.HODLRSolver)

    gp.compute(x.ravel(), yerr.ravel())
    
#     try:
    # optimize kernel parameters
    p0 = gp.get_parameter_vector()
    results = minimize(nll, p0, jac=grad_nll, method="L-BFGS-B", args = (gp, y))
    gp.set_parameter_vector(results.x)
#     print(results)
#     except:
#         print('couldn\'t optimize GP parameters.')

    x_pred = np.linspace(np.amin(x), np.amax(x), res)
    y_pred, pred_var = gp.predict(y.ravel(), x_pred, return_var=True)

    return x_pred, y_pred

def tuple_to_sfh_new(sfh_tuple, zval, 
                     interpolator = 'gp-george', 
                     scale='log', log_stretch = -2, 
                     cosmo = cosmo, 
                     sfr_Myr = 100, 
                     Nsfr_const_perc = 3, 
                     zerosfr_at_tbb = True, 
                     smooth_savgol = 9,
                     vb = False):
    
    Nparam = sfh_tuple.shape[0]-3
    mass = 10**sfh_tuple[0]
    sfr = 10**sfh_tuple[1]
    
    if sfr != -np.inf:
        mass_formed_sfr_Myr = sfr_Myr * 1e6 * (10**sfr)
        mass_frac_rem = 1 - (mass_formed_sfr_Myr / (10**mass))
    else:
        mass_frac_rem = 1.0
    
    if scale == 'log':
        txints = (1-np.flip(np.logspace(log_stretch,0,Nparam+2),0))/(1-10**log_stretch)*mass_frac_rem
        txints = txints[0:-1]
    elif scale == 'lin':
        txints = np.linspace(0,mass_frac_rem,nparams+2)[0:-1]
    
    mass_quantiles = txints
    time_quantiles = np.zeros_like(mass_quantiles)
    time_quantiles[1:] = sfh_tuple[3:]
    
    mass_quantiles = np.append(mass_quantiles,[1.0])
    time_quantiles = np.append(time_quantiles,[1.0])
    
    # ================== SFR at t=0 ========================
    
    if zerosfr_at_tbb == True:
        # SFR smoothly increasing from 0 at the big bang
        mass_quantiles = np.insert(mass_quantiles,1,[0.0000])
        time_quantiles = np.insert(time_quantiles,1,[1e-3])
    
    # ================ sSFR at t = t_obs ====================
    
    #correct for mass loss? (small at small sfr_Myr so ignore for now)
    mass_formed_sfr_Myr = sfr_Myr * 1e6 * (10**sfh_tuple[1])
    mass_remaining = (10**sfh_tuple[0]) - mass_formed_sfr_Myr
    
    mass_frac = 1 - mass_formed_sfr_Myr / (10**sfh_tuple[0])
    sfr_tuniv = 1 - (sfr_Myr / (cosmo.age(zval).value*1e3))
    
    mass_quantiles = np.insert(mass_quantiles, -1, [mass_frac])
    time_quantiles = np.insert(time_quantiles, -1, [sfr_tuniv])
    
    # ========= removing duplicate values ==============
    
    tqmask = (np.diff(time_quantiles)>0)
    tqmask = np.insert(tqmask,0,[True])
    time_quantiles = time_quantiles[tqmask]
    mass_quantiles = mass_quantiles[tqmask]
   
    if not np.all(np.isfinite(mass_quantiles)):
        print('Error in interpolation...')
        return np.nan, np.nan
    
    # ========= splining the time-cumulative SFR array ==============
    
    if interpolator == 'linear':
        time_arr_interp, mass_arr_interp = linear_interpolator(time_quantiles, mass_quantiles)
    elif interpolator == 'pchip':
        try:
            time_arr_interp, mass_arr_interp = Pchip_interpolator(time_quantiles, mass_quantiles)
        except:
            print((np.diff(time_quantiles)<=0))
            print(time_quantiles, np.diff(time_quantiles))
            time_arr_interp, mass_arr_interp = Pchip_interpolator(time_quantiles, mass_quantiles)
    elif interpolator == 'gp-george':
        
        time_arr_interp, mass_arr_interp = gp_george_interpolator(time_quantiles, mass_quantiles)
        time_arr_interp2, mass_arr_interp2 = Pchip_interpolator(time_quantiles, mass_quantiles)
        
        # dealing with -ve SFR in a bin
        for k in range(1,len(mass_quantiles)-1):
            trange = (time_arr_interp >= time_quantiles[k]) & (time_arr_interp <= time_quantiles[k+1])
            if np.nanmedian(np.diff(mass_arr_interp[trange])) < 1e-4:
                mass_arr_interp[trange] = mass_arr_interp2[trange]
    else:
        raise Exception('specified interpolator does not exist: {}. \n use one of the following: gp_george, gp_sklearn, linear, and pchip '.format(interpolator))
    
    # ========= derivative to compute SFR array ==============
    
    dx = np.mean(np.diff(time_arr_interp))
    sfh = np.gradient(mass_arr_interp, dx, edge_order=2)
    if smooth_savgol > 0:
        sfh = savgol_filter(sfh, smooth_savgol, 1)
    sfh[sfh<0] = 0
    timeax = time_arr_interp * cosmo.age(zval).value
    
    
    # ======== scale SFH to match input mass & SFR =============
    
    massformed = np.trapz(x=timeax*1e9, y=sfh)
    sfh = sfh / massformed * (mass)
    temp = np.argmin(np.abs(timeax-sfr_Myr/1e3))
    sfh[-temp:] = sfr

    if vb == True:

        plt.figure(figsize=(7,7))
        plt.plot([0,1],[0,1],'k--',alpha=0.3)
        plt.plot(time_quantiles, mass_quantiles,'o')
        plt.plot(time_arr_interp, mass_arr_interp)
        plt.xlabel('time / t$_{univ}$'); plt.ylabel('fraction of mass formed')
        plt.show()

        plt.plot(timeax, sfh)
        plt.xlabel('time'); plt.ylabel('SFR')
        plt.show()
    
    return sfh, timeax

def compress_sfh(galid, 
                 nparam = 30, 
                 sfhist_vars = [],
                 sfr_Myr = 10, dustval = 0.2, 
                 interpolator = 'gp-george'):
    
    galz, sfharrs, tbins, Zbins = sfhist_vars
    zval = galz[galid]
    z_index = np.argmin(np.abs(tbins - cosmo.age(zval).value)) + 1
    sfhZt, tvals, Zvals = np.log10(sfharrs)[galid].copy()[0:z_index,0:], tbins.copy()[0:z_index], Zbins
    sfhZtuple = calctuple_sfhZt(sfhZt, tvals, sfr_Myr=sfr_Myr, nparam = nparam,scale='log', log_stretch=-2)
    return sfhZtuple




def makespec(specdetails, priors, sp = db.mocksp, cosmo = cosmo, filter_list = [], filt_dir = [], 
             return_spec = False, peraa = False, input_sfh = False, interpolator = 'gp-george', addlines=True):
    """
    makespec() works in two ways to create spectra or SEDs from an input list of physical paramters. 
    with input_sfh = False, specdetails = [sfh_tuple, dust, met, zval]
    with input_sfh = True, specdetails = [sfh, timeax, dust, met, zval]
    
    return_spec can be true, false, or an array of wavelengths. in case 
    
    it uses the db.mocksp object for the underlying SPS calculations, so make sure it's set accordingly.
    it also uses the priors object for things like decouple_sfr. 
    
    """
    
#     tuple_to_sfh = db.tuple_to_sfh
    convert_to_microjansky = db.convert_to_microjansky
    make_filvalkit_simple = db.make_filvalkit_simple
    calc_fnu_sed_fast = db.calc_fnu_sed_fast
    convert_to_splined_spec = db.convert_to_splined_spec
        
    # hardcoded parameters - offload this to a seprarate function
    sp.params['sfh'] = 3
    sp.params['cloudy_dust'] = True
    sp.params['gas_logu'] = -2
    sp.params['add_igm_absorption'] = True
    sp.params['add_neb_emission'] = addlines
    sp.params['add_neb_continuum'] = addlines
    sp.params['imf_type'] = 1 # Chabrier
    
    # variable parameters
    if input_sfh == True:
        [sfh, timeax, dust, met, zval] = specdetails
#         print('inputsfh max timeax',np.amax(timeax))
    else:
        [sfh_tuple, dust, met, zval] = specdetails
        sfh, timeax = tuple_to_sfh_new(sfh_tuple, zval, interpolator = interpolator)
#         print('tuple max timeax',np.amax(timeax))
    sp.set_tabular_sfh(timeax, sfh)
    # sp.params['dust_type'] = 2
    # sp.params['dust1'] = dust1_rand
    sp.params['dust2'] = dust
    sp.params['logzsol'] = met
    sp.params['gas_logz'] = met # matching stellar to gas-phase metallicity
    sp.params['zred'] = zval
    
    #lam, spec = sp.get_spectrum(tage = cosmo.age(zval).value, peraa = peraa)
    # adding 0.1 Myr here to get the last couple of FSPS SSPs
    lam, spec = sp.get_spectrum(tage = cosmo.age(zval).value+1e-4, peraa = peraa)
    spec_ujy = convert_to_microjansky(spec, zval, cosmo)
    
    if type(return_spec) == type(True):
        
        if return_spec == True:
            return lam, spec_ujy
    
        elif return_spec == False:
            filcurves, _, _ = make_filvalkit_simple(lam, zval, fkit_name = filter_list, filt_dir = filt_dir)
            sed = calc_fnu_sed_fast(spec_ujy, filcurves)
            return sed
    
    elif len(return_spec) > 10:
        return convert_to_splined_spec(spec, lam, return_spec, zval)
    
    else:
        raise('Unknown argument for return_spec. Use True or False, or pass a wavelength grid.')
        
    return 0

def makespec_from_compressed(galid, compressed_sfhs = [],
                             sfhist_vars = [],
                             sfr_Myr = 100, dustval = 0.2, 
                             interpolator = 'gp-george', vb = False):
    
    spec_all = []
    galz, sfharrs, tbins, Zbins = sfhist_vars
    Zvals = Zbins.copy()
    zval = galz[galid]
    for i in range(compressed_sfhs[galid].shape[1]):
        
        if compressed_sfhs[galid][0,i] > -np.inf:
            sfh_tuple = compressed_sfhs[galid][0:,i]
            dust = dustval
            met = Zvals[i]
            tempsfh, temptime = tuple_to_sfh_new(sfh_tuple, zval, interpolator = interpolator, 
                                                     sfr_Myr=sfr_Myr, vb=False)
            specdetails = [tempsfh, temptime, dust, met, zval]
            lam, spec = makespec(specdetails, priors = db.Priors(), sp = db.mocksp, cosmo = cosmo, 
                                    input_sfh = True, return_spec = True, interpolator=interpolator)
            spec_all.append(spec)
            if vb == True:
                print('Sfh for Z=',met)
        else:
            spec_all.append(np.zeros((5994,)))
            
            if vb == True:
                print('no sfh for Z=',met)
    return lam, np.array(spec_all)

def makespec_from_SFHZt(galid,
                         sfhist_vars = [],
                         sfr_Myr = 100, dustval = 0.2, addlines=True,
                         interpolator = 'gp-george', vb = False):
    
    galz, sfharrs, tbins, Zbins = sfhist_vars
    zval = galz[galid]
    z_index = np.argmin(np.abs(tbins - cosmo.age(zval).value)) + 1
    sfhZt, tvals, Zvals = (sfharrs)[galid].copy()[0:z_index,0:], tbins.copy()[0:z_index], Zbins
    spec_all = []

    
    for i in range(sfharrs[galid].shape[1]):
        if sum(sfhZt[0:,i]) > 0:
            dust = dustval
            met = Zvals[i]
            specdetails = [sfhZt[0:,i], tvals, dust, met, zval]
            lam, spec = makespec(specdetails, priors = db.Priors(), sp = db.mocksp, cosmo = cosmo, 
                                    input_sfh = True, return_spec = True, interpolator=interpolator, addlines=addlines)
            spec_all.append(spec)
            if vb == True:
                print('Sfh for Z=',met)
        else:
            spec_all.append(np.zeros((5994,)))
            if vb == True:
                print('no sfh for Z=',met)
    return lam, np.array(spec_all)


def compress_sfh_v2(galid, 
                 nparam = 30, 
                 sfhist_vars = [],
                 sfr_Myr = 10, dustval = 0.2, 
                 scale='log', log_stretch = -2,
                 interpolator = 'gp-george'):
    
    galz, sfharrs, tbins, Zbins = sfhist_vars
    zval = galz[galid]    
    z_index = np.argmin(np.abs(tbins - cosmo.age(zval).value)) + 1
    sfhZt, tvals, Zvals = (sfharrs)[galid].copy()[0:z_index,0:], tbins.copy()[0:z_index], Zbins

    sfht = np.log10(np.sum(sfhZt,1))
    Zt = np.sum(sfhZt.T*Zvals.reshape(-1,1),0) / 10**sfht
    metval = np.nanmedian(Zt)
    Zt[np.isnan(Zt)] = -2.1
    m, sfr, tx = calctimes(tvals, 10**sfht, nparam, sfr_Myr = sfr_Myr, scale=scale, log_stretch=log_stretch)
#     print(m)
    sfhtuple = db.calctimes_to_tuple([m, sfr, tx])
    m, sfr, tx = calctimes(tvals, 10**Zt, nparam, sfr_Myr = sfr_Myr, scale=scale, log_stretch=log_stretch)
    Ztuple = db.calctimes_to_tuple([m, sfr, tx])
#     sfhZtuple = calctuple_sfhZt(sfhZt, tvals, sfr_Myr=sfr_Myr, nparam = nparam,scale='log', log_stretch=-2)
#     return sfhtuple, Zt, tvals, Ztuple
    return sfhtuple, Ztuple, metval



def makespec_v2(specdetails, priors, sp = db.mocksp, cosmo = cosmo, filter_list = [], filt_dir = [], 
             return_spec = False, peraa = False, input_sfh = False, interpolator = 'gp-george', nebem = True, methist = False):
    """
    makespec() works in two ways to create spectra or SEDs from an input list of physical paramters. 
    with input_sfh = False, specdetails = [sfh_tuple, dust, met, zval]
    with input_sfh = True, specdetails = [sfh, timeax, dust, met, zval]
    
    return_spec can be true, false, or an array of wavelengths. in case 
    
    it uses the db.mocksp object for the underlying SPS calculations, so make sure it's set accordingly.
    it also uses the priors object for things like decouple_sfr. 
    
    """
    
#     tuple_to_sfh = db.tuple_to_sfh
    convert_to_microjansky = db.convert_to_microjansky
    make_filvalkit_simple = db.make_filvalkit_simple
    calc_fnu_sed_fast = db.calc_fnu_sed_fast
    convert_to_splined_spec = db.convert_to_splined_spec
    
    
    # hardcoded parameters - offload this to a seprarate function
    sp.params['sfh'] = 3
    sp.params['cloudy_dust'] = True
    sp.params['gas_logu'] = -2
    sp.params['add_igm_absorption'] = True
    sp.params['add_neb_emission'] = nebem
    sp.params['add_neb_continuum'] = nebem
    sp.params['imf_type'] = 1 # Chabrier
    
    # variable parameters
    if input_sfh == True:
        [sfh, timeax, dust, met, zval] = specdetails
#         print('inputsfh max timeax',np.amax(timeax))
    else:
        [sfh_tuple, dust, met, zval] = specdetails
        sfh, timeax = tuple_to_sfh_new(sfh_tuple, zval, interpolator = interpolator)
#         print('tuple max timeax',np.amax(timeax))
    if methist == True:
        sp.set_tabular_sfh(timeax, sfh, Z=met)
    else:
        sp.set_tabular_sfh(timeax, sfh)
        sp.params['logzsol'] = met
        sp.params['gas_logz'] = met # matching stellar to gas-phase metallicity

    # sp.params['dust_type'] = 2
    # sp.params['dust1'] = dust1_rand
    sp.params['dust2'] = dust
    sp.params['zred'] = zval
    
    #lam, spec = sp.get_spectrum(tage = cosmo.age(zval).value, peraa = peraa)
    # adding 0.1 Myr here to get the last couple of FSPS SSPs
    lam, spec = sp.get_spectrum(tage = cosmo.age(zval).value+1e-4, peraa = peraa)
    spec_ujy = convert_to_microjansky(spec, zval, cosmo)
    
    if type(return_spec) == type(True):
        
        if return_spec == True:
            return lam, spec_ujy
    
        elif return_spec == False:
            filcurves, _, _ = make_filvalkit_simple(lam, zval, fkit_name = filter_list, filt_dir = filt_dir)
            sed = calc_fnu_sed_fast(spec_ujy, filcurves)
            return sed
    
    elif len(return_spec) > 10:
        return convert_to_splined_spec(spec, lam, return_spec, zval)
    
    else:
        raise('Unknown argument for return_spec. Use True or False, or pass a wavelength grid.')
        
    return 0


def makespec_from_compressed_v2(galid, compressed_sfhs = [],
                             sfhist_vars = [],
                             sfr_Myr = 100, dustval = 0.2, addlines=True,
                             interpolator = 'gp-george', vb = False):
    
    spec_all = []
    galz, sfharrs, tbins, Zbins = sfhist_vars
    Zvals = Zbins.copy()
    zval = galz[galid]
    
    temp = compressed_sfhs[galid]
    ntemp = temp.shape[0]
    ntemp2 = int((ntemp-1)/2)
    sfhtuple = temp[0:ntemp2]
    Ztuple = temp[ntemp2:2*ntemp2]
    metval = temp[-1]
#     print(sfhtuple, Ztuple, metval)
    sfhtuple, Ztuple, metval
    
    if sfhtuple[0] > -np.inf:
        sfh_tuple = sfhtuple
        dust = dustval
#         met = metval
        tempsfh, temptime = tuple_to_sfh_new(sfh_tuple, zval, interpolator = interpolator, 
                                                     sfr_Myr=sfr_Myr, vb=False)
        tempZt, _ = tuple_to_sfh_new(Ztuple, zval, interpolator = interpolator, 
                                                     sfr_Myr=sfr_Myr, vb=False)
        tempZt = np.log10(tempZt)
        tempZt[tempZt<-2.5] = -2.5
        
        # plt.plot(temptime, tempsfh);plt.show()
        # plt.plot(temptime, 0.019 * 10**tempZt); plt.show()
        
        specdetails = [tempsfh, temptime, dust, 0.0142 * 10**tempZt, zval]
        db.mocksp._zcontinuous = 3
        lam, spec = makespec_v2(specdetails, priors = db.Priors(), sp = db.mocksp, cosmo = cosmo, 
                                input_sfh = True, return_spec = True, interpolator=interpolator, nebem=False, methist=True)
        
        if addlines == True:
            db.mocksp._zcontinuous = 1
            specdetails = [tempsfh, temptime, dust, metval, zval]
            lam, spec2 = makespec_v2(specdetails, priors = db.Priors(), sp = db.mocksp, cosmo = cosmo, 
                                    input_sfh = True, return_spec = True, interpolator=interpolator, nebem=True, methist=False)
            lam, spec3 = makespec_v2(specdetails, priors = db.Priors(), sp = db.mocksp, cosmo = cosmo, 
                                    input_sfh = True, return_spec = True, interpolator=interpolator, nebem=False, methist=False)
            spec = spec + (spec2-spec3)

        return lam, spec
    else:
        print('no SFH data')
        return

def make_seds(spec, lam, galz, filter_list, filt_dir):

    seds = []

    for galid in range(len(spec)):
        zval = galz[galid]
        filcurves, _, _ = db.make_filvalkit_simple(lam, zval, fkit_name = filter_list, filt_dir = filt_dir)
        galseds = []
        for j in range(spec[0].shape[0]):
            galseds.append(db.calc_fnu_sed_fast(spec[galid][j,0:], filcurves))
        seds.append(np.array(galseds))
    return seds

def make_seds_v2(spec, lam, galz, filter_list, filt_dir):

    seds = []

    for galid in range(len(spec)):
        zval = galz[galid]
        filcurves, _, _ = db.make_filvalkit_simple(lam, zval, fkit_name = filter_list, filt_dir = filt_dir)
        galseds = db.calc_fnu_sed_fast(spec[galid], filcurves)
        seds.append(np.array(galseds))
    return seds


def makespec_fast(sfht, Zt, timeax, zval, 
                  sp = db.mocksp, cosmo=db.cosmo,
                  interpolator = 'pchip', 
                  peraa=False, use_med_Z = True):

    if use_med_Z == True:
        sp._zcontinuous = 1
        sp.params['logzsol'] = np.nanmedian(Zt)
        sp.set_tabular_sfh(timeax, sfht)
    else:
        sp._zcontinuous = 3
        sp.set_tabular_sfh(timeax, sfht, Z=(0.0142 * 10**Zt))
    
    sp.params['zred'] = zval
    lam, spec = sp.get_spectrum(tage = cosmo.age(zval).value+1e-4, peraa = peraa)
    spec_ujy = convert_to_microjansky(spec, zval, cosmo)


    return lam, spec_ujy


def uncompress_sfh(sfhtuple, Ztuple, zval, cosmo, sfr_Myr = 10, interpolator=None):
    tempsfh, temptime = tuple_to_sfh_new(sfhtuple, zval, cosmo = cosmo, interpolator = interpolator, 
                                                     sfr_Myr=sfr_Myr, vb=False)
    tempZt, _ = tuple_to_sfh_new(Ztuple, zval, cosmo = cosmo, interpolator = interpolator, 
                                                 sfr_Myr=sfr_Myr, vb=False)
    tempZt = np.log10(tempZt)
    tempZt[tempZt<-2.5] = -2.5
    return tempsfh, tempZt, temptime