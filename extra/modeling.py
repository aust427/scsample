

def uncompress_sfh(sfhtuple, Ztuple, zval, sfr_Myr = 10, interpolator=None):
    tempsfh, temptime = tuple_to_sfh_new(sfhtuple, zval, interpolator = interpolator, 
                                                     sfr_Myr=sfr_Myr, vb=False)
    tempZt, _ = tuple_to_sfh_new(Ztuple, zval, interpolator = interpolator, 
                                                 sfr_Myr=sfr_Myr, vb=False)
    tempZt = np.log10(tempZt)
    tempZt[tempZt<-2.5] = -2.5
    return tempsfh, tempZt, temptime


def makespec_fast(sfht, Zt, timeax, zval, 
                  sp = db.mocksp, cosmo=cosmo,
                  interpolator = interpolator, 
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
