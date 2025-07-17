# for offreso (loopback) analysis

def data_set(dir_path, lo):
    ## SG FREQ
    #lo = "5.51GHz"
    ## TODDATASET
    tfp = glob.glob(dir_path+'/tod*.rawdata')
    print(tfp)
    todset = klib.multreadfile_tod('rhea',tfp,lo)
    return todset


def raw_rwm(todset):
    
    rate = [itod.rate for itod in todset[0]]
    freq= todset[0][0].f
    if len(todset) == 1:
        iq = [itod.iq for itod in todset[0]]
        # scale
        siq = [iiq / np.mean(np.abs(iiq)) for iiq in iq]
        # rotation
        sita= np.average(np.arctan2(siq[0].imag,siq[0].real))
        rsiq = [isiq* np.exp(-1j*sita) for isiq in siq]
        dic = {'data':iq, 'sdata' : siq, 'rsdata' : rsiq}
    
    
    
#if len(todset) == 2:
    else:
        iq = [itod.iq for itod in todset[0]]
        biq = [itod.iq for itod in todset[1]]
        #rotation
        rawsita = np.average(np.arctan2(iq[0].imag,iq[0].real))
        bsita= np.average(np.arctan2(biq[0].imag,biq[0].real))
        riq = [iiq* np.exp(-1j*rawsita) for iiq in iq]
        rbiq = [ibiq* np.exp(-1j*bsita) for ibiq in biq]
        # modify
        miq = [None]*len(todset[0])
        for i, (itod, citod) in enumerate(zip(todset[0], todset[1])):
            miq[i] = klib.modIQ.subtract_blindtone(itod.iq, citod.iq)
        mriq = [None]*len(todset[0])
        for i, (iriq, irbiq) in enumerate(zip(riq,rbiq)):
            mriq[i] = klib.modIQ.subtract_blindtone(iriq, irbiq)
            
        # scale amp
        smiq = [imiq / np.mean(np.abs(imiq)) for imiq in miq]
        siq = [iiq / np.mean(np.abs(iiq)) for iiq in iq]
        bsiq = [iiq / np.mean(np.abs(iiq)) for iiq in biq]
        # rotation
        sita= np.average(np.arctan2(smiq[0].imag,smiq[0].real))
        rsiq = [isiq* np.exp(-1j*sita) for isiq in siq]
        rsmiq = [ismiq* np.exp(-1j*sita) for ismiq in smiq]

        dic = {'data':iq,'rdata':riq, 'rbdata':rbiq, 'mrdata':mriq, 'bdata':biq,
               'bsdata':bsiq, 'mdata':miq, 'sdata' : siq, 'smdata' : smiq,
               'rsdata' : rsiq, 'rsmdata' : rsmiq, 'freq': freq, 'rate':rate}
    return dic

def con_psd(psd):
    for i, (ifreq, iamp, iphase) in enumerate(zip(psd['freq'], psd['amp'], psd['phase'])):
        if i == 0:
            freq = ifreq
            amp = iamp
            phase = iphase
        else :
            index = np.where(freq[-1]<ifreq)
            freq = np.append(freq, ifreq[index])
            amp = np.append(amp, iamp[index])
            phase = np.append(phase, iphase[index])
    ret = {'freq':freq[1:],'amp':amp[1:], 'phase':phase[1:]}
    return ret

def cal_psd(iq, todset):
    frq =[None]*len(todset[0])
    amp= [None]*len(todset[0])
    phase= [None]*len(todset[0])

    for i, iq in enumerate(iq):
        frq[i], amp[i], phase[i] = klib.psd.calc_psd_welch(np.abs(iq), np.arctan2(iq.imag, iq.real), rate = todset[0][i].rate)
    dic = {'freq':frq, 'amp':amp, 'phase' : phase}
    ret = con_psd(dic)
    return ret

'''
def plot_psd(psd, title = None):
    plt.figure()
    frq = psd['freq']
    amp = psd['amp']
    phase = psd['phase']
    labels = ['amp', 'phase']
    for ifrq, iamp, iphase in zip(frq, amp, phase):
        plt.plot(ifrq[1:], 10*np.log10(iamp[1:]), '.:', c='red', label = labels[0])
        plt.plot(ifrq[1:], 10*np.log10(iphase[1:]), '.:', c= 'blue', label = labels[1])
    plt.xscale('log')
    if title is not None:
        plt.title(title)
    plt.legend(labels, bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0, fontsize=18)
'''

def plot_psd(psd, title = None):
    plt.plot(psd['freq'], 10*np.log10(psd['amp']), '.:', c='red', label = 'amp')
    plt.plot(psd['freq'], 10*np.log10(psd['phase']), '.:', c= 'blue', label = 'phase')
    plt.xscale('log')
    if title is not None:
        plt.title(title)
#    plt.legend(bbox_to_anchor=(1, 1), loc='upper right', borderaxespad=0, fontsize=18)
    plt.legend(loc='best')
        
def plot_iq(iq):
#    fig, ax = plt.subplots(nrows=1, ncols=4, figsize=(20,4))
 #   sps = ['1kSPS','10kSPS', '100kSPS', '1000kSPS']
    fig, ax = plt.subplots(nrows=1, ncols=3, figsize=(16,4))
    sps = ['1kSPS','10kSPS', '100kSPS']
    for i, (iiq,isps) in enumerate(zip(iq,sps)):
#        ax[i].plot(iiq.real, iiq.imag)
        ax[i].plot(iiq.real[::5], iiq.imag[::5])
        ax[i].set_title(isps)
        ax[i].axis('equal') 
    plt.subplots_adjust(wspace=0.4)
        
def ana(dir_path, lo, kind = 'rsmdata'):
    todset = data_set(dir_path, lo)
    iq = raw_rwm(todset)
    if 'rsmdata' == kind:
        psd = cal_psd(iq['rsmdata'], todset)
    elif 'smdata' == kind:
        psd = cal_psd(iq['smdata'], todset)
    elif 'rsdata' == kind:
        psd = cal_psd(iq['rsdata'], todset)
    elif 'sdata' == kind:
        psd = cal_psd(iq['sdata'], todset)
    elif 'rdata' == kind:
        psd = cal_psd(iq['rdata'], todset)
    elif 'rbdata' == kind:
        psd = cal_psd(iq['rbdata'], todset)
    elif 'mrdata' == kind:
        psd = cal_psd(iq['mrdata'], todset)
    elif 'mdata' == kind:
        psd = cal_psd(iq['mdata'], todset)
    elif 'bsdata' == kind:
        psd = cal_psd(iq['bsdata'], todset)
    elif 'bdata' == kind:
        psd = cal_psd(iq['bdata'], todset)
    elif 'data' == kind:
        psd = cal_psd(iq['data'], todset)    
    plot_psd(psd, title = kind)
    dic = {'iq':iq, 'psd': psd}
    return dic

def plot_2d(dir_path, tod_data):
    #swp_data = read_rhea_mulswp(dir_path + 'swp.rawdata')
    swp_data = read_rhea_swp(glob.glob(dir_path+'/swp*.rawdata')[0])
    plt.figure()
    #plt.plot(swp_data[0]['I'], swp_data[0]['Q'], zorder=1)
    plt.plot(swp_data['I'], swp_data['Q'], zorder=1)
    #plt.plot(data1['iq']['data'][2].real, data1['iq']['data'][2].imag)
    for i, iq in enumerate(tod_data['iq']['data']):
        plt.plot(iq.real, iq.imag, zorder = 5-i)
    plt.axis('equal')
    print('-------amp_rad-------')
    print(np.average(swp_data['amp_rad']))
    
def plot_amp(dir_path, tod_data, lo):
    swp_data = read_rhea_mulswp(dir_path + 'swp.rawdata')
    plt.figure()
    plt.plot(swp_data[0]['freq']+lo, swp_data[0]['amp_rad'])
    #plt.plot(data1['iq']['data'][2].real, data1['iq']['data'][2].imag)
#    for i, iq in enumerate(tod_data['iq']):
    plt.axvline(x=tod_data['iq']['freq'], ymin=0.1, ymax=0.9, c='r')
    #plt.axis('equal')
    print('-------amp_rad-------')
    print(np.average(swp_data[0]['amp_rad']))
    
def plot_sweep(path):
    swppath = glob.glob(path + '/swp*data')
    swpdata = read_rhea_mulswp(swppath[0])
    plt.plot(swpdata[0]['freq'], swpdata[0]['amp_rad'])
