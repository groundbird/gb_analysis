import matplotlib.pyplot as plt
import numpy as np
from . import misc

misc.def_style()

def plotSwp(swpdata, fig=None, ax=None, kind=None, legend=True, title='', is_compact=True, **kws):
    '''
    Plot Sweep IQ data.
    ax[0]: freq vs amplitude
    ax[1]: freq vs phase
    ax[2]: I vs Q
    '''
    if ax is None:
        fsz = (14,6)
        nc = 4
        nr = 4
        if legend == 'outside':
            fsz = (16,6)
            nc+=1
        if not is_compact:
            fsz = (20,6)
            nc=nc+2
            nr=1
        gsz = (nr,nc)
        fig,ax = misc.make_figure(fsz=fsz,nc=1,nr=1,dpi=100)
        ax = [0]*3
        if is_compact:
            ax[0] = plt.subplot2grid(gsz, (0, 0), colspan=2, rowspan=3)
            ax[1] = plt.subplot2grid(gsz, (3, 0), colspan=2, rowspan=1, sharex=ax[0])
            ax[2] = plt.subplot2grid(gsz, (0, 2), colspan=2, rowspan=4)
        else:
            ax[0] = plt.subplot2grid(gsz, (0, 0), colspan=2, rowspan=1)
            ax[1] = plt.subplot2grid(gsz, (0, 2), colspan=2, rowspan=1, sharex=ax[0])
            ax[2] = plt.subplot2grid(gsz, (0, 4), colspan=2, rowspan=1)
        if legend == 'outside':
            ax.append(plt.subplot2grid(gsz, (0, nc-1), colspan=1, rowspan=2))
            ax[3].axis('off')
    elif len(ax)<3:
        raise RuntimeError(f'Not eonugh # of axis ({len(ax)}) < 3')

    x = swpdata.f/1e9
    y,kind = misc.get_data_type(swpdata,kind)

    ax[0].plot(x, y.amplitude, **kws)
    ax[1].plot(x, y.corphase, **kws)
    ax[2].plot(y.i, y.q, **kws)

    ax[0].set_xlabel('frequency [GHz]')
    ax[0].set_ylabel('amplitude')
    ax[1].set_xlabel('frequency [GHz]')
    ax[1].set_ylabel('phase')
    ax[2].set_xlabel('I')
    ax[2].set_ylabel('Q')

    if legend and 'label' in kws:
        if legend == 'outside':
            ax[2].legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        elif legend:
            ax[0].legend(loc='best')
            ax[1].legend(loc='best')
            ax[2].legend(loc='best')

    xmin,xmax = ax[2].get_xlim()
    ymin,ymax = ax[2].get_ylim()
    ax[2].set_xlim(min(xmin,ymin),max(xmax,ymax))
    ax[2].set_ylim(min(xmin,ymin),max(xmax,ymax))
    ax[2].set_aspect('equal')

    if fig is not None:
        misc.adj_figure(fig,title=title)
        if is_compact:
            fig.subplots_adjust(hspace=0.001)
            plt.setp(ax[0].get_xticklabels(), visible=False)
            ax[0].set_xlabel('')
        #else:
        #    fig.tight_layout()

    return fig,ax

def plotPSD(data, ft=None, fig=None, ax=None, ampcol='r', phscol='b', fitcol='k', title='', label='', leg=True, **kws):
    if ax is None:
        fig,ax = misc.make_figure()
    if label != '':
        label = label+': '

    if 'ls' not in kws:
        kws['ls'] = ':'
    if 'alpha' not in kws:
        kws['alpha'] = 0.7
    if 'marker' not in kws:
        kws['marker'] = '.'

    if ampcol is not None:
        ax.semilogx(data.f, np.log10(data.amplitude)*10, lw=2, label=f'{label}amp', color=ampcol, **kws)
    if phscol is not None:
        ax.semilogx(data.f, np.log10(data.phase)*10, lw=2, label=f'{label}phase', color=phscol, **kws)
    if ft and fitcol is not None:
        kws['marker'] = ''
        kws['ls'] = '-'
        ax.semilogx(data.f,ft.eval(data.f), lw=2, label=f'{label}fit', color=fitcol, **kws)
    if leg:
        ax.legend(loc='best',fontsize='x-large')

    ax.set_ylabel('PSD [dBc/Hz]')
    ax.set_xlabel('frequency [Hz]')
    minf = data.f.min()
    minf = 10**int(np.log10(minf)-2) * 7
    maxf = data.f.max()
    maxf = np.min([10**int(np.log10(maxf) + 2) * 7 , 7e5])
    ax.set_xlim(minf,maxf)

    misc.adj_figure(fig,title=title)
    return fig, ax

def plotNEP(data, ft=None, fig=None, ax=None, ampcol='r', phscol='b', title='', label='', leg=True, **kws):
    if ax is None:
        fig,ax = misc.make_figure()
    if label != '':
        label = label+': '

    if ampcol is not None:
        ax.loglog(data.f, data.amplitude, marker='.', lw=1, ls=':', label=f'{label}amp', color=ampcol, **kws)
    if phscol is not None:
        ax.loglog(data.f, data.phase, marker='.', lw=1, ls=':', label=f'{label}phase', color=phscol, **kws)
    if leg:
        ax.legend(loc='best',fontsize='x-large')

    ax.set_ylabel('NEP [/Hz]')
    ax.set_xlabel('Frequency [Hz]')
    minf = data.f.min()
    minf = 10**int(np.log10(minf)) * 7
    maxf = data.f.max()
    maxf = np.min([10**int(np.log10(maxf) + 2) * 7 , 7e5])
    ax.set_xlim(minf,maxf)
    ax.grid(color='grey',linestyle=':')

    misc.adj_figure(fig,title=title)
    return fig, ax

def plotTOD(toddata, fig=None, ax=None, kind=None, trgft=None, swpdata=None, x_time=True, legend=True, title='', **kws):

    if ax is None:
        fsz = (26,3)
        nc = 14
        nr = 2
        gsz = (nr,nc)
        fig,ax = misc.make_figure(fsz=fsz,nc=1,nr=1,dpi=100)
        ax = [0]*5
        ax[0] = plt.subplot2grid(gsz, (0,   0), colspan=4, rowspan=1)
        ax[1] = plt.subplot2grid(gsz, (1,   0), colspan=4, rowspan=1, sharex=ax[0])
        ax[2] = plt.subplot2grid(gsz, (0,   4), colspan=4, rowspan=2, sharex=ax[0])
        ax[3] = plt.subplot2grid(gsz, (0,   8), colspan=4, rowspan=2, sharex=ax[0])
        ax[4] = plt.subplot2grid(gsz, (0,   12), colspan=2, rowspan=2)
    elif len(ax)<4:
        raise RuntimeError(f'Not eonugh # of axis ({len(ax)}) < 4')

    if x_time:
        x = toddata.time
    else:
        x = np.arange(len(toddata.x))
    y,kind = misc.get_data_type(toddata,kind)
    default_col = True
    if 'color' in kws:
        default_col = False

    if default_col: kws['color'] = 'black'
    ax[0].plot(x, y.i, **kws)
    ax[0].set_ylabel('I')
    plt.setp(ax[0].get_xticklabels(), visible=False)

    ax[1].plot(x, y.q, **kws)
    ax[1].set_ylabel('Q')
    ax[1].set_xlabel('time[s]')

    if default_col: kws['color'] = 'red'
    ax[2].plot(x, y.amplitude, **kws)
    ax[2].set_ylabel('amplitude')
    ax[2].set_xlabel('time[s]')

    if default_col: kws['color'] = 'blue'
    ax[3].plot(x, y.corphase, **kws)
    if trgft:
        if default_col: kws['color'] = 'orange'
        ax[3].plot(x, trgft.eval(x), ls=':', lw=4, color='orange', marker='')
    ax[3].set_ylabel('phase')
    ax[3].set_xlabel('time[s]')

    if swpdata is not None:
        if 'rw' in kind: iq = swpdata.fitdata.rwdata.iq
        else:            iq = swpdata.fitdata.iq
        ax[4].plot(iq.real, iq.imag, color='orange', ls=':', lw=2, marker='')
    if default_col: kws['color'] = 'blue'
    ax[4].plot(y.i, y.q, **kws)
    ax[4].set_ylabel('Q')
    ax[4].set_xlabel('I')
    xmin,xmax = ax[4].get_xlim()
    ymin,ymax = ax[4].get_ylim()
    ax[4].set_xlim(min(xmin,ymin),max(xmax,ymax))
    ax[4].set_ylim(min(xmin,ymin),max(xmax,ymax))
    ax[4].set_aspect('equal')

    misc.adj_figure(fig,title=title)
    fig.subplots_adjust(hspace=0.001)

    return fig,ax

def plotTOD2(todlist,psd=None,psdft=None,fig=None,ax=None,kind=None,legend=True,title='', ampcol='r', phscol='b', fitcol='k',label='',skimmedlegend=True):

    import matplotlib.pyplot as plt
    if ax is None:
        fsz = (20,6)
        nc = 5
        nr = len(todlist)*3
        if legend=='outside':
            fsz = (24,6)
            nc=nc+1
        gsz = (nr,nc)
        fig,ax = misc.make_figure(fsz=fsz,nc=1,nr=1,dpi=100)
        ax = [0]*(len(todlist)*2+1)
        tmpdict = {}
        tmpdict['colspan'] = 3
        tmpdict['rowspan'] = 1
        for ii in range(len(todlist)):
            if ii != 0: tmpdict['sharey'] = ax[0]
            ax[ii*2]   = plt.subplot2grid(gsz, (ii*3,0), **tmpdict)
            if ii != 0: tmpdict['sharey'] = ax[1]
            ax[ii*2+1]   = plt.subplot2grid(gsz, (ii*3+1,0), sharex=ax[ii*2], **tmpdict)

        ax[-1] = plt.subplot2grid(gsz, (0, tmpdict['colspan']), colspan=2, rowspan=len(todlist)*3)
        if legend == 'outside':
            ax.append(plt.subplot2grid(gsz, (0, 5), colspan=1, rowspan=len(todlist)*3))
            ax[-1].axis('off')
    elif len(ax) < (len(todlist)*2+1):
        raise RuntimeError(f'ERROR:: Need to have {(len(todlist)*2+1)} elements for subplots.')

    if label != '': label+=': '

    ii = 0
    for toddata in todlist:
        x = toddata.time
        y,kind = misc.get_data_type(toddata,kind)
        ax[ii*2].plot(x,y.amplitude,color=ampcol,label=label+f'{misc.tosi(toddata.rate,ret_int=True)}SPS, {misc.get_datatype_description(kind)}')
        ax[ii*2].set_ylabel('amp')
        ax[ii*2+1].plot(x,y.corphase,color=phscol,label=label+f'{misc.tosi(toddata.rate,ret_int=True)}SPS, {misc.get_datatype_description(kind)}')
        ax[ii*2+1].set_ylabel('phs')

        ax[ii*2].set_xlim(toddata.time.max()*(-0.05),toddata.time.max()*1.5)
        ax[ii*2].legend(loc='best',fontsize='small')
        ax[ii*2+1].legend(loc='best',fontsize='small')
        ii += 1

    for ii in range(len(todlist)):
        plt.setp(ax[ii*2].get_xticklabels(), visible=False)
    ax[(len(todlist)-1)*2+1].set_xlabel('time[s]')

    plotPSD(psd,psdft,fig,ax[(len(todlist)-1)*2+2],ampcol,phscol,fitcol,label=label,leg=True)

    misc.adj_figure(fig,title=title)
    fig.subplots_adjust(hspace=0.001)

    return fig,ax

    # if skimmedlegend:
    #     lines=[]
    #     labels=[]
    #     skimmedarray=[]
    #     for a in [ax[-2]]:
    #         h, l = a.get_legend_handles_labels()
    #         for hi, li in zip(h,l):
    #             if skimmedlegend:
    #                 keystr = li.split(':')[0]
    #                 if keystr in skimmedarray:
    #                     continue
    #                 else:
    #                     skimmedarray.append(keystr)
    #             lines.append(hi)
    #             labels.append(li)
    #     ax[-1].legend(handles=lines, labels=labels, loc='best')

def plotSwpTOD(swpdata, toddata, psddata=None, psdfitresult=None,
               onlyswp=False, todav=False, notfit=False,
               fig=None, ax=None, kind=None, legend=True, title='', **kws):

    if ax is None:
        fsz = (28,6)
        nc = 14
        nr = 4
        if legend == 'outside':
            fsz = (16,6)
            nc+=1
        gsz = (nr,nc)
        fig,ax = misc.make_figure(fsz=fsz,nc=1,nr=1,dpi=100)
        ax = [0]*4
        ax[0] = plt.subplot2grid(gsz, (0, 0),  colspan=4, rowspan=3)
        ax[1] = plt.subplot2grid(gsz, (3, 0),  colspan=4, rowspan=1, sharex=ax[0])
        ax[2] = plt.subplot2grid(gsz, (0, 4),  colspan=4, rowspan=4)
        ax[3] = plt.subplot2grid(gsz, (0, 8),  colspan=6, rowspan=4)
        if psddata is None and toddata is not None:
            ax += [0]*len(toddata)
            sharey = None
            for ii in range(len(toddata)):
                if ii!=0: sharey = ax[4]
                ax[ii+4] = plt.subplot2grid((len(toddata)*4,nc), (ii*4, 8),  colspan=6, rowspan=3, sharey=sharey)
        if legend == 'outside':
            ax.append(plt.subplot2grid(gsz, (0, nc-1), colspan=1, rowspan=2))
            ax[-1].axis('off')
    elif len(ax)<4:
        raise RuntimeError(f'Not eonugh # of axis ({len(ax)}) < 3')

    defaultlabel=False
    if 'label' not in kws:
        defaultlabel = True
    defaultcolor = False
    if 'color' not in kws:
        defaultcolor = True

    if defaultlabel: kws['label'] = 'data'
    if defaultcolor: kws['color'] = 'blue'

    ax[0].plot(swpdata.fGHz, swpdata.data.amplitude,marker='.',ls=':', **kws)
    ax[1].plot(swpdata.fGHz, swpdata.data.phase, marker='.',ls=':', **kws)
    ax[2].plot(swpdata.rwdata.i, swpdata.rwdata.q, marker='.',ls=':', **kws)

    if not notfit:
        if defaultlabel: kws['label'] = 'fitted'
        if defaultcolor: kws['color'] = 'red'
        ax[0].plot(swpdata.fGHz, swpdata.fitdata.data.amplitude, marker='', lw=2, **kws)
        ax[1].plot(swpdata.fGHz, swpdata.fitdata.data.phase, marker='', **kws)
        ax[2].plot(swpdata.fitdata.rwdata.i, swpdata.fitdata.rwdata.q, marker='', lw=2, **kws)

    ax[0].set_xlabel('frequency [GHz]')
    ax[0].set_ylabel('amplitude')
    ax[1].set_xlabel('frequency [GHz]')
    ax[1].set_ylabel('phase')
    ax[2].set_xlabel('I (rewind)')
    ax[2].set_ylabel('Q (rewind)')
    xmin,xmax = ax[2].get_xlim()
    ymin,ymax = ax[2].get_ylim()
    ax[2].set_xlim(min(xmin,ymin),max(xmax,ymax))
    ax[2].set_ylim(min(xmin,ymin),max(xmax,ymax))
    ax[2].set_aspect('equal')

    cmap = plt.get_cmap("tab10")
    if toddata is not None:
        for ii,tod in enumerate(reversed(toddata)):
            if not onlyswp:
                if not todav:
                    if defaultlabel: kws['label'] = f'{int(tod.rate/1e3)}kSPS'
                    if defaultcolor: kws['color'] = 'orange'
                    ax[0].plot([tod.fGHz]*len(tod.x),tod.data.amplitude,ls=':', **kws)
                    ax[2].plot(tod.rwdata.i,tod.rwdata.q,ls=':', **kws)
                    #if defaultcolor: kws['color'] = 'pink'
                    #ax[2].plot(tod.rwmdata.i,tod.rwmdata.q,ls=':', **kws)
                    ax[1].plot([tod.fGHz]*len(tod.x),tod.data.phase,ls=':', **kws)
                else:
                    if defaultlabel: kws['label'] = f'{int(tod.rate/1e3)}kSPS'
                    if defaultcolor: kws['color'] = cmap(ii)
                    ax[0].plot([tod.fGHz],[tod.data.amplitude.mean()],marker='o', lw=0, **kws)
                    ax[2].plot([tod.rwmdata.i.mean()],[tod.rwmdata.q.mean()],marker='o', lw=0, **kws)
                    if defaultlabel: kws['label'] = None
                    ax[1].plot([tod.fGHz],[tod.data.phase.mean()],marker='o', lw=0, **kws)
            if psddata is None:
                if defaultlabel: kws['label'] = f'{int(tod.rate/1e3)}kSPS'
                if defaultcolor: kws['color'] = cmap(ii)
                ax[ii+4].plot(tod.time,tod.rwmdata.phase, **kws)
                ax[ii+4].legend(loc='best')
                ax[ii+4].set_ylabel('phase')
        if psddata is None:
            ax[len(toddata)-1+4].set_xlabel('time [sec]')
        else:
            if defaultlabel: kws['label'] = ''
            if defaultcolor: kws['color'] = cmap(ii)
            plotPSD(psddata,ft=psdfitresult,fig=fig,ax=ax[3], label=kws['label'], ampcol='r' if defaultcolor else None, phscol='b' if defaultcolor else kws['color'], fitcol='k' if defaultcolor else kws['color'])

    ax[0].legend(loc='best')
    ax[2].legend(loc='best')

    if fig is not None:
        misc.adj_figure(fig,title=title)
        fig.subplots_adjust(hspace=0.001)
        plt.setp(ax[0].get_xticklabels(), visible=False)
        ax[0].set_xlabel('')

    return fig,ax

