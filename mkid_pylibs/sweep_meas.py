import numpy as np
import os

import matplotlib.pyplot as plt
import matplotlib.cm as cm
from lmfit.parameter import Parameter as pm

from . import fit
from . import misc

class SweepMeas(object):
    def __init__(self):
        self._filenamelist = None
        self._paramlist = None
        pass

    def __len__(self):
        return len(self._paramlist)
    def __getitem__(self,key):
        return self.get_vals(key)

    def _adjustlen(self):
        if self._paramlist is None: return
        maxlen = np.max([len(v) for v in self._paramlist.values()])
        for k in self._paramlist.keys():
            if len(self._paramlist[k]) < maxlen:
                for _ in range(maxlen - len(self._paramlist[k])):
                    self._paramlist[k] += [pm(name=k,value=None)]
        if not all([len(v) == maxlen for v in self._paramlist.values()]):
            raise RuntimeError(f'Unmatched length of list --> {[len(v) for v in self._paramlist.values()]}')
        while all([v[-1].value == -1*float('inf') or v[-1].value == 0.0 for v in self._paramlist.values()]):
            for k in self._paramlist.keys():
                self._paramlist[k] = self._paramlist[k][:-1]

    def print_params(self):
        for kk in self._paramlist.keys():
            print(kk)

    def sweepVar(self, vallist, name, renew=False):
        '''
        Register parameter directly from list.
        :param list[] vallist: 1D list of sweeping variables.
        :param str name: Variable name.
        :param bool renew: If True, old list will be removed.
        '''
        if self._paramlist is None:
            self._paramlist = {}

        if renew or name not in self._paramlist: self._paramlist[name] = []
        self._paramlist[name] += [pm(name=name,value=vv) for vv in vallist]
        self._adjustlen()

    def sweepVar_readfromfile(self, ftype, inpath_swp, frq, index, cont_index=None,
                              inpath_tod=None, dopsd=False, dotrg=False,
                              renew=False, plotted=False, label='', nmax=None,**kws):
        '''
        Register paramters from IQ data
        :param ftype: File format of infile. One of ['rhea', 'riken', 'vna'].
        :param inpath_swp: 1D list of swp file names.
        :param frq: LO-offset. This should be str with ending "GHz"/"MHz"/"kHz"/"Hz" or float [Hz].
        :param index: Index-number for on-resonance.
        :param cont_index: Index-number for off-resonance.
        :param todtype: IQdata type to be calclated as the PSD.
        :param renew: If True, old list will be removed.
        :param plotted: If not False, fitting results in all files will be plotted, and saved as [plotted].pdf.
        :param label: If not None and plotted is not False, each plot curve is labelled with self._paramlist[label]. If unit is needed, it should be 'LABEL[UNIT]'.
        :param nmax: If specified, the number of data is up to `nmax` for TODdata.
        :param **kws: Fitting options to be put in kidfit.fit_onepeak() function.
        '''
        if self._paramlist is None: self._paramlist = {}
        if self._filenamelist is None: self._filenamelist = {}

        if type(inpath_swp) is str:
            inpath_swp = [inpath_swp]
        if type(inpath_tod) is str:
            inpath_tod = [inpath_tod]
        if len(inpath_swp) > len(inpath_tod):
            inpath_tod += [None]*(len(inpath_swp)-len(inpath_tod))
        if len(inpath_tod) > len(inpath_swp):
            inpath_swp += [None]*(len(inpath_tod)-len(inpath_swp))

        flen = len(inpath_swp)

        if not hasattr(index,'__iter__'):
            index = [index] * flen
        if not hasattr(cont_index,'__iter__'):
            cont_index = [cont_index] * flen

        ulabel=''
        vlabel=[0]*flen
        if '[' in label and ']' in label:
            ulabel = label[label.rindex('[')+1:label.rindex(']')]
            label = label[:label.rindex('[')]
            vlabel = self.get_vals(label)

        if plotted:
            colors = cm.rainbow(np.linspace(0, 1, flen))
            fig = None
            ax = None

        swp = None 
        cswp = None
        tod = None
        ctod = None
        from .readfile import readfile_swp,readfile_tod
        from .kidana import KidAnalyzer
        for ii in range(flen):
            fname_swp = inpath_swp[ii]
            fname_tod = inpath_tod[ii]
            ind  = index[ii]
            cind = cont_index[ii]
            if fname_swp is not None:
                print(f'{ii:02d} ({vlabel[ii]:8.3f} {ulabel}) SWP: {os.path.basename(fname_swp) if fname_swp is not None else None}')
                datalist = readfile_swp(ftype,fname_swp,-1,frq,**kws)
                swp = datalist[ind]
            if fname_tod is not None:
                print(f'{ii:02d} ({vlabel[ii]:8.3f} {ulabel}) TOD: {os.path.basename(fname_tod) if fname_tod is not None else None}')
                datalist = readfile_tod(ftype,fname_tod,-1,frq,**kws)
                tod = datalist[ind]
                if cind is not None:
                    ctod = datalist[cind]

            kid = KidAnalyzer(swp=swp, tod=tod, ctod=ctod)
            kid.fitIQ(silent=True)

            if fname_swp is not None:
                if not 'swp' in self._filenamelist: self._filenamelist['swp'] = []
                if ii==0 and renew: self._filenamelist['swp'] = []
                self._filenamelist['swp'].append(fname_swp)
                for pp,vv in kid.swp.fitresult.params.items():
                    pp = 'swp_'+pp
                    if pp not in self._paramlist: self._paramlist[pp] = []
                    self._paramlist[pp].append(vv)

            if fname_tod is not None and (dopsd or dotrg):
                ## PSD
                if dopsd:
                    varname = ""
                    if 'phs' in dopsd or 'phase' in dopsd:
                        varname = "phs"
                    if 'amp' in dopsd or 'amplitude' in dopsd:
                        varname = "amp"
                    kid.calcpsd(dofit=True,varname=varname,fitsilent=True)

                    for pp,vv in kid.psdfitresult.params.items():
                        pp = 'psd'+varname+'_'+pp
                        if pp not in self._paramlist: self._paramlist[pp] = []
                        if ii==0 and renew: self._paramlist[pp] = []
                        self._paramlist[pp].append(vv)

                if dotrg:
                    try:
                        kid.fittrig(silent=True)
                    except:
                        print("fit failure. skip this.")
                        continue
                    for pp,vv in kid.trgfitresult.params.items():
                        pp = 'trg_'+pp
                        if pp not in self._paramlist: self._paramlist[pp] = []
                        if ii==0 and renew: self._paramlist[pp] = []
                        self._paramlist[pp].append(vv)

            if plotted:
                if vlabel is not None:
                    label = f'{ii:02d}: {vlabel[ii]:.2f}'
                    if ulabel != '':
                        label += f' {ulabel},'
                from . import plotter
                fig,ax = plotter.plotSwpTOD(kid.swp,kid.tod,psddata = kid.psd[1:] if kid.tod is not None else None,
                                            onlyswp=True,todav=True,notfit=True,color=colors[ii], label=label, fig=fig, ax=ax)

        if plotted:
            return fig,ax

    def get_vals(self, varname, invert=False, log=False, valrange=None, params=None):
        vs,evs = self.get_vals_errs(varname,invert,log,valrange,params)
        return vs

    def get_vals_errs(self, var, invert=False, log=False, valrange=None, params=None):
        if params is None:
            params = self._paramlist

        vs = None
        evs = None
        if type(var) is list:
            vs = np.array(var)
        elif not var in params:
            return None,None
        else:
            vs = np.array([pp.value for pp in params[var]])
            evs = np.array([pp.stderr for pp in params[var]])
            if any([ee is None for ee in evs]):
                evs = None

        if invert:
            if evs is not None:
                evs = evs/(vs*vs)
            vs = 1./vs
        if log:
            if evs is not None:
                evs = np.abs(evs/vs)
            vs = np.log10(vs)
        if valrange is not None:
            if evs is not None:
                evs = evs[valrange]
            vs = vs[valrange]

        return vs,evs

    def plot(self, x_var, y_var, noerrs=False,
             x_invert=False, x_log=False, x_range=None, xname = 'x',
             y_invert=False, y_log=False, y_range=None, yname = 'y'):
        xs,exs = self.get_vals_errs(x_var,invert=x_invert,log=x_log,valrange=x_range)
        if xname=='x' and type(x_var) is str:
            xname = x_var

        ys,eys = self.get_vals_errs(y_var,invert=y_invert,log=y_log,valrange=y_range)
        if yname=='y' and type(y_var) is str:
            yname = y_var

        fig, axs = plt.subplots()
        fig.set_size_inches(10,6)
        if noerrs:
            axs.plot(xs,ys,'bo')
        else:
            axs.errorbar(xs,ys,xerr=exs,yerr=eys,fmt='bo',ecolor='b')
        axs.set_xlabel(xname)
        axs.set_ylabel(yname)
        axs.grid(ls=':')
        fig.tight_layout()

        return fig, axs

    def hist_var(self, var, invert=False, log=False, valrange=None, name='x', **kws):

        xs,exs = self.get_vals_errs(var,invert=invert,log=log,valrange=valrange)
        if name=='x' and type(var) is str:
            name = var

        fig, axs = plt.subplots()
        fig.set_size_inches(10,6)
        axs.hist(xs,**kws)
        axs.set_xlabel(name)
        axs.set_ylabel('# of samples')
        axs.grid(ls=':')
        fig.tight_layout()

        return fig, axs

    def fit_linear(self, x_var,y_var,
                   x_invert=False, x_log=False, x_range=None, x_name = None,
                   y_invert=False, y_log=False, y_range=None, y_name = None,
                   **kws):
        '''
        Fit linear function.
        :param x_var: Fitting variable name on x-axis.
        :param y_var: Fitting variable name on y-axis.
        :param x_name: x name used if plotted is not False.
        :param x_invert: If True, x_var --> 1/x_var
        :param x_log: If True, x_var --> np.log(x_var)
        :param x_range: fit range for x_var
        :param y_name: y name used if plotted is not False.
        :param y_invert: If True, y_var --> 1/y_var
        :param y_log: If True, y_var --> np.log(y_var)
        :param y_range: fit range for y_var
        :param kws: fitting options. see misc.linear_fit()
        '''

        xs,exs = self.get_vals_errs(x_var,invert=x_invert,log=x_log,valrange=x_range)
        if x_name is None:
            x_name = x_var
            kws['xname'] = x_name

        ys,eys = self.get_vals_errs(y_var,invert=y_invert,log=y_log,valrange=y_range)
        if y_name is None:
            y_name = y_var
            kws['yname'] = y_name

        return misc.linear_fit(xs,ys,exs,eys,**kws)

def anasweep(ftype, inpath_swp, lo, index=0, cont_index=None, inpath_tod=None, temps=None, V=None):
    '''
    Wrapper for temperature sweep measurements.
    :param str ftype: one of ['rhea', 'riken', 'vna']
    :param list[str] inpath_swp: filename for Swpdata (one filename for each temperature)
    :param list[str] inpath_tod: filename for TODdata (one filename for each temperature)
    :param int index: index number for on-resonance
    :param int cont_index: index number for off-resonance
    :param list[float] temps: 1D list of temperature
    :param float V: KID volume [um3]
    '''
    tswp = SweepMeas()
    if temps is not None:
        tswp.sweepVar(temps, 'temp')
        tswp.sweepVar(misc.convert_Kelvin_nqp(np.array(temps)), 'nqp')
        if V is not None:
            tswp.sweepVar(misc.convert_Kelvin_nqp(np.array(temps)) * V, 'Nqp')
    tswp.sweepVar_readfromfile('rhea',inpath_swp,lo,index,cont_index,inpath_tod,'psd_phs',label='temp[K]',nmax=1e5, plotted=True)
    tswp.sweepVar(misc.convert_tqp_nqp(tswp.get_vals('psd_phs_tqp')), 'rolloff_nqp')
    tswp.sweepVar(misc.convert_tqp_nqp(tswp.get_vals('psd_phs_tqp'))*V, 'rolloff_Nqp')
    return tswp
