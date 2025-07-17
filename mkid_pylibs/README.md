mkid_pylibs
====
MKID analysis libraries

## Description

This module includes:
* rawfile to sweep/tod data list
* peak finding in sweep data
* sweep fitting
* calculating PSD from tod data list (available with two tones)
* PSD fitting
* calculating NEP from PSD and given dfr/dNqp and d(1/Qr)/dNqp

The sweep fitting and data strage classes were derived from tomita_adapter.

## Required Version

Python 3.6.5+

## Version

If you want to old mkid_pylibs, please use tag`1.0`.


## Used Python Libraries/Modules (lowest version)

* ad
* builtins
* collections
* copy
* datetime
* functools
* glob
* json
* lmfit (0.9.11)
* math
* matplotlib (2.2.2)
* numpy (1.14.3)
* operator
* os
* pan
* pandas (0.23.0)
* parsefile
* pickle
* pprint
* re
* scipy　(1.1.0)
* struct
* sympy (1.1.1)
* sys
* time
* warnings

## Setup Version (to have a consistency with shonda)

## Short explanation for each file

* `__init__.py`

### Data class
* `Drawer.py`: basic class for supporting plot function for 2D-array
* `GeneralArray.py`: function to define the custom named_array
* `FreqTimeArray.py`: frequency and time data class
* `KidGlobalResponse.py`: KID IQ data class
* `BaseMultiData.py`: basic class for 2D-array
* `Swpdata.py`: frequency + KID IQ class for KID sweep
* `TODdata.py`: time + KID IQ class for KID TOD
* `PSDdata.py`: frequency + KID amp/phs named_array for KID PSD

### IQ fitting
* `kidfit.py`: IQ fitting interface
* `fitters.py`: IQ fitting functions
* `peak_search.py`: simple peak search tool with Lorenzian

### File treatment
* `readfile.py`: getting one KID data from a given file
* `multreadfile.py`: getting all KID data from a given file
* `parsefile.py`: written how to treat each file type

* `logreader.py`: getting information from logger file (under development)
* `rotary_reader.py`: (under development)

### Analysis
* `modIQ.py`: IQ modifications (on-off subtraction)
* `psd.py`: PSD calculation/analysis
* `nep.py`: NEP calculation
* `trig.py`: trigger analysis

* `kidana.py`: integrated class to analyze KIDs
* `sweep_meas.py`: analysis of variable sweeping measurement (such as temperature)

### Several plots
* `plotter.py`

### Others
* `misc.py`: miscellaneous useful functions

### Directories
* `rhea_comm/` : submodule to handle rawfiles taken by FPGA.
* `fit/` : general fitting interface using lmfit

## Example usage

### jupyter notebook

```bash
# you are now in [HOGE]/mkid_pylibs/
mv exam_code.ipynb ../
cd ../
git clone ssh://git@emperor.riken.jp:52222/shonda/data_exams_mkid_pylibs.git
jupyter notebook exam_code.ipynb
```

### python

```bash
import mkid_pylibs as klib

## SWPDATASET
swp_fp = "data_exams_mkid_pylibs/data_ex_optSRONKIDS/data_2ch__2018-1011-140813/swp__0.0uW_RF1amp22.0_RF2amp25.0_KID06_5.31GHz_+053.779MHz_+059.779MHz_+000.010MHz_2018-1011-140826.rawdata"

## TODDATASET
tod_fp = glob.glob(os.path.dirname(swp_fp)+'/tod_*.rawdata')

## SG FREQ
lo = "5.31GHz" # OR 4.99e9

## EACH DATA HAS MULTIPLE IQ DATA
## ASSUMING even# data (0,2,...) = KID-tone, odd# data (1,3,...) = off-tone of one before KID

kid = klib.kidana.KidAnalyzer()

kid.swpfromfile('rhea',swp_fp,lo,0)
kid.swp.fitIQ(nfwhm=3, fitter='gaolinbg')

kid.swp.fitresult.report()

kid.todfromfile('rhea',tod_fp,lo,1,0)
kid.moddata()

kid.calcpsd()
kid.psdfitresult.report()

## TEMPERATURE SWEEP

T = 0.20968781221779434 #[K] for kid

dn = 'data_exams_mkid_pylibs/data_ex_optSRONKIDS/'
Tswp_files = [
    dn+'data_2ch__2018-1011-170744/swp__0.0uW_RF1amp22.0_RF2amp25.0_KID06_5.31GHz_+053.772MHz_+059.772MHz_+000.010MHz_2018-1011-170757.rawdata',
    dn+'data_2ch__2018-1011-165242/swp__0.0uW_RF1amp22.0_RF2amp25.0_KID06_5.31GHz_+053.758MHz_+059.758MHz_+000.010MHz_2018-1011-165255.rawdata',
    dn+'data_2ch__2018-1011-142339/swp__0.0uW_RF1amp22.0_RF2amp25.0_KID06_5.31GHz_+053.778MHz_+059.778MHz_+000.010MHz_2018-1011-142352.rawdata',
    dn+'data_2ch__2018-1011-155116/swp__0.0uW_RF1amp22.0_RF2amp25.0_KID06_5.31GHz_+053.652MHz_+059.652MHz_+000.010MHz_2018-1011-155128.rawdata',
    dn+'data_2ch__2018-1011-150954/swp__0.0uW_RF1amp22.0_RF2amp25.0_KID06_5.31GHz_+053.680MHz_+059.680MHz_+000.010MHz_2018-1011-151007.rawdata',
    dn+'data_2ch__2018-1011-163736/swp__0.0uW_RF1amp22.0_RF2amp25.0_KID06_5.31GHz_+053.492MHz_+059.492MHz_+000.010MHz_2018-1011-163749.rawdata',
    dn+'data_2ch__2018-1011-160640/swp__0.0uW_RF1amp22.0_RF2amp25.0_KID06_5.31GHz_+053.579MHz_+059.579MHz_+000.010MHz_2018-1011-160652.rawdata',
    dn+'data_2ch__2018-1011-162209/swp__0.0uW_RF1amp22.0_RF2amp25.0_KID06_5.31GHz_+053.561MHz_+059.561MHz_+000.010MHz_2018-1011-162222.rawdata',
    dn+'data_2ch__2018-1011-140813/swp__0.0uW_RF1amp22.0_RF2amp25.0_KID06_5.31GHz_+053.779MHz_+059.779MHz_+000.010MHz_2018-1011-140826.rawdata',
    dn+'data_2ch__2018-1011-143903/swp__0.0uW_RF1amp22.0_RF2amp25.0_KID06_5.31GHz_+053.775MHz_+059.775MHz_+000.010MHz_2018-1011-143916.rawdata',
    dn+'data_2ch__2018-1011-145429/swp__0.0uW_RF1amp22.0_RF2amp25.0_KID06_5.31GHz_+053.733MHz_+059.733MHz_+000.010MHz_2018-1011-145442.rawdata',
]

Tswp_temps = [
    0.25074814317359834,
    0.2914201404975535,
    0.21767800139278862,
    0.3381137580080448,
    0.3296872909109462,
    0.3630744075579845,
    0.3530849082626251,
    0.3531586061770769,
    0.20968781221779434,
    0.2633075187992045,
    0.314944882291794,
]

tmpary = zip(Tswp_temps,Tswp_files)
Tswp_temps,Tswp_files = zip(*sorted(tmpary,reverse=True))

Ttod_files = [max(glob.glob(os.path.dirname(ff)+'/tod__*1000kSPS*.rawdata'), key=os.path.getctime) for ff in Tswp_files]

# input KID property
V = 135 #[um3]

tswp = klib.SweepMeas()
tswp.sweepVar(Tswp_temps, 'temp', renew=True)
tswp.sweepVar(klib.misc.convert_Kelvin_nqp(np.array(Tswp_temps)), 'nqp', renew=True)
tswp.sweepVar(klib.misc.convert_Kelvin_nqp(np.array(Tswp_temps)) * V, 'Nqp', renew=True)
tswp.sweepVar_readfromfile('rhea',Tswp_files,lo,0,1,inpath_tod=Ttod_files,label='temp[K]',dopsd='phs',plotted=False,renew=True)
tswp.sweepVar(klib.misc.convert_tqp_nqp(tswp.get_vals('psdphs_tqp')), 'rolloff_nqp', renew = True)
tswp.sweepVar(klib.misc.convert_tqp_nqp(tswp.get_vals('psdphs_tqp')) * V, 'rolloff_Nqp', renew = True)
tswp.print_params()

fitres__dfr_dnqp = tswp.fit_linear('Nqp','swp_fr',plotted=False,fitname='c+lin')
fitres__d1overQi_dnqp = tswp.fit_linear('Nqp','swp_Qi',y_invert=True,plotted=False,fitname='c+lin')

nepdata = klib.nep.calc_nep(kid.psd,
                            fitres__dfr_dnqp.params['a'].value,
                            fitres__d1overQi_dnqp.params['a'].value,
                            kid.swp.fitresult,
                            tqp=kid.psdfitresult.params['tqp'].value,
                            Nqp=klib.misc.convert_Kelvin_nqp(np.array(T)) * V)

```
