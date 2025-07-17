# mkid_pylibs
MKID analysis libraries with python3

This module includes:
* read rawfile with `rhea_comm` and register them as sweep/tod data classes
* peak finding in sweep data
* sweep fitting
* rewinding (removing baseline offset) referred as `rdata`
* two-tone subtraction referred as `mdata`
* calculating PSD from tod
* PSD curve fitting to extract MKID properties
* analyzing results with temperture variations
* analyzing triggered tod by fitting the shape
* NEP from PSD and given dfr/dNqp and d(1/Qr)/dNqp 

## Setup

`rhea_comm` should be in `mkid_pylibs` directories to treat rhea files.

```shell
git clone git@github.com:groundbird/mkid_pylibs.git
cd mkid_pylibs
ln -s [DIRECTORY]/rhea_comm
cd ../
```

## How to use

### Data with `PsdAnalyzer` class
```python
import mkid_pylibs as klib

lo = '4.99GHz'
swpset = klib.readfile_swp('rhea','SWP_FILE_NAME',-1,lo) # All tone data will be read with `-1` otherwise set integer >=0 to select the index number to be read.
todset = klib.multreadfile_tod('rhea',['TOD_FILE_1KSPS','TOD_FILE_100kSPS', 'TOD_FILE_1000kSPS', ...],lo) # Register list of TOD files.

klist = [
    klib.kidana_psd.PsdAnalyzer(swp=swpset[0], tod=todset[0], ctod=todset[1]),
    klib.kidana_psd.PsdAnalyzer(swp=swpset[2], tod=todset[2], ctod=todset[3]),
    klib.kidana_psd.PsdAnalyzer(swp=swpset[4], tod=todset[4], ctod=todset[5]),
]

for ik in klist:
    ik.fitIQ(nfwhm=3)
    ik.calcpsd(dofit=False)
    klib.plotter.plotSwpTOD(ik.swp,ik.todlist)
    klib.plotter.plotPSD(ik.psd)
```

### Data with `KidAnalyzer` class
```python
import mkid_pylibs as klib

lo = '4.99GHz'
swpset = klib.readfile_swp('rhea','SWP_FILE_NAME',-1,lo) # All tone data will be read with `-1` otherwise set integer >=0 to select the index number to be read.
todset = klib.readfile_tod('rhea','TOD_FILE_NAME',-1,lo) # All tone data will be read with `-1` otherwise set integer >=0 to select the index number to be read.

klist = [
    klib.kidana_psd.PsdAnalyzer(swp=swpset[0], tod=todset[0], ctod=todset[1]),
    klib.kidana_psd.PsdAnalyzer(swp=swpset[2], tod=todset[2], ctod=todset[3]),
    klib.kidana_psd.PsdAnalyzer(swp=swpset[4], tod=todset[4], ctod=todset[5]),
]

fig,ax = klib.misc.make_figure(fsz=(20,6))
for i,ik in enumerate(klist):
    ik.fitIQ(nfwhm=3)
    ax.plot(ik.tod.time, ik.tod.rwmdata.corphase - np.mean(ik.tod.rwmdata.corphase[:100]), label=f'KID{i}', marker='', ls='-')

ax.set_xlabel('time [sec]')
ax.set_ylabel('phase [rad]')
klib.misc.adj_figure(fig)
plt.show()
```

### How to use classes

#### freq
- `ik.tod.f`

#### rate
- `ik.tod.rate`

#### time
- `ik.tod.time`

#### phase
- `ik.tod.data.phase` # raw data
- `ik.tod.rwdata.phase` # rewinded data
- `ik.tod.mdata.phase` # modified data (on/off subtractions for removing common noise)
- `ik.tod.rwmdata.phase` # rewined modified data

##### Corrected phase

phase is better to be re-defined to avoid the jump btw -pi and pi.

--> Use corphase instead of phase in rewinded data.

- `ik.tod.rwdata.corphase` # rewinded data + redefined phase
- `ik.tod.rwmdata.corphase` # rewinded modified data + redefined phase <-- This is nominal data to be used in the analysis

##### Linearization
If large signals are detected in KIDs, they could be non-linear responces.

This could be corrected by the linearization function

- `ik.tod.rwmdata.lincorphase` # rewinded modified data + redefined phase + linearization


#### amplitude
- `ik.tod.data.amplitude` # raw data
- `ik.tod.rwdata.amplitude` # rewinded data
- `ik.tod.mdata.amplitude` # modified data (on/off subtractions for removing common noise)
- `ik.tod.rwmdata.amplitude` # rewined modified data

##### Normalization

After rewinding, the radius of amplitude is 0.5.

It might be better to define normalized amplitude (such as for the PSD calculation)

- `ik.tod.rwdata.coramplitude` # rewinded data + normalized
- `ik.tod.rwmdata.coramplitude` # rewinded data + normalized

### Detailed information

See here: [UNDER CONSTRUCTION]

## Short description for each file

* `__init__.py`

### Main class to be used
* `kidana.py`: integrated class to analyze KIDs for observations
* `kidana_psd.py`: integrated class to analyze KIDs for PSD evaluations
* `sweep_meas.py`: analysis with some parameter variations (such as temperature)

### Data class
* `drawer.py`: basic class `DrawData` for supporting plot function
* `base.py`: basic class `BaseData` for integration of two arrays into one class such as `SwpData`, `TodData`, and `PsdData`.
* `generalarray.py`: function to define the custom named_array
* `freqarray.py`: frequency data class
* `timearray.py`: time data class
* `kiddata.py`: KID IQ data class
* `swpdata.py`: frequency + KID IQ class for sweep
* `toddata.py`: time + KID IQ class for TOD
* `psddata.py`: frequency + KID amp/phs named_array for PSD

### IQ fitting
* `kidfit.py`: IQ fitting interface
* `fitters.py`: IQ fitting functions
* `peak_search.py`: simple peak search tool with Lorenzian

### File treatment
* `readfile.py`: getting one KID data from an input file with `rhea_comm`
* `multreadfile.py`: getting all KID data from an input file with `rhea_comm`
* `parsefile.py`: how to treat each file type

### Analysis
* `modIQ.py`: IQ modifications (two-tone subtraction and rewinding)
* `psd.py`: PSD calculation/analysis
* `nep.py`: NEP calculation
* `trig.py`: trigger analysis

### Others
* `plotter.py`: example plotting functions
* `misc.py`: miscellaneous useful functions
* `fit/` : general fitting interface using lmfit

