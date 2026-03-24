# gb_analysis
# MKID Data Processing System
Based on Yoshinori's gb_cal implementation.

## Overview
Automatically calculates phase data from I-Q data acquired from multiple detector chips (1A, 1B, 2A, 2B, 3A, 3B, 220).

## Required Environment
### Required Libraries
- `gbirdproc` - C++ backend
- `mkid_pylibs` - KID analysis library
- `rhea_comm` - RHEA communication

## Usage

### Single Measurement Processing
```python
from read_rawdata_cpp import read_rawdata_cpp
# Process measurement ID 8448 (log enabled, save raw data)
data = read_rawdata_cpp(8448, log=True, saveraw=True)
```

### Batch Processing
1. Edit measurement IDs in `run_use_raw_cpp.sh`:
```bash
MEAS_IDS=(
    8448
    8449
)
```

2. Execute:
```bash
chmod +x process_measurements.sh
./run_use_raw_cpp.sh
```

### Command Line
```bash
python use_raw_cpp.py 8448
```

## Fit Configuration
### Manual Fitconfig (default)
Chip-specific fitting parameters are defined in `config.py` → `get_fitconfig(chip)`.  
Each chip has manually tuned settings for:
- Frequency range (`rangeind`, `freqranges`)
- Double resonance detection (`twokidind`, `twokidfitter`)
- Initial parameters (`initind`, `fitinit`)
- Fitting depth (`depind`, `dep`)
- Guess skip (`skipind`, `guessskip`)

### Auto Fitconfig (`auto_fitconfig=True`)
Defined in `auto_fitconfig.py`. Fitting parameters are automatically determined from sweep data without any manual configuration.

#### Double Dip Detection
Dips in the sweep amplitude are detected using `scipy.signal.find_peaks` on the smoothed (uniform filter, size=20) inverted amplitude.  
A KID is classified as a double resonance if exactly **2 peaks** are found above the prominence threshold (default: 5% of the amplitude range).

- `gaolinbg2f`: TOD tone is closer to **fr1** (the lower frequency dip)
- `gaolinbg2l`: TOD tone is closer to **fr2** (the higher frequency dip)

#### Fitting Range
Auto fitconfig does **not** set a frequency range (`rangeind`/`freqranges` are empty).  
The full sweep range is used for all KIDs (`frqrange=[None, None]`).  
If fitting fails due to a wide sweep range, switch to manual fitconfig and set `rangeind`/`freqranges` explicitly in `config.py`.

#### Initial Parameters
For each KID, initial parameters with constraints are automatically set as `lmfit.Parameter` objects:

| Parameter | Initial Value | Range |
|-----------|--------------|-------|
| `fr` / `fr1` / `fr2` | Detected dip frequency | ±2×10⁵ Hz |
| `Qr` / `Qr1` / `Qr2` | 1×10⁴ | 1×10³ – 1×10⁵ |
| `Qc` / `Qc1` / `Qc2` | 2×10⁴ | 1×10³ – 1×10⁵ |
| `phi0` / `phi01` / `phi02` | 0.0 | −π – π |
| `absa` | max(sweep amplitude) | — |

> **Note**: `dep`, `guess_skip` are not auto-detected and fall back to defaults (`dep=3`, `guess_skip=False`).

#### Limitations
- Detection may fail for noisy or asymmetric sweep profiles. In such cases, adjust `prominence_threshold` in `_has_double_dip()` or switch to manual fitconfig.
- Frequency range is not constrained, so nearby resonances from other KIDs may interfere with the fit.


## Output Data

### Save Location
- Directory: `./raw_data/`(config.py `SAVEDIR`)
- Filename: `{chip}_{meas_id}.pkl`

### Data Structure
```python
# Dictionary format
{
    'swp_param': DataFrame,     # Resonance parameters (fr, Qr, Qc, Qi)
    'utime': array,            # Time data
    'el': array,               # Elevation
    'az': array,               # Azimuth
    'phase': {                 # Phase data for each KID
        'kid00': array,
        'kid01': array,
        ...
    }
}
```

## Plot Output

### Save Location
Images are saved under `SAVEDIR/plot/` by category (configured by `SAVEDIR` in config.py).
```
SAVEDIR/
├── raw_data/                        # pickle data
│   └── {chip}_{meas_id}.pkl
└── plot/
    ├── swpamp/                      # Sweep amplitude
    │   └── swpamp_{meas_id}_{daq}_{chip}.jpg
    ├── bswpamp/                     # Blind tone sweep amplitude
    │   └── bswpamp_{meas_id}_{daq}_{chip}.jpg
    ├── swpiq/                       # Sweep I-Q circle
    │   └── swpiq_{meas_id}_{daq}_{chip}.jpg
    ├── psd/                         # Power spectral density
    │   └── psd_{meas_id}_{daq}_{chip}.jpg
    ├── tod/                         # Time ordered data
    │   └── tod_{meas_id}_{daq}_{chip}.jpg
    └── log/                         # Environmental log (temperature, PWV, humidity)
        └── log_{meas_id}_{daq}_{chip}.jpg
```

### Plot Methods
```python
data = read_rawdata_cpp(8448, log=True, saveraw=True)

data.plot_swpamp(save=True)   # Sweep amplitude for each KID
data.plot_bswpamp(save=True)  # Blind tone sweep amplitude
data.plot_swpiq(save=True)    # I-Q circle with fit result
data.plot_psd(save=True)      # Power spectral density
data.plot_tod(save=True)      # Time ordered data (KID vs blind tone)
data.plot_log(save=True)      # Environmental log (requires log=True)
```

> **Note**: `plot_log` is only available when initialized with `log=True` (default). If an error occurs, it will be skipped automatically.

## Main Features
- **Automatic fitting**: Resonance fitting with chip-specific configurations
- **Quality control**: Glitch detection and environmental condition filtering
- **Data output**: Save in pickle format
- **Batch processing**: Process multiple measurements at once

## Chip-DAQ Mapping
```
3A←→GB01, 2A←→GB02, 3B←→GB03, 1A←→GB04
2B←→GB05, 1B←→GB06, 220←→GB07
```

## TODO (Contributions Welcome!)
1. **Blind tone quality verification**
2. **Improve sweep fitting accuracy**: Currently uses initial values in fittingconfig, but it's sensitive to variations and some detectors are unusable
3. **Add data flags**
4. **Review data storage methods**
