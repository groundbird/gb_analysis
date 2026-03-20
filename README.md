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
