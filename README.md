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
1. Edit measurement IDs in `process_measurements.sh`:
```bash
MEAS_IDS=(
    8448
    8449
)
```

2. Execute:
```bash
chmod +x process_measurements.sh
./process_measurements.sh
```

### Command Line
```bash
python use_raw_cpp.py 8448
```

## Output Data

### Save Location
- Directory: `./raw_data/`
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
