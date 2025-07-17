# anayzer_db

Observation libraries with python3.
This module includes:

- read rawfile with `rhea_comm`
- analysis rawdata with `mkid_pylibs`
- getting sweep data and resonant frequency from fitting and TOD data with 1kSPS for any duration from configuration file.
- getting sweep data and resonant frequency from fitting and TOD data with some sampling rate condition from configuration file.

## Set up
`rhea_comm` should be in `mkid_pylibs` directory.

```
git clone --recursive git@github.com:groundbird/analyzer_db.git
cd mkid_pylibs
ln -s ../rhea_comm .
cd ..
```

