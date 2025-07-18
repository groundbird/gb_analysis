import sys
from read_raw_cpp import read_rawdata_cpp as rcpp

def process_meas_id(meas_id):
    meas_id = int(meas_id)
    rcpp(meas_id=meas_id, log=True, saveraw=True)


def main(meas_ids):
    for meas_id in meas_ids:
        process_meas_id(meas_id)

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python process_rawdata.py <meas_id1> [meas_id2] [meas_id3] ...")
        sys.exit(1)

    meas_ids = sys.argv[1:]
    main(meas_ids)
