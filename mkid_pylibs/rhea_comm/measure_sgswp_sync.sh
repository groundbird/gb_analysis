#!/bin/bash

# environmental setting
## for dynabook
#rhea_comm=/mnt/c/Users/kucmb/ownCloud/ikemitsu/FPGA/src/rhea_comm
## for Monarch
rhea_comm=/mnt/d/Users/kucmb/takuji/ownCloud/ikemitsu/FPGA/src/rhea_comm
## for dodo
#rhea_comm=/home/taku/src/rhea_comm
## for mac
#rhea_comm=/Users/takuji/ownCloud/ikemitsu/FPGA/src/rhea_comm

measure_sgswp_path=$rhea_comm/measure_sgswp_sync.py
sg_manager_path=$rhea_comm/sg_manager.py
sg_sweep_path=$rhea_comm/sg_sweep_sync.py

# rhea_dac frequency setting
rhea_freq="20"
amps="1"

function measure_sgswp(){
    measure_arg=$rhea_freq" -a "$amps" -r "$rate_ksps" -fc "$f_center" -fw "$f_width_half" -fs "$f_step" -d "$dwell
    python3 $measure_sgswp_path $measure_arg & >> $tmp 2>/dev/null
    pid=$!
    sleep 2s
}

function monitor_sgswp(){
    measure_arg=$rhea_freq" -a "$amps" -r "$rate_ksps" -fc "$f_center" -fw "$f_width_half" -fs "$f_step" -d "$dwell" --mon "True" -f "sgswp_mon.rawdata
    python3 $measure_sgswp_path $measure_arg & >> $tmp 2>/dev/null
    pid=$!
    sleep 2s
}

# function measure_sgswp_fast(){
#     measure_arg=$rhea_freq" -a "$amps" -r "$rate_ksps" -fc "$f_center" -fw "$f_width_half" -fs "$f_step" -d "$dwell
#     python3 $measure_sgswp_path $measure_arg --fastswp & >> $tmp 2>/dev/null
#     pid=$!
#     sleep 1s
# }

function sg_normal_status(){
    python3 $sg_manager_path -f $f_center -p on --port $compath > /dev/null
    echo
}

function normal_swp(){
    python3 $sg_sweep_path normal -fc $f_center -fw $f_width_half -fs $f_step -d $dwell -r $run --port $compath
}

# function fast_swp(){
#     python3 $sg_sweep_path fast -fc $f_center -fw $f_width_half -p $points -d $dwell -r $run --port $compath
# }

function usage(){
    echo "Usage: $PROGNAME [OPTIONS]"
    echo "Options:"
    echo " -h, --help"
    echo " -r, --rate sample rate(kSPS)"
    echo "     default :" $rate_ksps
    echo " -fc, --fcenter f_center(with unit)"
    echo "     default :" $f_center
    echo " -fw, --fwidth f_width_half(with unit)"
    echo "     default :" $f_width_half
    echo " -fs, --fstep f_step(with unit)"
    echo "     default :" $f_step
    echo " -d, --dwell dwell time in us"
    echo "     default :" $dwell
    echo " -rt, --run # of run swp times"
    echo "     default :" $run
    echo " -c, --comport [com-port path]"
    echo "     default :" $compath
    echo " -m, --mon [ 0 or 1 ]"
    echo "     default :" $mon
}
#    echo " -f, --fastswp [no arg]:use fastswp mode"
#    echo " -p, --point [# of swp points]:for fastswp mode"
#    echo "     default :" $points
#}


PROGNAME=$(basename $0)

# default setting for measure rhea
rate_ksps=10 #kSps

# default setting for swp
mode=normal
f_center=4GHz
f_width_half=2MHz
f_step=1kHz
dwell=1000 # us
run=1
points=2000

# if use as monitor, mon=1
mon=0

# default com port path
compath=/dev/ttyS5

# output message
tmp=./tmp

for OPT in "$@"
do
    case "$OPT" in
        '-h'|'--help' )
	    usage
	    exit 1
	    ;;
        '-r'|'--rate' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
                echo "$PROGNAME: option requires an argument -- $1" 1>&2
                exit 1
	    fi
	    rate_ksps="$2"
	    shift 2
	    ;;
        '-fc'|'--fcenter' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
                echo "$PROGNAME: option requires an argument -- $1" 1>&2
                exit 1
	    fi
	    f_center="$2"
	    shift 2
	    ;;
        '-fw'|'--fwidth' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
                echo "$PROGNAME: option requires an argument -- $1" 1>&2
                exit 1
	    fi
	    f_width_half="$2"
	    shift 2
	    ;;
        '-fs'|'--fstep' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
                echo "$PROGNAME: option requires an argument -- $1" 1>&2
                exit 1
	    fi
	    f_step="$2"
	    shift 2
	    ;;
        '-d'|'--dwell' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
                echo "$PROGNAME: option requires an argument -- $1" 1>&2
                exit 1
	    fi
	    dwell="$2"
	    shift 2
	    ;;
        '-rt'|'--run' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
                echo "$PROGNAME: option requires an argument -- $1" 1>&2
                exit 1
	    fi
	    run="$2"
	    shift 2
	    ;;
        '-c'|'--comport' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
                echo "$PROGNAME: option requires an argument -- $1" 1>&2
                exit 1
	    fi
	    compath="$2"
	    shift 2
	    ;;
        '-m'|'--mon' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
                echo "$PROGNAME: option requires an argument -- $1" 1>&2
                exit 1
	    fi
	    mon="$2"
	    shift 2
	    ;;
#        '-f'|'--fastswp' )
#	    $mode=fast
#	    shift 1
#	    ;;
        -*)
	    echo "$PROGNAME: illegal option -- '$(echo $1 | sed 's/^-*//')'" 1>&2
	    exit 1
	    ;;
    esac
done

sg_normal_status
touch $tmp

# if test $mode = normal; then
#     measure_sgswp
#     normal_swp
# elif test $mode = fast; then
#     measure_sgswp_fast
#     fast_swp
# fi

if test $mon = 1; then
    monitor_sgswp
elif test $mon = 0; then
    measure_sgswp
fi

run=$((run+1))
normal_swp

cat $tmp

while ps -p $pid > /dev/null
do
    sleep 1s
    tail -n1 $tmp
done

sg_normal_status

tail -n1 $tmp
\rm $tmp
