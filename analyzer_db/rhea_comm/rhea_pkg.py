#!/usr/bin/env python3
'''RHEA firmware information'''

FREQ_CLK_HZ = 200_000_000
RATE_KSPS_DEFAULT   = 1 # 1 kSPS ~ downsample:200000

PHASE_WIDTH      = 32
SIN_COS_WIDTH    = 16
ADC_DATA_WIDTH   = 14
IQ_DATA_WIDTH    = SIN_COS_WIDTH + ADC_DATA_WIDTH + 1
IQ_DS_DATA_WIDTH = 56
DS_RATE_MIN      = 20
DS_RATE_MAX      = 200000
DDS_AMP_BW       = 8

IP_ADDRESS_DEFAULT = '192.168.10.16'
RBCP_PORT_DEFAULT = 4660
TCP_PORT_DEFAULT = 24

def main():
    '''Nothing to do'''
    print('Hi.')

if __name__ == '__main__':
    main()
