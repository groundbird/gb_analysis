#!/usr/bin/env python3
'''Common functions used in rhea_comm'''

def two_div(val):
    '''Rounde-up value of log2(val).'''
    if val < 2:
        return 0

    return len(bin(val-1)) -2


def packet_size(read_width):
    '''Packet length in bytes.
    Header, footer (1 byte each)
    Timestamp (5 bytes)
    Data (7 bytes * 2 (I and Q) * read_width)

    Parameter
    ---------
    read_length : int
        Number of readout channels to be read.

    Returns
    -------
    packet_size : int
        Length of packets in bytes.
    '''
    return 2 + 5 + 7*2*read_width
