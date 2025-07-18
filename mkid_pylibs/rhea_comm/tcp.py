#!/usr/bin/env python3
'''TCP communication handling codes.'''
import socket
from time import perf_counter
from rhea_pkg import IP_ADDRESS_DEFAULT, TCP_PORT_DEFAULT

TCP_N_TRY = 10
TCP_BUFFSIZE = 2**18
TCP_TIMEOUT = 0.1


class TCPError(Exception):
    '''TCP error'''


class TCP:
    '''TCP communication handler.'''
    def __init__(self, ip_address=IP_ADDRESS_DEFAULT, port_num=TCP_PORT_DEFAULT):
        self.error_try_count = TCP_N_TRY
        self.buff_size = TCP_BUFFSIZE

        self.sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self.sock.connect((ip_address, port_num))
        self.sock.settimeout(TCP_TIMEOUT)
        self.read_buff = b''

    def __pull(self, length):
        error_cnt = 0
        tmplen = min(length, self.buff_size)
        while error_cnt < self.error_try_count:
            if len(self.read_buff) >= length:
                break
            error_cnt += 1
            try:
                self.read_buff += self.sock.recv(tmplen)
                error_cnt = 0
            except socket.timeout:
                pass

        if error_cnt >= self.error_try_count:
            raise TCPError('tcp.__pull: receive error')
        buff = self.read_buff[0:length]
        self.read_buff = self.read_buff[length:]
        return buff

    def clear(self):
        '''Clear the TCP buffer.'''
        try:
            while True:
                self.__pull(self.buff_size)
        except TCPError:
            pass
        self.read_buff = b''

    def read(self, length=1):
        '''Read TCP stream.

        Parameter
        ---------
        length : int
            Read length.

        Returns
        -------
        data : bytes
            Data bytes.
        '''
        try:
            buff = self.__pull(length)
        except TCPError as err:
            print(err)
            buff = self.read_buff
            self.read_buff = b''
        except TypeError as err:
            print(err)
            buff = self.read_buff
            self.read_buff = b''

        return buff

    def send(self, data):
        '''Write data.'''
        self.sock.send(data)


def main():
    '''Main function.'''
    client = TCP(ip_address='127.0.0.1', port_num=32409)

    test_length = 10_000_000
    lenbytes = bytes(f'{test_length}', encoding='utf8')

    client.send(lenbytes)
    start = perf_counter()
    data = client.read(test_length)
    stop = perf_counter()
    print(f'Length: {len(data)}.')
    print(f'Time: {stop - start} s.')
    print(f'Speed: {len(data)*8/(stop - start)/1e6} Mbps')


if __name__ == '__main__':
    main()
