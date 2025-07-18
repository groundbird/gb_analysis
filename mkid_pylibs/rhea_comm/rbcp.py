#!/usr/bin/env python3
'''RBCP communication.
'''
from sys import stderr
from time import time, sleep
import socket

from rhea_pkg import IP_ADDRESS_DEFAULT, RBCP_PORT_DEFAULT

RBCP_VER_TYPE   = 0xff
RBCP_CMD_FLAG_R = 0xc0
RBCP_CMD_FLAG_W = 0x80
RBCP_FLAG_MASK  = 0x0f
RBCP_CMD_MASK = 0xf0
RBCP_FLAG_CHECK = 0x8
RBCP_ERROR_MASK = 0x1

RBCP_HEADER_LENGTH = 8

class RBCPPacket:
    '''RBCP packet class'''
    def __init__(self, is_read, packet_id, length, address, data=b'', flag=0):
        self.is_read = is_read
        self.packet_id = packet_id
        self.length = length
        self.address = address
        self.data = data
        self.flag = flag

        if self.is_write and (self.length != len(self.data)):
            raise RBCPError('Inconsistency in length.')

    @property
    def is_write(self):
        '''True if the packet is in write mode.'''
        return not self.is_read

    @property
    def packet_length(self):
        '''Calculate packet length.

        Returns
        -------
        packet_length : int
            Length of RBCP packet.
        '''
        return 8 + (0 if self.is_read else self.length)

    def repr(self):
        '''Byte representation of the packet

        Returns
        -------
        repr_byte : bytearray
            Representation of the packet.
        '''
        repr_byte = bytearray(self.packet_length)
        repr_byte[0] = RBCP_VER_TYPE
        repr_byte[1] = RBCP_CMD_FLAG_R if self.is_read else RBCP_CMD_FLAG_W
        repr_byte[1] += self.flag
        repr_byte[2] = self.packet_id
        repr_byte[3] = self.length
        repr_byte[4:8] = self.address.to_bytes(4, byteorder='big')
        repr_byte[8:] = self.data

        return repr_byte

    @classmethod
    def interpret(cls, packet):
        '''Interpret the packet and return RBCPPacket object.

        Parameter
        ---------
        packet : bytes or bytearray
            RBCP packet.
        '''
        repr_byte = bytearray(packet)
        assert repr_byte[0] == RBCP_VER_TYPE

        cmd = repr_byte[1] & RBCP_CMD_MASK
        if cmd == RBCP_CMD_FLAG_R:
            is_read = True
        elif cmd == RBCP_CMD_FLAG_W:
            is_read = False
        else:
            raise RBCPError('Wrong cmd.')

        flag = repr_byte[1] & RBCP_FLAG_MASK

        packet_id = repr_byte[2]
        length = repr_byte[3]
        address = int.from_bytes(repr_byte[4:8], byteorder='big')
        data = repr_byte[8:]

        return cls(is_read, packet_id, length, address, data, flag)

    @property
    def is_ack(self):
        '''Acknowledge bit status

        Returns
        -------
        is_ack : bool
            True when ack bit is high.
        '''
        return bool(self.flag & RBCP_FLAG_CHECK)

    @property
    def bus_error(self):
        '''Bus error status.

        Returns
        -------
        bus_error : bool
            True when bus error bit is high.
        '''
        return bool(self.flag & RBCP_ERROR_MASK)

    def __str__(self) -> str:
        srepr =  f'is_read  : {self.is_read}\n'
        srepr += f'flag     : {self.flag:02x}\n'
        srepr += f'packet_id: {self.packet_id:02x}\n'
        srepr += f'length   : {self.length}\n'
        srepr += f'Address  : 0x{self.address:08x}'

        for i, d_int in enumerate(self.data):
            srepr += f'\nData{i:03d}: {d_int:02x}'

        return srepr

class RBCPError(Exception):
    '''RBCP error.'''

class RBCP:
    '''Configuration of the readout with protocol called RBCP.'''
    def __init__(self, ip_address=IP_ADDRESS_DEFAULT, port_num=RBCP_PORT_DEFAULT,
                 retry_max=30):
        # setting
        self._error_try_count = 10
        self._buff_size       = 4096
        self._retry_max       = retry_max

        self.packet_id  = int(time()) %  256

        self.read_buff  = b''
        # make socket
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        self.sock.connect((ip_address, port_num))
        self.sock.settimeout(1)

    def __push(self, buff):
        error_cnt = 0

        while error_cnt < self._error_try_count:
            error_cnt += 1
            send_byte = self.sock.send(buff)
            buff = buff[send_byte:]

            if len(buff) == 0:
                break

            sleep(0.1)

        if error_cnt >= self._error_try_count:
            raise RBCPError('rbcp.__push: send error')

    def __pull(self, length):
        error_cnt = 0
        while error_cnt < self._error_try_count:
            if len(self.read_buff) >= length:
                break
            error_cnt += 1

            try:
                self.read_buff += self.sock.recv(self._buff_size)
            except socket.timeout:
                pass

        if error_cnt >= self._error_try_count:
            raise RBCPError('rbcp.__pull: receive error')

        buff = self.read_buff[0:length]
        self.read_buff = self.read_buff[length:]

        return buff

    def clear(self):
        '''Flush UDP'''
        error_cnt = 0
        while error_cnt < self._error_try_count:
            try:
                self.sock.recv(self._buff_size)
            except socket.timeout:
                error_cnt += 1

    def _send(self, address, data=1):
        '''Create and send RBCP packet to the firmware.

        Parameters
        ----------
        address : int
            Register address. Should be less than 2**32
        data : bytes, bytearray or int
            To write data to register, specify data here with bytes or bytearray.
            To read register, put length here in integer.
        '''
        if isinstance(data, int):
            is_read = True
            payload = b''
            length = data
        else:
            is_read = False
            payload = data
            length = len(payload)

        self.packet_id = (self.packet_id + 1) % 256

        packet = RBCPPacket(is_read, self.packet_id, length, address,
                            data=payload)

        self.__push(packet.repr())

    def _recv(self):
        '''Receive data according to the packet information.'''
        while True:
            header = self.__pull(RBCP_HEADER_LENGTH)
            length = header[3]
            payload = self.__pull(length)

            packet = RBCPPacket.interpret(header + payload)

            if packet.packet_id == self.packet_id:
                break

        return packet.data

    def read(self, address, length=1):
        '''RBCP read.

        Parameters
        ----------
        address : int
            RBCP address. Should be less than 2**32.
        length : int
            Number of bytes to be read.

        Returns
        -------
        payload : bytearray
            Data.
        '''
        retry_cnt = 0
        while True:
            try:
                self._send(address, length)
                payload = self._recv()
                break
            except RBCPError as err:
                retry_cnt += 1
                if retry_cnt > self._retry_max:
                    raise RBCPError from err
                print(err)
                print(f"rbcp.write: error! retry: {retry_cnt}", file=stderr)
                sleep(0.1)
                continue

        return payload

    def write(self, address, data):
        '''Write data to RBCP register.

        Parameters
        ----------
        address : int
            RBCP register address.
        data : bytes or bytearray
            Data.

        Returns
        -------
        payload : bytearray
            Data.
        '''
        retry_cnt = 0
        while True:
            try:
                self._send(address, data)
                ret = self._recv()
                break
            except RBCPError as err:
                retry_cnt += 1
                if retry_cnt > self._retry_max:
                    raise RBCPError from err
                print(err)
                print(f"rbcp.write: error! retry: {retry_cnt}", file=stderr)
                sleep(0.1)
                continue

        return ret

    def read_intn(self, address, length):
        '''Read n bytes and return data interpreted as integer.

        Parameters
        ----------
        address : int
            RBCP address.
        length : int
            Length to be read.

        Returns
        -------
        data : int
            Data interpreted as a big-endian integer.
        '''
        buff = self.read(address, length=length)
        return int.from_bytes(buff, byteorder='big')

    def read_int1(self, address):
        '''Read 1 byte.'''
        return self.read_intn(address, 1)

    def read_int2(self, address):
        '''Read 2 bytes.'''
        return self.read_intn(address, 2)

    def read_int4(self, address):
        '''Read 4 bytes.'''
        return self.read_intn(address, 4)

    def write_intn(self, address, data_int, length):
        '''Write n bytes by intepreting given integer
        into a byte array of a specified length.

        Parameters
        ----------
        address : int
            RBCP address.
        data : int
            Data integer.
        length : int
            Interpretation length.

        Returns
        -------
        ret_data : bytearray
            Returned byte array.
        '''
        assert isinstance(data_int, int)

        data = data_int.to_bytes(length, byteorder='big')

        return self.write(address, data)

    def write_int1(self, address, data_int):
        '''Write 1 byte from int.'''
        return self.write_intn(address, data_int, 1)

    def write_int2(self, address, data_int):
        '''Write 2 bytes from int.'''
        return self.write_intn(address, data_int, 2)

    def write_int4(self, address, data_int):
        '''Write 4 bytes from int.'''
        return self.write_intn(address, data_int, 4)


def main():
    '''Main function.'''
    print('Nothing')

if __name__ == '__main__':
    main()
