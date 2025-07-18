#!/usr/bin/env python3
'''Class to manage `clock_man` module in the firmware.
   Clocking wizard can be reconfigured using the RBCP commands shown here.
'''

from raw_setting import RawSetting

OFFSET_CLOCK_MAN = 0x1300_0000

# Register address
CLOCK_MAN_SSR = 0x00 # Software reset register
CLOCK_MAN_SR = 0x04 # Status Register
CLOCK_MAN_ER = 0x08 # Clock monitor error status register
CLOCK_MAN_IS = 0x0C # Intterupt Status
CLOCK_MAN_IE = 0x10 # Intterupt Enable

CLOCK_MAN_CR0 = 0x200
CLOCK_MAN_CR1 = 0x204
# CLKOUT0 - 6
CLOCK_MAN_CRDIV = lambda _ch: 0x208 + 3*4*_ch
CLOCK_MAN_CRPHA = lambda _ch: 0x20C + 3*4*_ch
CLOCK_MAN_CRDUT = lambda _ch: 0x210 + 3*4*_ch
CLOCK_MAN_NCH = 7
# bit0: load/sen, bit1: saddr
CLOCK_MAN_CR23 = 0x25C

# Data masks
CLOCK_MAN_LOCKEDM = 0x00_00_00_01
CLOCK_MAN_SSR_KEY = 0x00_00_00_0A


class ClockMan(RawSetting):
    '''Class to manage `clock_man` module.
       `clock_man` module provides the way to reconfigure Xilinx's
       clocking wizard via RBCP-AXI interconnect module.
    '''
    def __init__(self, rbcp_inst, verbose=True, offset=OFFSET_CLOCK_MAN):
        RawSetting.__init__(self, rbcp_inst, verbose, label='clock_man')
        self._offset = offset


    def _write_axi(self, addr, data):
        self._write4(self._offset + addr, data, byteorder='little')

    def _read_axi(self, addr):
        return self._read4(self._offset + addr, byteorder='little')

    @property
    def _sr(self):
        return self._read_axi(CLOCK_MAN_SR)

    @property
    def locked(self):
        '''Clock status

        Returns
        -------
        locked : bool
            True if locked.
        '''
        return (self._sr & CLOCK_MAN_LOCKEDM) == CLOCK_MAN_LOCKEDM

    def soft_reset(self):
        '''Software reset'''
        self._write_axi(CLOCK_MAN_SSR, CLOCK_MAN_SSR_KEY)

    def phase_reconf(self, channel, angle):
        '''Reconfigure clock phase.

        Parameter
        ---------
        angle : float
            Angle of phase in degree.
        '''
        reg_addrs = [0x200 + 4*i for i in range(2 + 3*CLOCK_MAN_NCH)]
        reg_dict = {addr : self._read_axi(addr) for addr in reg_addrs}
        self.soft_reset()

        reg_dict[CLOCK_MAN_CRPHA(channel)] = int(angle * 1000)

        for addr, val in reg_dict.items():
            self._write_axi(addr, val)

        if self.locked:
            self._write_axi(CLOCK_MAN_CR23, 7)
            self._write_axi(CLOCK_MAN_CR23, 2)
        else:
            raise Exception('Clocking wizard not working.')


def main():
    '''Main function'''
    print('Nothing to do')


if __name__ == '__main__':
    main()
