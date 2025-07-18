#!/usr/bin/env python3
'''Contains definition of a class to manage countup_man module.
   `countup_man` module in the firmware enables debugging ADCs.
   This class provides functionalities to operate the module.
'''

from time import sleep

from raw_setting import RawSetting
from rhea_pkg import ADC_DATA_WIDTH, FREQ_CLK_HZ

# Register address
CMAN_STATUS = 0
CMAN_COUNTER_LSB = 1
CMAN_COUNTER_MSB = 2
CMAN_WRONG_DATA_LSB = 3
CMAN_WRONG_DATA_MSB = 4

# Data masks
CMAN_IRQ_MASK = 0b0000_0010
CMAN_SRST = 0b0000_0001


class CountupMan(RawSetting):
    '''Class to manage countup_man module.
       `countup_man` module in the firmware enables debugging ADCs.
       This class provides functionalities to operate the module.
    '''
    def __init__(self, rbcp_inst, verbose=True, label='cman0',
                 offset=0x1400_0000):
        RawSetting.__init__(self, rbcp_inst, verbose, label)
        self._offset = offset

    def _write_reg(self, addr, data):
        self._write(self._offset + addr, data)

    def _read_reg(self, addr):
        return self._read(self._offset + addr)

    def soft_reset(self):
        '''Perform software reset.'''
        self._write_reg(CMAN_STATUS, CMAN_SRST)

    def get_irq(self):
        '''Check whether interrupt request has been publised or not.
        If there are discrepancy between the internal counter
        and the data after the software reset, irq becomes high.

        Returns
        -------
        irq : bool
            True iff irq has been published.
        '''
        return (self._read_reg(CMAN_STATUS) & CMAN_IRQ_MASK) != 0

    @property
    def irq(self):
        '''Alias for `irq` property

        Returns
        -------
        irq : bool
            True iff irq has been published.
        '''
        return self.get_irq()

    def get_counter(self):
        '''Get `counter` value.

        Returns
        -------
        counter : int
            Counter value in the `countup_man` core
        '''
        lsb = self._read_reg(CMAN_COUNTER_LSB)
        msb = self._read_reg(CMAN_COUNTER_MSB)

        return (msb << 8) + lsb

    def get_wrong_data(self):
        '''Get `wrong_data` value.
        Returns
        -------
        wrong_data : int
            Value of `wrong_data` in the `countup_man` core
        '''

        lsb = self._read_reg(CMAN_WRONG_DATA_LSB)
        msb = self._read_reg(CMAN_WRONG_DATA_MSB)

        return (msb << 8) + lsb

    def line_diff(self):
        '''Bitwize difference between `counter` and `wrong_data`.
        Meaningful only if irq is high.

        Returns
        -------
        line_diff : int
            `counter` xor `wrong_data`
        '''
        return self.get_counter() ^ self.get_wrong_data()

    def health_check(self):
        '''Check whether irq happens during one cycle of ramp.

        Returns
        -------
        health : bool
            True if healthy.
        '''
        duration = 2**(ADC_DATA_WIDTH)/FREQ_CLK_HZ
        self.soft_reset()
        sleep(duration)

        return not self.irq


def main():
    '''Main function'''
    print('Nothing to do')


if __name__ == '__main__':
    main()
