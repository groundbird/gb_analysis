#!/usr/bin/env python3
# coding: utf-8
'''Contorll IDELAY and ODELAY for the ADC/DAC interfaces.
'''

from raw_setting import RawSetting

OFFSET_DAC = 0x2200_0000
OFFSET_ADC = 0x1200_0000


class DelayController(RawSetting):
    '''
    Class to handle delay configuration.
    '''
    def __init__(self, rbcp_ins, verbose=True, offset_dac=OFFSET_DAC, offset_adc=OFFSET_ADC):
        RawSetting.__init__(self, rbcp_ins, verbose, 'DELAY')
        self._offset_dac = offset_dac
        self._offset_adc = offset_adc

    def _dac_addr_lsb(self, lane):
        return self._offset_dac + 0x01_00 + lane

    def _dac_addr_msb(self, lane):
        return self._offset_dac + 0x02_00 + lane

    def _adc_addr_lsb(self, channel, lane):
        return self._offset_adc + 0x01_00 + ((2*channel) << 8) + lane

    def _adc_addr_msb(self, channel, lane):
        return self._offset_adc + 0x02_00 + ((2*channel) << 8) + lane

    def set_dac_delay(self, lane, value):
        '''
        Set delay value to the DAC ODELAYE3.

        Parameters
        ----------
        lane : int
            Lanes to be configured. Should be 0 to 7 for the data lanes.
            8 is assigned to the frame lane.
        value : int
            Delay count. Should be within 0 to 511
        '''
        assert 0 <= lane < 9
        assert 0 <= value < 512

        lsb = value & 0xff
        msb = (value - lsb) >> 8
        address_lsb = self._dac_addr_lsb(lane)
        address_msb = self._dac_addr_msb(lane)

        self._write(address_lsb, lsb)
        self._write(address_msb, msb)

    def get_dac_delay(self, lane):
        '''
        Get delay value of the DAC ODELAYE3.

        Parameters
        ----------
        lane : int
            Lanes to be configured. Should be 0 to 7 for the data lanes.
            8 is assigned to the frame lane.

        Returns
        -------
        delay : int
            Delay taps.
        '''

        assert 0 <= lane < 9
        lsb = self._read(self._dac_addr_lsb(lane))
        msb = self._read(self._dac_addr_msb(lane))

        return lsb + (msb << 8)

    def set_adc_delay(self, channel, lane, value):
        '''
        Set delay value to the ADC IDELAYE3.

        Parameters
        ----------
        ch : int
            0: ch A, 1: ch B
        lane : int
            Lanes to be configured. Should be 0 to 6.
        value : int
            Delay count. Should be within 0 to 511
        '''

        assert 0 <= lane < 7
        assert 0 <= value < 512

        lsb = value & 0xff
        msb = (value - lsb) >> 8
        address_lsb = self._adc_addr_lsb(channel, lane)
        address_msb = self._adc_addr_msb(channel, lane)

        self._write(address_lsb, lsb)
        self._write(address_msb, msb)

    def get_adc_delay(self, channel, lane):
        '''
        Get delay value of the ADC IDELAYE3.

        Parameters
        ----------
        ch : int
            0: ch A, 1: ch B
        lane : int
            Lanes to be configured. Should be 0 to 6.

        Returns
        -------
        delay : int
            Delay taps.
        '''

        assert 0 <= lane < 7
        lsb = self._read(self._adc_addr_lsb(channel, lane))
        msb = self._read(self._adc_addr_msb(channel, lane))

        return lsb + (msb << 8)

def main():
    '''Main function'''
    print('Nothing.')

if __name__ == '__main__':
    main()
