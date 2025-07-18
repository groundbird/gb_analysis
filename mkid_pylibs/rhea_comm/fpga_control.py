#!/usr/bin/env python3
'''Readout controller.'''
import platform
import sys
import os
import stat

from rbcp          import RBCP
from info_setting  import InfoSetting
from adc_setting   import AdcSetting
from dac_setting   import DacSetting
from dds_setting   import DdsSetting
from iq_setting    import IQSetting
from ds_setting    import DsSetting
from snap_setting  import SnapSetting
from trg_setting   import TrgSetting
from debug_setting import DebugSetting

from countup_man import CountupMan
from clock_man import ClockMan

from tcp           import TCP


from rhea_pkg import IP_ADDRESS_DEFAULT, TCP_PORT_DEFAULT, RBCP_PORT_DEFAULT

IS_WINDOWS = platform.system() == 'Windows'


class FPGAControl:
    '''FPGA integrated controller.

    Parameters
    ----------
    verbose : bool
        Verbosity.
    ip_address : str
        IP address.
    '''
    def __init__(self, verbose=False, ip_address=IP_ADDRESS_DEFAULT):
        self._verbose = verbose
        self.__lock_path = '/tmp/.'+ip_address+'.lock'
        self._vprint(f'lock file: {self.__lock_path}')

        ## lock file
        try:
            if not IS_WINDOWS:
                import fcntl
                if not os.path.isfile(self.__lock_path):
                    self.__lockf = open(self.__lock_path, 'a', encoding='utf-8')
                    self.__lockf.close()
                    os.chmod(self.__lock_path,
                            mode=stat.S_IRUSR |\
                                stat.S_IWUSR |\
                                stat.S_IRGRP |\
                                stat.S_IWGRP |\
                                stat.S_IROTH |\
                                stat.S_IWOTH) # chmod 666

                self.__lockf = open(self.__lock_path, 'a', encoding='utf-8')
                fcntl.lockf(self.__lockf.fileno(), fcntl.LOCK_EX | fcntl.LOCK_NB)
        except IOError as err:
            print(f"{self.__lock_path} is locked.")
            print(err)
            sys.exit(1)


        self.rbcp = RBCP(ip_address=ip_address, port_num=RBCP_PORT_DEFAULT)
        self.tcp = TCP(ip_address=ip_address, port_num=TCP_PORT_DEFAULT)

        self.info = InfoSetting(self.rbcp, verbose=verbose)
        self.max_ch = self.info.max_ch
        self.en_snap = self.info.en_snap
        self.trig_ch = self.info.trig_ch

        self.adc_setting = AdcSetting(self.rbcp, verbose=verbose)
        self.dac_setting = DacSetting(self.rbcp, verbose=verbose)
        self.dds_setting = DdsSetting(self.rbcp, self.max_ch, verbose=verbose)
        self.iq_setting = IQSetting(self.rbcp, verbose=verbose)
        self.ds_setting = DsSetting(self.rbcp, verbose=verbose)

        self.count0 = CountupMan(self.rbcp, verbose=verbose, offset=0x1400_0000)
        self.count1 = CountupMan(self.rbcp, verbose=verbose, offset=0x1400_0100)

        self.clock = ClockMan(self.rbcp, verbose=verbose)

        if self.en_snap:
            self.snap_setting  = SnapSetting(self.rbcp, verbose=verbose)
        if self.trig_ch > 0:
            self.trg_setting = TrgSetting(self.rbcp, self.trig_ch, verbose=verbose)

        self.debug = DebugSetting(self.rbcp, verbose=verbose)


    def _vprint(self, *args, **kwargs):
        if self._verbose:
            if 'file' in kwargs:
                print(*args, **kwargs)
            else:
                print(*args, **kwargs, file=sys.stderr)


    def init(self):
        '''Initialize FPGA'''
        self._vprint(f'Firmware version: {self.info.version:d}')

        self.adc_setting.reset()
        self.dac_setting.reset()

        self.dac_setting.channel_swap()
        # self.adc.channel_swap()
        self.dac_setting.txenable_on()

        self.ds_setting.set_accum(200000)

        self.dds_setting.reset()

        self.iq_setting.set_read_width(0)
        self.iq_setting.iq_off()
        self.iq_setting.time_reset()
        self.iq_setting.clear_fifo_error()

        if self.en_snap:
            self.snap_setting.snap_off()
            self.snap_setting.time_reset()

        if self.trig_ch > 0:
            self.trg_setting.init()
