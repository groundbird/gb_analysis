#!/usr/bin/env python3
'''Communication health check / reconfiguration codes.'''
from argparse import ArgumentParser
import sys
import numpy as np

from fpga_control import FPGAControl
from adc_setting import ADSChannel, ADSTestPattern
from rhea_pkg import IP_ADDRESS_DEFAULT, ADC_DATA_WIDTH


def dac_check(fpga:FPGAControl):
    '''Test DAC communication.

    Returns
    -------
    cert : dict
        Dictionary of diagnostic results.
    '''
    results = []

    for i in range(256):
        a_in = i
        b_in = i

        fpga.dac_setting.dac_test_input(a_in, b_in, chmode=True)
        fpga.dac_setting.set_test_pattern(a_in, b_in)
        fpga.dac_setting.clear_test_result()
        result = fpga.dac_setting.get_test_result()

        results.append(result)

    fpga.dac_setting.dac_normal_input()

    result = np.bitwise_or.reduce(results)
    cert = {'result': result == 0,
            'wrong_bits': result}

    return cert


def adc_check(fpga:FPGAControl):
    '''Test ADC communication.

    Returns
    -------
    cert : dict
        Dictionary of diagnostic results.
    '''
    fpga.adc_setting.set_digital_on()
    fpga.adc_setting.set_test_pattern(ADSChannel.A, ADSTestPattern.RAMP)
    fpga.adc_setting.set_test_pattern(ADSChannel.B, ADSTestPattern.RAMP)

    result_a = fpga.count0.health_check()
    diff_a = fpga.count0.line_diff()
    result_b = fpga.count1.health_check()
    diff_b = fpga.count1.line_diff()

    cert = {'result': result_a & result_b,
            'result_a': result_a,
            'result_b': result_b,
            'diff_a': diff_a,
            'diff_b': diff_b}

    fpga.adc_setting.set_test_pattern(ADSChannel.A, ADSTestPattern.NORMAL)
    fpga.adc_setting.set_test_pattern(ADSChannel.B, ADSTestPattern.NORMAL)
    fpga.adc_setting.set_digital_off()

    return cert


def eprint(*args, **kwargs):
    '''Error print.'''
    print(*args, **kwargs, file=sys.stderr)


def print_record(adc_result, dac_result):
    '''Print medical record.'''
    if adc_result['result'] & dac_result['result']:
        print('Healthy.')
        return

    # Medical record.
    if not adc_result['result']:
        eprint('Communication error in ADC.')
        eprint('===========================')

        if adc_result['result_a']:
            eprint(f'Ch A: {0:016b}')
        else:
            eprint(f'Ch A: {adc_result["diff_a"]:016b}')

        eprint()

        if adc_result['result_b']:
            eprint(f'Ch B: {0:016b}')
        else:
            eprint(f'Ch B: {adc_result["diff_b"]:016b}')

        eprint()

    if not dac_result['result']:
        eprint('Communication error in DAC.')
        eprint('===========================')

        eprint(f'{dac_result["wrong_bits"]:08b}')


def adc_score(adc_result):
    '''ADC health score value.

    Returns
    -------
    score : float
        0 to 1 value with 0 when healthy.
    '''
    count_a = f"{adc_result['diff_a']:016b}".count('1')
    score_a = 0 if adc_result['result_a'] else count_a/ADC_DATA_WIDTH/2
    count_b = f"{adc_result['diff_b']:016b}".count('1')
    score_b = 0 if adc_result['result_b'] else count_b/ADC_DATA_WIDTH/2

    return score_a + score_b


def dac_score(dac_result):
    '''DAC health score value.

    Returns
    -------
    score : float
        0 to 1 value with 0 when healthy.
    '''
    return f"{dac_result['wrong_bits']:08b}".count('1')/8


def main():
    '''Main function.'''
    parser = ArgumentParser()
    parser.add_argument('-ip', '--ip_address',
                    type=str,
                    default=IP_ADDRESS_DEFAULT,
                    help=f'IP-v4 address of target SiTCP. (default={IP_ADDRESS_DEFAULT})')

    parser.add_argument('command', choices=['diagnosis', 'adjust'],
                        default='diganosis')

    args = parser.parse_args()
    fpga = FPGAControl(ip_address=args.ip_address)

    adc_result = adc_check(fpga)
    dac_result = dac_check(fpga)

    print_record(adc_result, dac_result)

    if args.command == 'diagnosis':
        return

    score = adc_score(adc_result) + dac_score(dac_result)
    if score == 0:
        print('Completed.')
        return

    # Phase reconfiguration.
    scores = {}

    for angle in range(0, 360, 45):
        fpga.clock.phase_reconf(0, angle)
        fpga.init()

        adc_result = adc_check(fpga)
        dac_result = dac_check(fpga)

        print_record(adc_result, dac_result)
        score = adc_score(adc_result) + dac_score(dac_result)
        if score == 0:
            print('Completed.')
            return
        scores[angle] = score

    angle_best = max(angle, key=angle.get)
    fpga.clock.phase_reconf(0, angle_best)
    fpga.init()

    adc_result = adc_check(fpga)
    dac_result = dac_check(fpga)

    # Lane adjust
    # To be implemented.

    print_record(adc_result, dac_result)


if __name__ == '__main__':
    main()
