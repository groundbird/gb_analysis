#!/usr/bin/env python

dname = '/home/pi/logger/data/rotation_encoder/'
buff_size = 4096

from os import listdir
from os.path import join, getsize
from numpy import array, concatenate
from struct import unpack
from lzma import LZMAFile # for xz compressed file

def find_latest_file():
    path = dname
    for i in range(4):
        names = listdir(path)
        names.sort()
        path = join(path, names[-1])
        #print path
        pass
    return path

class rot_reader(object):
    def __init__(self, fname):
        self.fname = fname
        if fname.split('.')[-1] == 'xz':
            self.f = LZMAFile(fname, 'rb')
        else:
            self.f = open(fname, 'rb')
            pass

        self.hsize = int(self.f.readline()) # size of header
        self.psize = 8                      # size of packet

        # check file type
        filetype = self.f.readline().strip().decode('utf-8')
        if filetype != 'rotary_encoder_data_for_GB':
            raise Exception

        # read header
        self.header = {}
        for line in self.f:
            line = line.decode('utf-8')
            if line == '\n': break
            wds = line.split(':')
            self.header[wds[0].strip()] = wds[1].strip()
            pass

        # version check
        if not 'version' in self.header:
            raise Exception
        if not self.header['version'] in ['2017120501']:
            raise Exception

        # read header by raw data
        self.f.seek(0)
        self.header_raw = self.f.read(self.hsize)

        # reset of iterator
        self.seek(0)
        pass

    def __iter__(self):
        return self

    def get_size(self):
        if self.fname.split('.')[-1] == 'xz':
            ret = 0
            self.seek(0)
            rsize = buff_size * self.psize
            while True:
                bcnt = len(self.f.read(rsize))
                ret += bcnt
                if bcnt == 0: break
                pass
            self.seek(0)
            ret /= self.psize
        else:
            ret = (getsize(self.fname) - self.hsize) / self.psize
        return ret

    def _convert_packet(self, packets):
        assert self.psize == 8
        l = len(packets)
        if l == 0: return []
        if l % 8 != 0:
            raise Exception
        d = array(unpack('Q' * int(l/8), packets))
        t = d >> (8*2)
        v = d & 0xffff
        return list(zip(t, v))

    def get_data(self, start, end):
        self.seek(start)
        return self.read(end - start)

    def get_last_data(self, data_length):
        end = self.get_size()
        start = max([0, end - data_length])
        return self.get_data(start, end)

    def seek(self, pos = 0):
        self.f.seek(self.hsize + pos * self.psize)
        self.rawbuff = b''
        self.buff = []
        self.buff_it = 0
        return

    def read(self, dnum):
        ret = []
        while len(ret) < dnum:
            if len(self.buff) <= self.buff_it:
                self.rawbuff += self.f.read(buff_size * self.psize - len(self.rawbuff))
                rrnum = int(len(self.rawbuff) / self.psize * self.psize)
                self.buff = self._convert_packet(self.rawbuff[:rrnum])
                self.buff_it = 0
                self.rawbuff = self.rawbuff[rrnum:]
                pass
            if len(self.buff) == 0: break
            rnum = min([dnum - len(ret), len(self.buff) - self.buff_it])
            tmp = self.buff[self.buff_it : self.buff_it + rnum]
            self.buff_it += rnum
            if len(ret) is 0:
                ret = tmp
            else:
                concatenate((ret, tmp), axis=0)
            pass
        return ret

    def __next__(self):
        ret = self.read(1)
        if len(ret) == 0: raise StopIteration()
        return ret[0]

    pass


def test1():
    f = rot_reader(find_latest_file())
    print(f.get_size())
    print(f.header)
    d = f.get_last_data(10)
    #print d.shape
    print(d)
    return

def test2():
    f = rot_reader('data/2018/06/08/rot_2018-0608-151641.dat.xz')
    #f = rot_reader('data/2018/06/08/rot_2018-0608-141411.dat.xz')
    print(f.get_size())
    print(f.header)
    d = f.get_last_data(10)
    #print d.shape
    print(d)
    return

def test3():
    cnt = 10
    for d in rot_reader(find_latest_file()):
        print(d)
        cnt -= 1
        if cnt == 0: break
        pass
    return

def test4():
    f = rot_reader('data/2018/06/08/rot_2018-0608-151641.dat.xz')
    f.seek(1248000)
    for d in f:
        print(d)
        pass
    return

if __name__ == '__main__':
    test1()
    pass
