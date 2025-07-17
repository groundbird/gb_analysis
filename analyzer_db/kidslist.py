#!/usr/bin/env python
import json
import numpy as np

class KidsList:
    def __init__(self, path):
        with open(path) as f:
            self._kldict = json.load(f)

        self.kids_index   = self._kldict['kids']
        self.kids_power   = self._kldict['kids_power']
        self.kids_freqs   = self._kldict['kids_freqs']

        if 'kids_amps' in self._kldict.keys():
            self.kids_amps = self._kldict['kids_amps']
        else:
            self.kids_amps = [1.0] * len(self.kids_freqs)

        if 'kids_phases' in self._kldict.keys():
            self.kids_phases = self._kldict['kids_phases']
        else:
            self.kids_phases = [0.] * len(self.kids_freqs)


        self.blinds_index = self._kldict['blinds']
        self.blinds_power = self._kldict['blinds_power']

        if 'blinds_freqs' in self._kldict.keys():
            self.blinds_freqs = self._kldict['blinds_freqs']
        elif 'blinds_relfreqs' in self._kldict.keys():
            if len(self._kldict['blinds_relfreqs']) == 0:
                self.blinds_freqs = []
            else:
                if 'blinds_relindex' in self._kldict.keys():
                    self.blinds_relindex = self._kldict['blinds_relindex']
                else:
                    self.blinds_relindex = self.kids_index
                if len(self.blinds_relindex) != len(self._kldict['blinds_relfreqs']):
                    raise Exception('Invalid KidsList: different length: blinds_relindex <=> blinds_relfreqs')
                self.blinds_freqs = [self._kldict['kids_freqs'][np.where(np.array(self.kids_index) == ik)[0][0]]+fb for ik,fb in zip(self.blinds_relindex,self._kldict['blinds_relfreqs'])]

        if 'blinds_anaindex' in self._kldict.keys():
            self.blinds_anaindex = self._kldict['blinds_anaindex']
        else:
            self.blinds_anaindex = self.blinds_index

        if 'blinds_amps' in self._kldict.keys():
            self.blinds_amps = self._kldict['blinds_amps']
        else:
            self.blinds_amps = [1.0] * len(self.blinds_freqs)

        if 'blinds_phases' in self._kldict.keys():
            self.blinds_phases = self._kldict['blinds_phases']
        else:
            self.blinds_phases = [0.] * len(self.blinds_freqs)

        if 'kids_fitconfig' in self._kldict.keys():
            self.fit_conf = self._kldict['kids_fitconfig']
            for i,x in enumerate(self.fit_conf):
                if type(x) is not str: continue
                if x.lower() == "default":
                    self.fit_conf[i] = True
        else:
            self.fit_conf = [True] * len(self.kids_freqs)

        self.sg_freq      = self._kldict['sg_freq']

        self.array_name   = self._kldict['array_name']

        ############################################################### Health check
        ##################################################### Collision
        for bi in self.blinds_index:
            if bi in self.kids_index:
                raise Exception('Invalid KidsList: index collision: blinds <=> kids')

        if len(set(self.kids_index)) != len(self.kids_index):
            raise Exception('Invalid KidsList: index collision: kids')

        if len(set(self.blinds_index)) != len(self.blinds_index):
            raise Exception('Invalid KidsList: index collision: blinds')

        ##################################################### Completeness
        ind_max = max(self.kids_index)
        if len(self.blinds_index)>0:
            ind_max = max(max(self.kids_index), max(self.blinds_index))

        freqlist = [None]*(ind_max + 1)
        for ki, kf in zip(self.kids_index, self.kids_freqs):
            freqlist[ki] = kf

        for bi, bf in zip(self.blinds_index, self.blinds_freqs):
            freqlist[bi] = bf

        if None in freqlist:
            raise Exception('Invalid KidsList: incomplete.')

        self.tone_freqs = freqlist

        amplist = [None]*(ind_max + 1)
        for ki, kf in zip(self.kids_index, self.kids_amps):
            amplist[ki] = kf

        for bi, bf in zip(self.blinds_index, self.blinds_amps):
            amplist[bi] = bf

        if None in amplist:
            raise Exception('Invalid KidsList: incomplete.')

        self.tone_amps = amplist

        phaselist = [None]*(ind_max + 1)
        for ki, kf in zip(self.kids_index, self.kids_phases):
            phaselist[ki] = kf

        for bi, bf in zip(self.blinds_index, self.blinds_phases):
            phaselist[bi] = bf

        if None in phaselist:
            raise Exception('Invalid KidsList: incomplete.')

        self.tone_phases = phaselist

if __name__ == '__main__':
    import sys
    KidsList(sys.argv[1])
