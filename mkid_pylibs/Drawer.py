# -*- coding: utf-8 -*-

import warnings
import matplotlib.pyplot as plt

################################################################
# class for drawing data
################################################################

class Drawer(object):
    """
    A mix-in class for plot data
    """
    def _getdata(self,name):
        if hasattr(self, name):
            return getattr(self, name)
        raise RuntimeError(f'ERROR:: No {name} in data({self.__class__.__name__})')

    def _getlabel(self,name):
        if getattr(self.__class__, name).__doc__:
            return getattr(self.__class__, name).__doc__
        return ''

    def draw(self, xname='x', yname='amplitude', ax=None, *args, **kws):
        """
        plot data from `array_data`
        """
        x_data = self._getdata(xname)
        y_data = self._getdata(yname)

        if ax is None:
            ax = plt.gca()

        ax.plot(x_data, y_data, *args, **kws)

        label = self._getlabel(xname)
        prev = ax.get_xlabel()
        if not prev:
            ax.set_xlabel(label)
        elif prev != label:
            warnings.warn(f'xlabel "{label}" is not macthed to current label: "{prev}"')

        label = self._getlabel(yname)
        prev = ax.get_ylabel()
        if not prev:
            ax.set_ylabel(label)
        elif label != "" and prev != label:
            warnings.warn(f'ylabel "{label}" is not macthed to current label: "{prev}"')

        return ax

