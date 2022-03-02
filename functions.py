#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 25 10:47:20 2022

@author: pst019
"""


"""make the sine scale"""
import matplotlib.scale as mscale
import matplotlib.transforms as mtransforms
import matplotlib.ticker as ticker

import numpy as np

class SineScale(mscale.ScaleBase):
    """
    ScaleBase class for generating square root scale.
    """
 
    name = 'sine'
 
    def __init__(self, axis, **kwargs):
        # note in older versions of matplotlib (<3.1), this worked fine.
        # mscale.ScaleBase.__init__(self)

        # In newer versions (>=3.1), you also need to pass in `axis` as an arg
        mscale.ScaleBase.__init__(self, axis)
 
    def set_default_locators_and_formatters(self, axis):
        axis.set_major_locator(ticker.AutoLocator())
        axis.set_major_formatter(ticker.ScalarFormatter())
        axis.set_minor_locator(ticker.NullLocator())
        axis.set_minor_formatter(ticker.NullFormatter())
 
    def limit_range_for_scale(self, vmin, vmax, minpos):
        return  max(-90., vmin), min(90, vmax)
 
    class SineTransform(mtransforms.Transform):
        input_dims = 1
        output_dims = 1
        is_separable = True
 
        def transform_non_affine(self, a): 
            # return np.array(a)**0.5
            return np.sin(np.deg2rad(a))
 
        def inverted(self):
            return SineScale.InvertedSineTransform()
 
    class InvertedSineTransform(mtransforms.Transform):
        input_dims = 1
        output_dims = 1
        is_separable = True
 
        def transform(self, a):
            return np.rad2deg(np.arcsin(a))
 
        def inverted(self):
            return SineScale.SineTransform()
 
    def get_transform(self):
        return self.SineTransform()
 

mscale.register_scale(SineScale)




def sign_change_ds(ds1):
    #ds1 is the data array
    #dslat is just ds.lat
    from scipy.signal import argrelextrema
    # from scipy import ndimage
    latsmooth= len(ds1.lat)//10 #10% of the latitudes are smoothed in the filter

    if len(ds1.shape)== 2: ds1= ds1.mean(axis=1)
    ds1_s= np.sign(ds1)
    sign_change= (np.roll(ds1_s, 1)- ds1_s != 0).astype(int)
    sign_change[0]= 0
    # return sign_change

    # sign_change_mean= sign_change.mean(axis= 1)
    locmax= argrelextrema(sign_change.to_numpy(), np.greater)
    sign_change= sign_change.rolling(lat= latsmooth, center= True).max().fillna(0)
    # sign_change[sign_change<1E-7]=0
    return ds1.lat.values[locmax], sign_change
