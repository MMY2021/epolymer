#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function, division, absolute_import
from netCDF4 import Dataset
import numpy as np
from polymer.block import Block
from datetime import datetime
from polymer.ancillary import Ancillary_NASA
from polymer.common import L2FLAGS
from collections import OrderedDict
from polymer.utils import raiseflag
from polymer.level1 import Level1_base
from os.path import dirname, join
import pandas as pd


def filled(A, ok=None, fill_value=0):
    """
    Returns a filled from a filled or masked array, use fill_value
    modifies ok (if provided) to take this mask into account
    """
    if hasattr(A, 'filled'):
        # masked array: returns filled array
        if ok is not None:
            ok &= ~A.mask
        return A.filled(fill_value=fill_value)
    else:
        # filled array: does nothing
        return A


class Level1_NASA(Level1_base):
    '''
    Interface to NASA level-1C files
    Applies to sensors:
        - SeaWiFS
        - VIIRS
        - MODIS
    '''
    def __init__(self, filename, sensor=None, blocksize=(500, 400),
                 sline=0, eline=-1, scol=0, ecol=-1, ancillary=None,
                 altitude=0.):
        self.sensor = sensor
        self.filename = filename
        self.root = Dataset(filename)
        self.altitude = altitude
        lat = self.root.variables['latitude']
       
        totalheight, totalwidth = lat.shape
        self.blocksize = blocksize


        self.init_shape(
                totalheight=totalheight,
                totalwidth=totalwidth,
                sline=sline,
                eline=eline,
                scol=scol,
                ecol=ecol)

        # init dates
        self.__read_date()

        self.init_spectral_info()



    def init_spectral_info(self):
        # NOTE: central wavelengths are from SeaDAS

        if self.sensor == 'MODIS':
            bands = [412,443,469,488,531,547,555,645,667,678,748,859,869,1240,1640,2130]
            
        elif self.sensor == 'CASI':
            bands = [412,443,469,488,531,547,555,645,667,678,748,859,869,1240,1640,2130]
       
        elif self.sensor == 'SeaWiFS':
       
            bands = [412,443,490,510,555,670,765,865]
        elif self.sensor in ['VIIRS', 'VIIRSN']:
     
            bands = [410,443,486,551,671,745,862,1238,1601,2257]
        elif self.sensor == 'VIIRSJ1':
     
            bands = [411,445,489,556,667,746,868,1238,1604,2258]
        else:
            raise Exception('Invalid sensor "{}"'.format(self.sensor))

        self.central_wavelength = dict([(b, float(b)) for b in bands])


    def read_block(self, size, offset, bands):

        nbands = len(bands)
        size3 = size + (nbands,)
        (ysize, xsize) = size

        SY = slice(offset[0]+self.sline, offset[0]+self.sline+size[0])
        SX = slice(offset[1]+self.scol , offset[1]+self.scol+size[1])

        # initialize block
        block = Block(offset=offset, size=size, bands=bands)

        # read lat/lon
        block.latitude = filled(self.root.variables[
                'latitude'][SY, SX], fill_value=-999.)
        block.longitude = filled(self.root.variables[
                'longitude'][SY, SX], fill_value=-999.)

        ok = block.latitude > -90.

        # read geometry
        # note: we disactivate automasking because of bad formatting of SeaWiFS L1C, for which azimuth angles >180 are masked
        block.sza = filled(self.root.variables['solz'][SY, SX], ok=ok)
        block.vza = filled(self.root.variables['senz'][SY, SX], ok=ok)

        saa = self.root.variables['sola']
      
        saa.set_auto_mask(False)
        block.saa = filled(saa[SY, SX])
       
        vaa = self.root.variables['sena']
        vaa.set_auto_mask(False)
        block.vaa = filled(vaa[SY, SX])
        
        wind_speed = self.root.variables['wind_speed']
        wind_speed.set_auto_mask(False)
        block.wind_speed = filled(wind_speed[SY, SX])
        #print(block.vza)

        block.Rtoa = np.zeros(size3) + np.NaN
        block.t_no2 = np.zeros(size3) + np.NaN
        block.trans_o3 = np.zeros(size3) + np.NaN
        block.Rmol = np.zeros(size3) + np.NaN
        block.Rmolgli = np.zeros(size3) + np.NaN
        block.Tmol = np.zeros(size3,dtype='float32') + np.NaN
        
        for iband, band in enumerate(bands):
            Rtoa = filled(self.root.groups['rapp'].variables[
                    'rapp_{}'.format(band)][SY, SX], ok=ok)
            block.Rtoa[:,:,iband] = Rtoa
            t_no2 = filled(self.root.groups['t_no2'].variables[
                    't_no2_{}'.format(band)][SY, SX], ok=ok)
            block.t_no2[:,:,iband] = t_no2
            trans_o3 = filled(self.root.groups['t_o3'].variables[
                    't_o3_{}'.format(band)][SY, SX], ok=ok)
            block.trans_o3[:,:,iband] = trans_o3
            Rmol = filled(self.root.groups['Rmol'].variables[
                    'Rmol_{}'.format(band)][SY, SX], ok=ok)
            block.Rmol[:,:,iband] = Rmol
            Rmolgli = filled(self.root.groups['Rmolgli'].variables[
                    'Rmolgli_{}'.format(band)][SY, SX], ok=ok)
            block.Rmolgli[:,:,iband] = Rmolgli
            Tmol = filled(self.root.groups['Ray_tra'].variables[
                    'Tmol_{}'.format(band)][SY, SX], ok=ok)
            block.Tmol[:,:,iband] = Tmol
           
        '''    
        for iband, band in enumerate(bands):
            Rtoa = filled(self.root.variables[
                    'rhot_{}'.format(band)][SY, SX], ok=ok)
            block.Rtoa[:,:,iband] = Rtoa
            t_no2 = filled(self.root.variables[
                    't_no2_{}'.format(band)][SY, SX], ok=ok)
            block.t_no2[:,:,iband] = t_no2
            trans_o3 = filled(self.root.variables[
                    't_o3_{}'.format(band)][SY, SX], ok=ok)
            block.trans_o3[:,:,iband] = trans_o3
            Rmol = filled(self.root.variables[
                    'Rmol_{}'.format(band)][SY, SX], ok=ok)
            block.Rmol[:,:,iband] = Rmol
            Rmolgli = filled(self.root.variables[
                    'Rmolgli_{}'.format(band)][SY, SX], ok=ok)
            block.Rmolgli[:,:,iband] = Rmolgli
            Tmol = filled(self.root.variables[
                    'Tmol_{}'.format(band)][SY, SX], ok=ok)
            block.Tmol[:,:,iband] = Tmol
        '''   
   
       # print(block.Tmol)
        #print(block.Tmol.dtype)
         

        # bitmask
        block.bitmask = np.zeros(size, dtype='uint16')

        ok &= block.Rtoa[:,:,0] >= 0
        raiseflag(block.bitmask, L2FLAGS['L1_INVALID'], ~ok)
        

        block.jday = self.date().timetuple().tm_yday
        block.month = self.date().timetuple().tm_mon

        block.wavelen = np.zeros(size3, dtype='float32') + np.NaN
        block.cwavelen = np.zeros(nbands, dtype='float32') + np.NaN
        for iband, band in enumerate(bands):
            block.wavelen[:,:,iband] = self.central_wavelength[band]
            block.cwavelen[iband] = self.central_wavelength[band]

        return block

    def __read_date(self):
        try:
            dstart = datetime.strptime(self.root.getncattr('time_coverage_start'),
                                  '%Y-%m-%dT%H:%M:%S.%fZ')
        except ValueError: # try again without decimal part
            dstart = datetime.strptime(self.root.getncattr('time_coverage_start'),
                                  '%Y-%m-%dT%H:%M:%S')
        try:
            dstop = datetime.strptime(self.root.getncattr('time_coverage_end'),
                                  '%Y-%m-%dT%H:%M:%S.%fZ')
        except ValueError: # try again without decimal part
            dstop = datetime.strptime(self.root.getncattr('time_coverage_end'),
                                  '%Y-%m-%dT%H:%M:%S')

        self.dstart = dstart
        self.dstop = dstop

    def date(self):
        return self.dstart + (self.dstop - self.dstart)//2

    def attributes(self, datefmt):
        '''
        Returns level1 attributes

        dates are formatted to string using datefmt
        '''
        attr = OrderedDict()
        attr['l1_filename'] = self.filename
        attr['start_time'] = self.dstart.strftime(datefmt)
        attr['stop_time'] = self.dstop.strftime(datefmt)
        attr['central_wavelength'] = self.central_wavelength


        return attr

    def __enter__(self):
        return self

    def __exit__(self, *args):
        pass


class Level1_VIIRS(Level1_NASA):
    ''' Interface to VIIRS Level-1C '''
    def __init__(self, filename, **kwargs):
        root = Dataset(filename)
        platform = root.getncattr('platform')
        sensor = {
                'Suomi-NPP':'VIIRSN',
                'JPSS-1': 'VIIRSJ1',
                }[platform]
        super(self.__class__, self).__init__(
                filename, sensor=sensor, **kwargs)

class Level1_SeaWiFS(Level1_NASA):
    ''' Interface to SeaWiFS Level-1C '''
    def __init__(self, filename, **kwargs):
        super(self.__class__, self).__init__(
                filename, sensor='SeaWiFS', **kwargs)

class Level1_MODIS(Level1_NASA):
    ''' Interface to MODIS Level-1C '''
    def __init__(self, filename, **kwargs):
        super(self.__class__, self).__init__(
                filename, sensor='MODIS', **kwargs)
                
class Level1_CASI(Level1_NASA):
    ''' Interface to CASI Level-1C '''
    def __init__(self, filename, **kwargs):
        super(self.__class__, self).__init__(
                filename, sensor='CASI', **kwargs)

