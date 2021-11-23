"""
napari-STPT reads zarr files and displays them
"""

import os
import sys
import napari
import numpy as np
import xarray as xr
from qtpy import QtCore, QtWidgets
import SimpleITK as sitk
from scipy import ndimage, stats
import cv2
from stardist.models import StarDist2D
from napari_animation import AnimationWidget
from PIL import Image, ImageDraw
import random
# import skimage.io
# from stardist.models import StarDist3D
from csbdeep.utils import normalize
#from naparimovie import Movie
import math
import gc
import json


class NapariSTPT:

    def __init__(self):
        self.scroll = None
        self.image_slice = None
        self.normalize_value = None

        self.pixel_size = None
        self.viewer = None
        self.comboBoxPath = None
        self.comboBoxResolution = None
        self.cb_C1 = None
        self.cb_C2 = None
        self.cb_C3 = None
        self.cb_C4 = None
        self.search_folder = None
        self.aligned_1 = None
        self.aligned_2 = None
        self.aligned_3 = None
        self.aligned_4 = None

        self.spinN = None
        self.maxSizeN = None
        self.thresholdN = None
        self.slice_names = None

        self.align_x = None
        self.align_y = None
        self.corrected_align_x = None
        self.corrected_align_y = None
        self.ds1 = None
        self.im1 = None
        self.im2 = None
        self.im4 = None
        self.im8 = None
        self.im16 = None
        self.im32 = None

        #self.movie = None
        self.optical_slices = None
        self.nr_optical_slices = None
        self.optical_slices_available = None
        self.cb_correct_brightness_optical_section = None
        self.shape = None

        self.crop_start_ratio_x = None
        self.crop_size_ratio_x = None
        self.crop_start_ratio_y = None
        self.crop_size_ratio_y = None

        self.crop = False
        self.old_method = True

        self.start_slice = None
        self.end_slice = None

        self.image_translation = None
        self.m_slice_spacing = None
        self.slice_spacing = None
        
    def CropOriginal(self, volume):

        threshold = float(self.thresholdN.text())

        tissue = volume > threshold
        objs = ndimage.find_objects(tissue)
        maxX = int(objs[0][0].stop)
        minX = int(objs[0][0].start)
        maxY = int(objs[0][1].stop)
        minY = int(objs[0][1].start)
        maxZ = int(objs[0][2].stop)
        minZ = int(objs[0][2].start)

        print("cropping to {}-{},{}-{},{}-{}".format(minX,
              maxX, minY, maxY, minZ, maxZ))

        volume = volume[minX:maxX, minY:maxY, minZ:maxZ]

        return volume
        
    def Load(self, text):

        try:
            self.viewer.layers.remove('C1')
            del self.aligned_1
        except Exception:
            pass
        try:
            self.viewer.layers.remove('C2')
            del self.aligned_2
        except Exception:
            pass
        try:
            self.viewer.layers.remove('C3')
            del self.aligned_3
        except Exception:
            pass
        try:
            self.viewer.layers.remove('C4')
            del self.aligned_4
        except Exception:
            pass

        gc.collect()

        print(self.search_folder.text() +
              str(self.comboBoxPath.currentText()) + '/mos.zarr')

        try:
            self.ds1 = xr.open_zarr(self.search_folder.text(
            ) + str(self.comboBoxPath.currentText())+'/mos.zarr', consolidated=True)
            print("consolidated")
        except Exception:
            print("none-consolidated")
            self.ds1 = xr.open_zarr(self.search_folder.text() +
                                    str(self.comboBoxPath.currentText())+'/mos.zarr')

        # ds2 = xr.open_zarr(str(self.comboBoxPath.currentText())+'/mos.zarr', group='l.2')
        # ds4 = xr.open_zarr(str(self.comboBoxPath.currentText())+'/mos.zarr', group='l.4')
        # ds8 = xr.open_zarr(str(self.comboBoxPath.currentText())+'/mos.zarr', group='l.8')
        # ds32 = xr.open_zarr(str(self.comboBoxPath.currentText())+'/mos.zarr', group='l.32')
        
        # print("self.ds1.attrs")
        # print(self.ds1.attrs)
        # print("self.ds1.keys()")
        # print(self.ds1.keys()["Coordinates"])
        # print("self.ds1.attrs.keys()")
        # print(self.ds1.attrs.keys())
        # print("self.ds1['S002'].attrs")
        # print(self.ds1['S002'].attrs)
        # print(self.ds1[f"S{slice_no:03d}"].attrs)


        #print(self.ds1['S002'].attrs['offsets']['x'])
        # print(self.ds1['S002'].attrs['bscale'])
        # print(json.loads(self.ds1['S002'].attrs['scale'])["x"])
        # print(json.loads(self.ds1['S002'].attrs['scale'])["y"])
        # print(json.loads(self.ds1['S002'].attrs['scale'])["z"])

        #print(type(self.ds1['S002'].attrs['scale']))
        #print(self.ds1['S002'].attrs['scale']["y"])
        #print(self.ds1['S002'].attrs['scale']['z'])
        #print(self.ds1['S001'].attrs['raw_meta'][0][0]['zres'])
        #print(self.ds1['S001'].attrs.keys())
        #print(self.ds1.attrs["multiscale"]['datasets'][1]['level'])
        #print(self.ds1.attrs["multiscale"]['datasets'][1]['path'])
        # print(self.slice_names)

        #old_method = True

        if self.old_method:
            #Old method
            bscale = self.ds1.attrs['bscale']
            bzero = self.ds1.attrs['bzero']
            self.slice_names = self.ds1.attrs['cube_reg']['slice']
            self.align_x = self.ds1.attrs['cube_reg']['abs_dx']
            self.align_y = self.ds1.attrs['cube_reg']['abs_dy']
            #print(f"self.align_x length {len(self.align_x)}")
            #print("")
            #print(f"self.align_x {self.align_x}")
            #print(f"self.align_y {self.align_y}")
            #print("")
        else:
            bscale = self.ds1['S001'].attrs['bscale']
            bzero = self.ds1['S001'].attrs['bzero']

        # print("size: {}".format(self.ds1.dims))

        output_resolution = float(self.pixel_size.text())
        #resolution = int(self.comboBoxResolution.currentText())
        resolution = 32
        index = 5
        if (output_resolution / 0.5) < 32:
            resolution = 16
            index = 4
        if (output_resolution / 0.5) < 16:
            resolution = 8
            index = 3
        if (output_resolution / 0.5) < 8:
            resolution = 4
            index = 2
        if (output_resolution / 0.5) < 4:
            resolution = 2
            index = 1
        if (output_resolution / 0.5) < 2:
            resolution = 1
            index = 0
        
        if self.old_method:
            #print("loading at resolution {}".format(resolution))
            #index = self.comboBoxResolution.currentIndex()
            #print("group: " + self.ds1.attrs["multiscale"]['datasets'][index]['path'])
            gr = self.ds1.attrs["multiscale"]['datasets'][index]['path']
            ds = xr.open_zarr(self.search_folder.text() + str(self.comboBoxPath.currentText())+'/mos.zarr', group=gr)
        else:
            print("resolution")
            print(resolution)
            print(json.loads(self.ds1.attrs["multiscale"])['datasets'][index]['path'])
            gr = json.loads(self.ds1.attrs["multiscale"])['datasets'][index]['path']
            print(gr)
            ds = xr.open_zarr(self.search_folder.text() + str(self.comboBoxPath.currentText())+'/mos.zarr', group=gr)

        self.optical_slices_available = len(ds.z)
        self.optical_slices = self.optical_slices_available
        print(f"optical_slices available: {self.optical_slices}")
        
        self.slice_spacing = float(self.m_slice_spacing.text())

        if self.optical_slices > 1:
            optical_section_spacing = self.ds1['S001'].attrs['raw_meta'][0][0]['zres']
            print(f"optical_slices zres: {optical_section_spacing}")
            expected_nr_of_slices = round(self.slice_spacing / optical_section_spacing)
            if self.optical_slices > expected_nr_of_slices:
                self.optical_slices = expected_nr_of_slices

        print(f"optical_slices used: {self.optical_slices}")

        self.corrected_align_x = []
        for index, value in enumerate(self.align_x):
            if (index % self.optical_slices_available) < self.optical_slices:
                self.corrected_align_x.append(value)

        #self.align_x = self.corrected_align_x
        #print(f"self.align_x shape {len(self.corrected_align_x)}")

        self.corrected_align_y = []
        for index, value in enumerate(self.align_y):
            if (index % self.optical_slices_available) < self.optical_slices:
                self.corrected_align_y.append(value)

        #self.align_y = corrected_align_y
        print(f"align_y  {self.align_y}")
        print("")
        print(f"corrected_align_y {self.corrected_align_y}")
        print("")

        for i in range(0,self.optical_slices):
            self.corrected_align_y[i::self.optical_slices] = self.corrected_align_y[::self.optical_slices]
            
        print(f"corrected_align_y {self.corrected_align_y}")
        print("")

        with napari.gui_qt():
            if self.cb_C1.isChecked():
                print("loading C1")

                volume_1_temp = (ds.sel(channel=1, type='mosaic', z=0).to_array(
                ).data * bscale + bzero).astype(dtype=np.float32)

                if self.start_slice.text() == "":
                    start_slice_number = 0
                    chop_bottom = 0
                else:
                    start_slice_number = int(math.floor(float(self.start_slice.text())/float(self.optical_slices)))
                    chop_bottom = int(self.start_slice.text()) - (self.optical_slices * start_slice_number) 
                if self.end_slice.text() == "":
                    end_slice_number = volume_1_temp.shape[0]-1
                    chop_top = 0
                else:
                    end_slice_number = int(math.floor(float(self.end_slice.text())/float(self.optical_slices)))
                    chop_top = (self.optical_slices * (end_slice_number + 1) -1) - int(self.end_slice.text()) 

                number_of_slices = end_slice_number - start_slice_number + 1

                if self.crop:
                    size_y = int(math.floor(
                        volume_1_temp.shape[1]/self.crop_size_ratio_x))
                    size_z = int(math.floor(
                        volume_1_temp.shape[2]/self.crop_size_ratio_y))
                    start_y = int(math.floor(
                        volume_1_temp.shape[1]/self.crop_start_ratio_x))
                    start_z = int(math.floor(
                        volume_1_temp.shape[2]/self.crop_start_ratio_y))
                else:
                    size_y = int(math.floor(volume_1_temp.shape[1]))
                    size_z = int(math.floor(volume_1_temp.shape[2]))
                    start_y = 0
                    start_z = 0

                volume_1 = np.zeros((self.optical_slices*number_of_slices, size_y, size_z), dtype=np.float32)
                volume_1[0::self.optical_slices, :, :] = volume_1_temp[start_slice_number:end_slice_number+1, start_y:start_y+size_y, start_z:start_z+size_z]

                for optical_slice in range(1, self.optical_slices):
                    volume_1_temp = (ds.sel(channel=1, type='mosaic', z=optical_slice).to_array(
                    ).data * bscale + bzero).astype(dtype=np.float32)
                    volume_1[optical_slice::self.optical_slices, :, :] = volume_1_temp[start_slice_number:end_slice_number+1, start_y:start_y+size_y, start_z:start_z+size_z]

                # self.aligned_2 = volume_2
                if self.cb_correct_brightness_optical_section.isChecked():
                    print("correcting brightness of optical sections C2")
                    self.Normalize_slices(volume_1, self.optical_slices)

                print("aligning C1")
                self.aligned_1 = self.Align(volume_1, resolution, output_resolution, start_slice_number*self.optical_slices)
                #self.aligned_1 = volume_1
                self.aligned_1 = self.aligned_1[chop_bottom:self.aligned_1.shape[0]-chop_top,:,:]
                self.shape = self.aligned_1.shape
                self.viewer.add_image([self.aligned_1], name='C1', scale=(float(self.slice_spacing)/float(self.optical_slices), output_resolution, output_resolution), blending='additive', colormap='bop purple', contrast_limits=[0,30])#, translate=(0, int(start_y*output_resolution), int(start_z*output_resolution)))

            if self.cb_C2.isChecked():
                print("loading C2")

                volume_2_temp = (ds.sel(channel=2, type='mosaic', z=0).to_array(
                ).data * bscale + bzero).astype(dtype=np.float32)

                if self.start_slice.text() == "":
                    start_slice_number = 0
                    chop_bottom = 0
                else:
                    start_slice_number = int(math.floor(float(self.start_slice.text())/float(self.optical_slices)))
                    chop_bottom = int(self.start_slice.text()) - (self.optical_slices * start_slice_number) 
                if self.end_slice.text() == "":
                    end_slice_number = volume_2_temp.shape[0]-1
                    chop_top = 0
                else:
                    end_slice_number = int(math.floor(float(self.end_slice.text())/float(self.optical_slices)))
                    chop_top = (self.optical_slices * (end_slice_number + 1) -1) - int(self.end_slice.text()) 

                number_of_slices = end_slice_number - start_slice_number + 1

                if self.crop:
                    size_y = int(math.floor(
                        volume_2_temp.shape[1]/self.crop_size_ratio_x))
                    size_z = int(math.floor(
                        volume_2_temp.shape[2]/self.crop_size_ratio_y))
                    start_y = int(math.floor(
                        volume_2_temp.shape[1]/self.crop_start_ratio_x))
                    start_z = int(math.floor(
                        volume_2_temp.shape[2]/self.crop_start_ratio_y))
                else:
                    size_y = int(math.floor(volume_2_temp.shape[1]))
                    size_z = int(math.floor(volume_2_temp.shape[2]))
                    start_y = 0
                    start_z = 0

                volume_2 = np.zeros((self.optical_slices*number_of_slices, size_y, size_z), dtype=np.float32)
                volume_2[0::self.optical_slices, :, :] = volume_2_temp[start_slice_number:end_slice_number+1, start_y:start_y+size_y, start_z:start_z+size_z]

                for optical_slice in range(1, self.optical_slices):
                    volume_2_temp = (ds.sel(channel=2, type='mosaic', z=optical_slice).to_array(
                    ).data * bscale + bzero).astype(dtype=np.float32)
                    volume_2[optical_slice::self.optical_slices, :, :] = volume_2_temp[start_slice_number:end_slice_number+1, start_y:start_y+size_y, start_z:start_z+size_z]

                # self.aligned_2 = volume_2
                if self.cb_correct_brightness_optical_section.isChecked():
                    print("correcting brightness of optical sections C2")
                    self.Normalize_slices(volume_2, self.optical_slices)


                print("aligning C2")
                self.aligned_2 = self.Align(volume_2, resolution, output_resolution, start_slice_number*self.optical_slices)
                #self.aligned_2 = volume_2
                self.aligned_2 = self.aligned_2[chop_bottom:self.aligned_2.shape[0]-chop_top,:,:]
                self.shape = self.aligned_2.shape
                self.viewer.add_image([self.aligned_2], name='C2', scale=(float(self.slice_spacing)/float(self.optical_slices), output_resolution, output_resolution), blending='additive', colormap='red', contrast_limits=[0,30])#, translate=(0, int(start_y*output_resolution), int(start_z*output_resolution)))

            if self.cb_C3.isChecked():
                print("loading C3")

                volume_3_temp = (ds.sel(channel=3, type='mosaic', z=0).to_array(
                ).data * bscale + bzero).astype(dtype=np.float32)

                if self.start_slice.text() == "":
                    start_slice_number = 0
                    chop_bottom = 0
                else:
                    start_slice_number = int(math.floor(float(self.start_slice.text())/float(self.optical_slices)))
                    chop_bottom = int(self.start_slice.text()) - (self.optical_slices * start_slice_number) 
                if self.end_slice.text() == "":
                    end_slice_number = volume_3_temp.shape[0]-1
                    chop_top = 0
                else:
                    end_slice_number = int(math.floor(float(self.end_slice.text())/float(self.optical_slices)))
                    chop_top = (self.optical_slices * (end_slice_number + 1) -1) - int(self.end_slice.text())

                number_of_slices = end_slice_number - start_slice_number + 1

                if self.crop:
                    size_y = int(math.floor(
                        volume_3_temp.shape[1]/self.crop_size_ratio_x))
                    size_z = int(math.floor(
                        volume_3_temp.shape[2]/self.crop_size_ratio_y))
                    start_y = int(math.floor(
                        volume_3_temp.shape[1]/self.crop_start_ratio_x))
                    start_z = int(math.floor(
                        volume_3_temp.shape[2]/self.crop_start_ratio_y))
                else:
                    size_y = int(math.floor(volume_3_temp.shape[1]))
                    size_z = int(math.floor(volume_3_temp.shape[2]))
                    start_y = 0
                    start_z = 0

                volume_3 = np.zeros((self.optical_slices*number_of_slices, size_y, size_z), dtype=np.float32)
                volume_3[0::self.optical_slices, :, :] = volume_3_temp[start_slice_number:end_slice_number+1, start_y:start_y+size_y, start_z:start_z+size_z]

                for optical_slice in range(1, self.optical_slices):
                    volume_3_temp = (ds.sel(channel=3, type='mosaic', z=optical_slice).to_array(
                    ).data * bscale + bzero).astype(dtype=np.float32)
                    volume_3[optical_slice::self.optical_slices, :, :] = volume_3_temp[start_slice_number:end_slice_number+1, start_y:start_y+size_y, start_z:start_z+size_z]

                # self.aligned_3 = volume_3
                if self.cb_correct_brightness_optical_section.isChecked():
                    print("correcting brightness of optical sections C2")
                    self.Normalize_slices(volume_3, self.optical_slices)

                print("aligning C3")
                self.aligned_3 = self.Align(volume_3, resolution, output_resolution, start_slice_number*self.optical_slices)
                #self.aligned_3 = volume_3
                self.aligned_3 = self.aligned_3[chop_bottom:self.aligned_3.shape[0]-chop_top,:,:]
                self.shape = self.aligned_3.shape
                self.viewer.add_image([self.aligned_3], name='C3', scale=(float(self.slice_spacing)/float(self.optical_slices), output_resolution, output_resolution), blending='additive', colormap='green', contrast_limits=[0,30])#, translate=(0, int(start_y*output_resolution), int(start_z*output_resolution)))

            if self.cb_C4.isChecked():
                print("loading C4")

                volume_4_temp = (ds.sel(channel=4, type='mosaic', z=0).to_array(
                ).data * bscale + bzero).astype(dtype=np.float32)

                if self.start_slice.text() == "":
                    start_slice_number = 0
                    chop_bottom = 0
                else:
                    start_slice_number = int(math.floor(float(self.start_slice.text())/float(self.optical_slices)))
                    chop_bottom = int(self.start_slice.text()) - (self.optical_slices * start_slice_number) 
                if self.end_slice.text() == "":
                    end_slice_number = volume_4_temp.shape[0]-1
                    chop_top = 0
                else:
                    end_slice_number = int(math.floor(float(self.end_slice.text())/float(self.optical_slices)))
                    chop_top = (self.optical_slices * (end_slice_number + 1) -1) - int(self.end_slice.text())

                number_of_slices = end_slice_number - start_slice_number + 1

                if self.crop:
                    size_y = int(math.floor(
                        volume_4_temp.shape[1]/self.crop_size_ratio_x))
                    size_z = int(math.floor(
                        volume_4_temp.shape[2]/self.crop_size_ratio_y))
                    start_y = int(math.floor(
                        volume_4_temp.shape[1]/self.crop_start_ratio_x))
                    start_z = int(math.floor(
                        volume_4_temp.shape[2]/self.crop_start_ratio_y))
                else:
                    size_y = int(math.floor(volume_4_temp.shape[1]))
                    size_z = int(math.floor(volume_4_temp.shape[2]))
                    start_y = 0
                    start_z = 0

                volume_4 = np.zeros((self.optical_slices*number_of_slices, size_y, size_z), dtype=np.float32)
                volume_4[0::self.optical_slices, :, :] = volume_4_temp[start_slice_number:end_slice_number+1, start_y:start_y+size_y, start_z:start_z+size_z]

                for optical_slice in range(1, self.optical_slices):
                    volume_4_temp = (ds.sel(channel=4, type='mosaic', z=optical_slice).to_array(
                    ).data * bscale + bzero).astype(dtype=np.float32)
                    volume_4[optical_slice::self.optical_slices, :, :] = volume_4_temp[start_slice_number:end_slice_number+1, start_y:start_y+size_y, start_z:start_z+size_z]

                # self.aligned_4 = volume_4
                if self.cb_correct_brightness_optical_section.isChecked():
                    print("correcting brightness of optical sections C2")
                    self.Normalize_slices(volume_4, self.optical_slices)

                print("aligning C4")
                self.aligned_4 = self.Align(volume_4, resolution, output_resolution, start_slice_number*self.optical_slices)
                #self.aligned_4 = volume_4
                self.aligned_4 = self.aligned_4[chop_bottom:self.aligned_4.shape[0]-chop_top,:,:]
                self.shape = self.aligned_4.shape
                self.viewer.add_image([self.aligned_4], name='C4', scale=(float(self.slice_spacing)/float(self.optical_slices), output_resolution, output_resolution), blending='additive', colormap='bop blue', contrast_limits=[0,30])#, translate=(0, int(start_y*output_resolution), int(start_z*output_resolution)))

            # create naparimovie object
            #if self.movie is None:
                #self.movie = Movie(myviewer=self.viewer)
                #self.movie.inter_steps = 30

    def LoadInRegion(self, text):

        print(self.viewer.layers['Shapes'].data)
        data_length = len(self.viewer.layers['Shapes'].data[0])
        print(data_length)
        minX = 100000000
        maxX = 0
        minY = 100000000
        maxY = 0
        for i in range(0, data_length):
            print(self.viewer.layers['Shapes'].data[0][i])
            if minX > self.viewer.layers['Shapes'].data[0][i][1]:
                minX = self.viewer.layers['Shapes'].data[0][i][1]
            if maxX < self.viewer.layers['Shapes'].data[0][i][1]:
                maxX = self.viewer.layers['Shapes'].data[0][i][1]
            if minY > self.viewer.layers['Shapes'].data[0][i][2]:
                minY = self.viewer.layers['Shapes'].data[0][i][2]
            if maxY < self.viewer.layers['Shapes'].data[0][i][2]:
                maxY = self.viewer.layers['Shapes'].data[0][i][2]

        if self.aligned_1 is not None:
            self.crop_start_ratio_x = self.aligned_1.shape[1]/minX
            self.crop_size_ratio_x = self.aligned_1.shape[1]/(maxX-minX)
            self.crop_start_ratio_y = self.aligned_1.shape[2]/minY
            self.crop_size_ratio_y = self.aligned_1.shape[2]/(maxY-minY)
        elif self.aligned_2 is not None:
            self.crop_start_ratio_x = self.aligned_2.shape[1]/minX
            self.crop_size_ratio_x = self.aligned_2.shape[1]/(maxX-minX)
            self.crop_start_ratio_y = self.aligned_2.shape[2]/minY
            self.crop_size_ratio_y = self.aligned_2.shape[2]/(maxY-minY)
        elif self.aligned_3 is not None:
            self.crop_start_ratio_x = self.aligned_3.shape[1]/minX
            self.crop_size_ratio_x = self.aligned_3.shape[1]/(maxX-minX)
            self.crop_start_ratio_y = self.aligned_3.shape[2]/minY
            self.crop_size_ratio_y = self.aligned_3.shape[2]/(maxY-minY)
        elif self.aligned_4 is not None:
            self.crop_start_ratio_x = self.aligned_4.shape[1]/minX
            self.crop_size_ratio_x = self.aligned_4.shape[1]/(maxX-minX)
            self.crop_start_ratio_y = self.aligned_4.shape[2]/minY
            self.crop_size_ratio_y = self.aligned_4.shape[2]/(maxY-minY)

        self.viewer.layers.remove('Shapes')
        # print("crop to: {} {} {} {}".format(minX, maxX, minY, maxY))

        self.crop = True
        self.Load(text)
        self.crop = False

    def Load2D(self, text):

        print(self.search_folder.text() +
              str(self.comboBoxPath.currentText()) + '/mos.zarr')

        self.ds1 = xr.open_zarr(self.search_folder.text() +
                                str(self.comboBoxPath.currentText())+'/mos.zarr')
        ds2 = xr.open_zarr(self.search_folder.text(
        ) + str(self.comboBoxPath.currentText())+'/mos.zarr', group='l.2')
        ds4 = xr.open_zarr(self.search_folder.text(
        ) + str(self.comboBoxPath.currentText())+'/mos.zarr', group='l.4')
        ds8 = xr.open_zarr(self.search_folder.text(
        ) + str(self.comboBoxPath.currentText())+'/mos.zarr', group='l.8')
        ds16 = xr.open_zarr(self.search_folder.text(
        ) + str(self.comboBoxPath.currentText())+'/mos.zarr', group='l.16')
        ds32 = xr.open_zarr(self.search_folder.text(
        ) + str(self.comboBoxPath.currentText())+'/mos.zarr', group='l.32')

        
        z = self.scroll.value()

        if self.old_method:
            #Old method
            bscale = self.ds1.attrs['bscale']
            bzero = self.ds1.attrs['bzero']

            self.slice_names = self.ds1.attrs['cube_reg']['slice']

            self.scroll.setRange(0, len(self.slice_names)-1)
            if z >= len(self.slice_names):
                z = len(self.slice_names)-1
                self.scroll.setValue(z)

            slice_name = self.slice_names[z]
        else:
            bscale = self.ds1['S001'].attrs['bscale']
            bzero = self.ds1['S001'].attrs['bzero']
            slice_name = f"S{(z+1):03d}"




        # z = 0
        # slice_name = 'S00'+str(z+1)
        # if(z+1 > 9):
        #    slice_name = 'S0'+str(z+1)
        # if(z+1 > 99):
        #    slice_name = 'S'+str(z+1)


        if self.cb_C1.isChecked():
            im1 = (self.ds1[slice_name].sel(
                channel=1, type='mosaic', z=0).data * bscale + bzero)
            im2 = (ds2[slice_name].sel(
                channel=1, type='mosaic', z=0).data * bscale + bzero)
            im4 = (ds4[slice_name].sel(
                channel=1, type='mosaic', z=0).data * bscale + bzero)
            im8 = (ds8[slice_name].sel(
                channel=1, type='mosaic', z=0).data * bscale + bzero)
            im16 = (ds16[slice_name].sel(
                channel=1, type='mosaic', z=0).data * bscale + bzero)
            im32 = (ds32[slice_name].sel(
                channel=1, type='mosaic', z=0).data * bscale + bzero)

            with napari.gui_qt():
                self.viewer.add_image([im1, im2, im4, im8, im16, im32], multiscale=True,
                                      name='C1', blending='additive', colormap='bop purple', contrast_limits=[0,30])

        if self.cb_C2.isChecked():
            im1 = (self.ds1[slice_name].sel(
                channel=2, type='mosaic', z=0).data * bscale + bzero)
            im2 = (ds2[slice_name].sel(
                channel=2, type='mosaic', z=0).data * bscale + bzero)
            im4 = (ds4[slice_name].sel(
                channel=2, type='mosaic', z=0).data * bscale + bzero)
            im8 = (ds8[slice_name].sel(
                channel=2, type='mosaic', z=0).data * bscale + bzero)
            im16 = (ds16[slice_name].sel(
                channel=2, type='mosaic', z=0).data * bscale + bzero)
            im32 = (ds32[slice_name].sel(
                channel=2, type='mosaic', z=0).data * bscale + bzero)

            with napari.gui_qt():
                self.viewer.add_image([im1, im2, im4, im8, im16, im32], multiscale=True, 
                                        name='C2', blending='additive', colormap='red', contrast_limits=[0,30])
                

        if self.cb_C3.isChecked():
            im1 = (self.ds1[slice_name].sel(
                channel=3, type='mosaic', z=0).data * bscale + bzero)
            im2 = (ds2[slice_name].sel(
                channel=3, type='mosaic', z=0).data * bscale + bzero)
            im4 = (ds4[slice_name].sel(
                channel=3, type='mosaic', z=0).data * bscale + bzero)
            im8 = (ds8[slice_name].sel(
                channel=3, type='mosaic', z=0).data * bscale + bzero)
            im16 = (ds16[slice_name].sel(
                channel=3, type='mosaic', z=0).data * bscale + bzero)
            im32 = (ds32[slice_name].sel(
                channel=3, type='mosaic', z=0).data * bscale + bzero)

            with napari.gui_qt():
                self.viewer.add_image([im1, im2, im4, im8, im16, im32], multiscale=True,
                                      name='C3', blending='additive', colormap='green', contrast_limits=[0,30])
                #self.viewer.add_image(im4, multiscale=False, name='C3', blending='additive', colormap='green')

        if self.cb_C4.isChecked():
            im1 = (self.ds1[slice_name].sel(
                channel=4, type='mosaic', z=0).data * bscale + bzero)
            im2 = (ds2[slice_name].sel(
                channel=4, type='mosaic', z=0).data * bscale + bzero)
            im4 = (ds4[slice_name].sel(
                channel=4, type='mosaic', z=0).data * bscale + bzero)
            im8 = (ds8[slice_name].sel(
                channel=4, type='mosaic', z=0).data * bscale + bzero)
            im16 = (ds16[slice_name].sel(
                channel=4, type='mosaic', z=0).data * bscale + bzero)
            im32 = (ds32[slice_name].sel(
                channel=4, type='mosaic', z=0).data * bscale + bzero)

            with napari.gui_qt():
                self.viewer.add_image([im1, im2, im4, im8, im16, im32], multiscale=True,
                                      name='C4', blending='additive', colormap='bop blue', contrast_limits=[0,30])


    def Align(self, volume, resolution, output_resolution, start_slice_number):

        spacing = (self.ds1['S001'].attrs['scale'])
        # size_multiplier = (resolution*0.1*spacing[0])/7.5
        # print(resolution*0.1*spacing[0])
        size_multiplier = (resolution*0.1*spacing[0])/output_resolution
        size = (volume.shape[0], int(size_multiplier*volume.shape[1]), int(size_multiplier*volume.shape[2]))
        aligned = np.zeros(size, dtype=np.float32)
        size2D = (int(size_multiplier*volume.shape[2]), int(size_multiplier*volume.shape[1]))

        #print("spacing: {}".format(spacing))
        #print("volume.shape[0]: {}".format(volume.shape[0]))

        z_size = volume.shape[0]
        for z in range(0, z_size):

            fixed = sitk.GetImageFromArray(volume[z, :, :])
            fixed.SetOrigin((0, 0))

            # slice_name = 'S00'+str(z+1)
            # if(z+1 > 9):
            #    slice_name = 'S0'+str(z+1)
            # if(z+1 > 99):
            #    slice_name = 'S'+str(z+1)
            slice_name = self.slice_names[z+start_slice_number]

            spacing = (self.ds1[slice_name].attrs['scale'])
            # print(f"spacing {spacing}")
            # fixed.SetSpacing([1.6*spacing[0],1.6*spacing[1]])
            fixed.SetSpacing([resolution*0.1*spacing[1],resolution*0.1*spacing[0]])

            transform = sitk.Euler2DTransform()

            align_pos = z + start_slice_number #*self.optical_slices
            alignY = 0
            if not np.isnan(self.corrected_align_y[align_pos]):
                alignY = -self.corrected_align_y[align_pos]*0.1*spacing[1]

            alignX = 0
            if not np.isnan(self.corrected_align_x[align_pos]):
                alignX = -self.corrected_align_x[align_pos]*0.1*spacing[0]

            #print(f"alignX {alignX}")
            #print(f"align_pos {align_pos} alignY {self.corrected_align_y[align_pos]}")

            transform.SetTranslation([alignY, alignX])

            # transform.SetTranslation([-align_y[z]*0.1*spacing[1],-align_x[z]*0.1*spacing[0]])
            # transform.SetTranslation([0.5*-align_y[z],0.5*-align_x[z]])
            # transform.SetTranslation([0,0])
            
            #print('{}/{}, translate_x: {}, translate_y: {}'.format(z,z_size, alignX, alignY))
            #print('{}/{}, align_x[z]: {}, align_y[z]: {}'.format(z,z_size, self.align_x[z], self.align_y[z]))

            resampler = sitk.ResampleImageFilter()
            # resampler.SetReferenceImage(fixed)

            resampler.SetSize(size2D)
            resampler.SetOutputSpacing([output_resolution, output_resolution])
            resampler.SetOutputOrigin((0, 0))
            resampler.SetInterpolator(sitk.sitkLinear)
            resampler.SetDefaultPixelValue(0)
            resampler.SetTransform(transform)

            out = resampler.Execute(fixed)

            np_out = sitk.GetArrayFromImage(out)
            aligned[z, :, :] = np_out

        return aligned.astype(dtype=np.float32)

    def AlignNew(self, volume, resolution, output_resolution):
        
        if self.old_method:
            spacing = (self.ds1['S001'].attrs['scale'])
        else:
            spacing = [float,float,float]
            #print(self.ds1['S001'].attrs['scale'])
            #print(json.loads(self.ds1['S001'].attrs['scale']))
            #print(json.loads(self.ds1['S001'].attrs['scale'])["x"])
            #print(type(json.loads(self.ds1['S001'].attrs['scale'])["x"]))
            spacing[0] = float(json.loads(self.ds1['S001'].attrs['scale'])["x"])
            spacing[1] = float(json.loads(self.ds1['S001'].attrs['scale'])["y"])
            spacing[2] = float(json.loads(self.ds1['S001'].attrs['scale'])["z"])
            print(spacing[0])
            print(spacing[1])
            print(spacing[2])

        # size_multiplier = (resolution*0.1*spacing[0])/7.5
        # print(resolution*0.1*spacing[0])
        size_multiplier = (resolution*spacing[0])/output_resolution
        size = (volume.shape[0], int(size_multiplier*volume.shape[1]), int(size_multiplier*volume.shape[2]))
        aligned = np.zeros(size, dtype=np.float32)
        size2D = (int(size_multiplier*volume.shape[2]), int(size_multiplier*volume.shape[1]))

        #print("spacing: {}".format(spacing))
        #print("volume.shape[0]: {}".format(volume.shape[0]))

        z_size = volume.shape[0]
        for z in range(1, z_size+1):

            fixed = sitk.GetImageFromArray(volume[z-1, :, :])
            fixed.SetOrigin((0, 0))

            # slice_name = 'S00'+str(z+1)
            # if(z+1 > 9):
            #    slice_name = 'S0'+str(z+1)
            # if(z+1 > 99):
            #    slice_name = 'S'+str(z+1)
            if self.old_method:
                slice_name = self.slice_names[z-1]
                print(slice_name)
            else:
                slice_name = f"S{z:03d}"

                spacing[0] = float(json.loads(self.ds1[slice_name].attrs['scale'])["x"])
                spacing[1] = float(json.loads(self.ds1[slice_name].attrs['scale'])["y"])
                spacing[2] = float(json.loads(self.ds1[slice_name].attrs['scale'])["z"])

            #print(spacing)
            #spacing = (self.ds1[slice_name].attrs['scale'])
            #print(f"spacing {spacing}")
            # fixed.SetSpacing([1.6*spacing[0],1.6*spacing[1]])
            #print(type(resolution))
            #print(type(spacing[1]))
            
            fixed.SetSpacing([resolution*spacing[1],
                             resolution*spacing[0]])

            transform = sitk.Euler2DTransform()

            if self.old_method:
                alignY = 0
                if not np.isnan(self.align_y[z-1]):
                    alignY = -self.align_y[z-1]*0.1*spacing[1]

                alignX = 0
                if not np.isnan(self.align_x[z-1]):
                    alignX = -self.align_x[z-1]*0.1*spacing[0]
                print(alignY)
            else:
                alignX = self.ds1[slice_name].attrs['offsets']['x']*spacing[0]
                alignY = self.ds1[slice_name].attrs['offsets']['y']*spacing[1]
                print('{}/{}, offsets_x: {}, offsets_y: {}'.format(z, z_size, self.ds1[slice_name].attrs['offsets']['x'], self.ds1[slice_name].attrs['offsets']['y']))
                # print('{}/{}, spacing[0]: {}, spacing[1]: {}'.format(z, z_size, spacing[0], spacing[1]))
                
                #alignX = 0
                #alignY = 0

            transform.SetTranslation([alignY, alignX])

            # transform.SetTranslation([-align_y[z]*0.1*spacing[1],-align_x[z]*0.1*spacing[0]])
            # transform.SetTranslation([0.5*-align_y[z],0.5*-align_x[z]])
            # transform.SetTranslation([0,0])
            # print('{}/{}, alignX: {}, alignY: {}'.format(z, z_size, alignX, alignY))

            resampler = sitk.ResampleImageFilter()
            # resampler.SetReferenceImage(fixed)

            resampler.SetSize(size2D)
            resampler.SetOutputSpacing([output_resolution, output_resolution])
            resampler.SetOutputOrigin((0, 0))
            resampler.SetInterpolator(sitk.sitkLinear)
            resampler.SetDefaultPixelValue(0)
            resampler.SetTransform(transform)

            out = resampler.Execute(fixed)

            np_out = sitk.GetArrayFromImage(out)
            aligned[z-1, :, :] = np_out

        return aligned.astype(dtype=np.float32)

    def Remove_Regions(self, use_size):

        # output_resolution = float(self.pixel_size.text())
        threshold = float(self.thresholdN.text())

        if self.aligned_1 is not None and self.cb_C1.isChecked() and any(i.name == 'C1' for i in self.viewer.layers):

            threholded = self.aligned_1 > threshold
            threholded = ndimage.binary_fill_holes(threholded)
            threholded = threholded.astype(np.uint8)

            keep_n = self.spinN.value()

            for z in range(0, threholded.shape[0]):
                print('{}/{}'.format(z+1, threholded.shape[0]))

                threholded_z = threholded[z, :, :]
                volume_z = self.aligned_1[z, :, :]

                nb_components, output, stats, centroids = cv2.connectedComponentsWithStats(
                    threholded_z, connectivity=4)
                sizes = stats[:, -1]
                sizes_sorted = np.sort(sizes, axis=0)

                if use_size:
                    max_size = float(self.maxSizeN.text())
                else:
                    max_size = sizes_sorted[len(sizes_sorted)-1-keep_n]
                for i in range(1, nb_components):
                    if sizes[i] < max_size:
                        volume_z[output == i] = threshold

                self.aligned_1[z, :, :] = volume_z

            self.viewer.layers['C1'].visible = False
            self.viewer.layers['C1'].visible = True

        if self.aligned_2 is not None and self.cb_C2.isChecked() and any(i.name == 'C2' for i in self.viewer.layers):

            threholded = self.aligned_2 > threshold
            threholded = ndimage.binary_fill_holes(threholded)
            threholded = threholded.astype(np.uint8)

            # regions = threholded.copy()

            # with napari.gui_qt():
            #    self.viewer.add_image(threholded, name='C2 threholded', scale=(15,output_resolution,output_resolution), blending='additive', colormap='green')

            keep_n = self.spinN.value()

            for z in range(0, threholded.shape[0]):
                print('{}/{}'.format(z+1, threholded.shape[0]))

                threholded_z = threholded[z, :, :]
                volume_z = self.aligned_2[z, :, :]

                nb_components, output, stats, centroids = cv2.connectedComponentsWithStats(
                    threholded_z, connectivity=4)
                sizes = stats[:, -1]
                sizes_sorted = np.sort(sizes, axis=0)

                # regions[z,:,:] = output

                if use_size:
                    max_size = float(self.maxSizeN.text())
                else:
                    max_size = sizes_sorted[len(sizes_sorted)-1-keep_n]
                for i in range(1, nb_components):
                    if sizes[i] < max_size:
                        volume_z[output == i] = threshold

                self.aligned_2[z, :, :] = volume_z

            self.viewer.layers['C2'].visible = False
            self.viewer.layers['C2'].visible = True

            # with napari.gui_qt():
            #    self.viewer.add_image(regions, name='C2 regions', scale=(15,output_resolution,output_resolution), blending='additive', colormap='blue')

        if self.aligned_3 is not None and self.cb_C3.isChecked() and any(i.name == 'C3' for i in self.viewer.layers):

            threholded = self.aligned_3 > threshold
            threholded = ndimage.binary_fill_holes(threholded)
            threholded = threholded.astype(np.uint8)

            keep_n = self.spinN.value()

            for z in range(0, threholded.shape[0]):
                print('{}/{}'.format(z+1, threholded.shape[0]))

                threholded_z = threholded[z, :, :]
                volume_z = self.aligned_3[z, :, :]

                nb_components, output, stats, centroids = cv2.connectedComponentsWithStats(
                    threholded_z, connectivity=4)
                sizes = stats[:, -1]
                sizes_sorted = np.sort(sizes, axis=0)

                if use_size:
                    max_size = float(self.maxSizeN.text())
                else:
                    max_size = sizes_sorted[len(sizes_sorted)-1-keep_n]
                for i in range(1, nb_components):
                    if sizes[i] < max_size:
                        volume_z[output == i] = threshold

                self.aligned_3[z, :, :] = volume_z

            self.viewer.layers['C3'].visible = False
            self.viewer.layers['C3'].visible = True

        if self.aligned_4 is not None and self.cb_C4.isChecked() and any(i.name == 'C4' for i in self.viewer.layers):

            threholded = self.aligned_4 > threshold
            threholded = ndimage.binary_fill_holes(threholded)
            threholded = threholded.astype(np.uint8)

            keep_n = self.spinN.value()

            for z in range(0, threholded.shape[0]):
                print('{}/{}'.format(z+1, threholded.shape[0]))

                threholded_z = threholded[z, :, :]
                volume_z = self.aligned_4[z, :, :]

                nb_components, output, stats, centroids = cv2.connectedComponentsWithStats(
                    threholded_z, connectivity=4)
                sizes = stats[:, -1]
                sizes_sorted = np.sort(sizes, axis=0)

                if use_size:
                    max_size = float(self.maxSizeN.text())
                else:
                    max_size = sizes_sorted[len(sizes_sorted)-1-keep_n]
                for i in range(1, nb_components):
                    if sizes[i] < max_size:
                        volume_z[output == i] = threshold

                self.aligned_4[z, :, :] = volume_z

            self.viewer.layers['C4'].visible = False
            self.viewer.layers['C4'].visible = True

        # self.viewer.layers['C2'].contrast_limits=(0, 3)
        # with napari.gui_qt():
        #    if self.cb_C2.isChecked():
        #        self.viewer.add_image([self.aligned_2], name='C2 processed', scale=(15,output_resolution,output_resolution), blending='additive', colormap='red')

    def Remove_Small_Regions(self):
        self.Remove_Regions(True)

    def Keep_n_Regions(self):
        self.Remove_Regions(False)

    def Crop(self):

        output_resolution = float(self.pixel_size.text())
        threshold = float(self.thresholdN.text())

        maxX_final = 0
        minX_final = sys.maxsize
        maxY_final = 0
        minY_final = sys.maxsize
        maxZ_final = 0
        minZ_final = sys.maxsize

        if self.aligned_1 is not None and self.cb_C1.isChecked() and any(i.name == 'C1' for i in self.viewer.layers):

            tissue = self.aligned_1 > threshold
            objs = ndimage.find_objects(tissue)
            maxX = int(objs[0][0].stop)
            minX = int(objs[0][0].start)
            maxY = int(objs[0][1].stop)
            minY = int(objs[0][1].start)
            maxZ = int(objs[0][2].stop)
            minZ = int(objs[0][2].start)

            maxX_final = max(maxX_final, maxX)
            minX_final = min(minX_final, minX)
            maxY_final = max(maxY_final, maxY)
            minY_final = min(minY_final, minY)
            maxZ_final = max(maxZ_final, maxZ)
            minZ_final = min(minZ_final, minZ)

        if self.aligned_2 is not None and self.cb_C2.isChecked() and any(i.name == 'C2' for i in self.viewer.layers):

            print(self.aligned_2.shape)

            tissue = self.aligned_2 > threshold
            objs = ndimage.find_objects(tissue)
            maxX = int(objs[0][0].stop)
            minX = int(objs[0][0].start)
            maxY = int(objs[0][1].stop)
            minY = int(objs[0][1].start)
            maxZ = int(objs[0][2].stop)
            minZ = int(objs[0][2].start)

            maxX_final = max(maxX_final, maxX)
            minX_final = min(minX_final, minX)
            maxY_final = max(maxY_final, maxY)
            minY_final = min(minY_final, minY)
            maxZ_final = max(maxZ_final, maxZ)
            minZ_final = min(minZ_final, minZ)

        if self.aligned_3 is not None and self.cb_C3.isChecked() and any(i.name == 'C3' for i in self.viewer.layers):

            tissue = self.aligned_3 > threshold
            objs = ndimage.find_objects(tissue)
            maxX = int(objs[0][0].stop)
            minX = int(objs[0][0].start)
            maxY = int(objs[0][1].stop)
            minY = int(objs[0][1].start)
            maxZ = int(objs[0][2].stop)
            minZ = int(objs[0][2].start)

            maxX_final = max(maxX_final, maxX)
            minX_final = min(minX_final, minX)
            maxY_final = max(maxY_final, maxY)
            minY_final = min(minY_final, minY)
            maxZ_final = max(maxZ_final, maxZ)
            minZ_final = min(minZ_final, minZ)

        if self.aligned_4 is not None and self.cb_C4.isChecked() and any(i.name == 'C4' for i in self.viewer.layers):

            tissue = self.aligned_4 > threshold
            objs = ndimage.find_objects(tissue)
            maxX = int(objs[0][0].stop)
            minX = int(objs[0][0].start)
            maxY = int(objs[0][1].stop)
            minY = int(objs[0][1].start)
            maxZ = int(objs[0][2].stop)
            minZ = int(objs[0][2].start)

            maxX_final = max(maxX_final, maxX)
            minX_final = min(minX_final, minX)
            maxY_final = max(maxY_final, maxY)
            minY_final = min(minY_final, minY)
            maxZ_final = max(maxZ_final, maxZ)
            minZ_final = min(minZ_final, minZ)

        print("cropping to {}-{},{}-{}".format(minX_final,
              maxX_final, minY_final, maxY_final))

        if self.aligned_1 is not None and self.cb_C1.isChecked() and any(i.name == 'C1' for i in self.viewer.layers):
            self.aligned_1 = self.aligned_1[minX_final:maxX_final,
                                            minY_final:maxY_final, minZ_final:maxZ_final]
            self.viewer.layers.remove('C1')
            self.viewer.add_image([self.aligned_1], name='C1', scale=(
                15, output_resolution, output_resolution), blending='additive', colormap='bop purple', contrast_limits=[0,30])
        if self.aligned_2 is not None and self.cb_C2.isChecked() and any(i.name == 'C2' for i in self.viewer.layers):
            self.aligned_2 = self.aligned_2[minX_final:maxX_final,
                                            minY_final:maxY_final, minZ_final:maxZ_final]
            self.viewer.layers.remove('C2')
            self.viewer.add_image([self.aligned_2], name='C2', scale=(
                15, output_resolution, output_resolution), blending='additive', colormap='red', contrast_limits=[0,30])
        if self.aligned_3 is not None and self.cb_C3.isChecked() and any(i.name == 'C3' for i in self.viewer.layers):
            self.aligned_3 = self.aligned_3[minX_final:maxX_final,
                                            minY_final:maxY_final, minZ_final:maxZ_final]
            self.viewer.layers.remove('C3')
            self.viewer.add_image([self.aligned_3], name='C3', scale=(
                15, output_resolution, output_resolution), blending='additive', colormap='green', contrast_limits=[0,30])
        if self.aligned_4 is not None and self.cb_C4.isChecked() and any(i.name == 'C4' for i in self.viewer.layers):
            self.aligned_4 = self.aligned_4[minX_final:maxX_final,
                                            minY_final:maxY_final, minZ_final:maxZ_final]
            self.viewer.layers.remove('C4')
            self.viewer.add_image([self.aligned_4], name='C4', scale=(
                15, output_resolution, output_resolution), blending='additive', colormap='bop blue', contrast_limits=[0,30])

        self.viewer.layers.remove('Shapes')
        
    def CropToRegion(self):

        output_resolution = float(self.pixel_size.text())
        #print(self.viewer.layers['Shapes'])
        #print(self.viewer.layers['Shapes'].data)
        #print(self.viewer.layers.selection[0])
        data_length = len(self.viewer.layers['Shapes'].data[0])
        print(data_length)
        minX = 100000000
        maxX = 0
        minY = 100000000
        maxY = 0
        for i in range(0, data_length):
            print(self.viewer.layers['Shapes'].data[0][i])
            if minX > self.viewer.layers['Shapes'].data[0][i][1]:
                minX = self.viewer.layers['Shapes'].data[0][i][1]
            if maxX < self.viewer.layers['Shapes'].data[0][i][1]:
                maxX = self.viewer.layers['Shapes'].data[0][i][1]
            if minY > self.viewer.layers['Shapes'].data[0][i][2]:
                minY = self.viewer.layers['Shapes'].data[0][i][2]
            if maxY < self.viewer.layers['Shapes'].data[0][i][2]:
                maxY = self.viewer.layers['Shapes'].data[0][i][2]

        print("crop to: {} {} {} {}".format(minX, maxX, minY, maxY))

        if self.aligned_1 is not None and self.cb_C1.isChecked() and any(i.name == 'C1' for i in self.viewer.layers):
            self.aligned_1 = self.aligned_1[:, int(
                minX):int(maxX), int(minY):int(maxY)]
            self.viewer.layers.remove('C1')
            self.image_translation = (int(minX*output_resolution), int(minY*output_resolution))
            self.viewer.add_image([self.aligned_1], name='C1', scale=(float(self.slice_spacing)/float(self.optical_slices), output_resolution, output_resolution), blending='additive', colormap='bop purple')#, translate=(0, int(minX*output_resolution), int(minY*output_resolution)))
        if self.aligned_2 is not None and self.cb_C2.isChecked() and any(i.name == 'C2' for i in self.viewer.layers):
            self.aligned_2 = self.aligned_2[:, int(
                minX):int(maxX), int(minY):int(maxY)]
            self.viewer.layers.remove('C2')
            self.image_translation = (int(minX*output_resolution), int(minY*output_resolution))
            self.viewer.add_image([self.aligned_2], name='C2', scale=(float(self.slice_spacing)/float(self.optical_slices), output_resolution, output_resolution), blending='additive', colormap='red')#, translate=(0, int(minX*output_resolution), int(minY*output_resolution)))
        if self.aligned_3 is not None and self.cb_C3.isChecked() and any(i.name == 'C3' for i in self.viewer.layers):
            self.aligned_3 = self.aligned_3[:, int(
                minX):int(maxX), int(minY):int(maxY)]
            self.viewer.layers.remove('C3')
            self.image_translation = (int(minX*output_resolution), int(minY*output_resolution))
            self.viewer.add_image([self.aligned_3], name='C3', scale=(float(self.slice_spacing)/float(self.optical_slices), output_resolution, output_resolution), blending='additive', colormap='green')#, translate=(0, int(minX*output_resolution), int(minY*output_resolution)))
        if self.aligned_4 is not None and self.cb_C4.isChecked() and any(i.name == 'C4' for i in self.viewer.layers):
            self.aligned_4 = self.aligned_4[:, int(
                minX):int(maxX), int(minY):int(maxY)]
            self.viewer.layers.remove('C4')
            self.image_translation = (int(minX*output_resolution), int(minY*output_resolution))
            self.viewer.add_image([self.aligned_4], name='C4', scale=(float(self.slice_spacing)/float(self.optical_slices), output_resolution, output_resolution), blending='additive', colormap='bop blue')#, translate=(0, int(minX*output_resolution), int(minY*output_resolution)))
            
            
        self.viewer.layers.remove('Shapes')

                
    def add_polygon_simple(self):

        pos = int(self.viewer.cursor.position[0]/15)

        data_length = len(self.viewer.layers['Shapes'].data[0])
        print(data_length)

        new_shape = []
        for i in range(0, data_length):
            
            x = pos
            y = self.viewer.layers['Shapes'].data[0][i][1]
            z = self.viewer.layers['Shapes'].data[0][i][2]
            new_shape.append((x,y,z))

        new_shape = np.array(new_shape)
        
        shapes_layer = self.viewer.add_shapes(new_shape, shape_type='polygon', name = "Shapes", scale=(15, 15, 15),)
    
    
    def add_polygon(self):
        self.viewer.window.qt_viewer.update()

        #print("pos {}".format(self.viewer.dims.point[0]))


        

        #pos = float(self.viewer.cursor.position[0]/(float(self.slice_spacing)/float(self.optical_slices)))
        pos = self.viewer.dims.point[0] / (float(self.slice_spacing)/float(self.optical_slices))
        print("pos {}".format(pos))
        
        data_length = len(self.viewer.layers['Shapes'].data[0])
        print("data_length {}".format(data_length))

        contour_list = []
        z_pos = self.viewer.layers["Shapes"].data[0][0][0] # z pos
        contour = self.viewer.layers["Shapes"].data[0]
        contour_list.append([z_pos,contour])

        i = 1
        while True:
            layer_name = "Shapes [{}]".format(i)
            try:
                z_pos = self.viewer.layers[layer_name].data[0][0][0] # z pos
                contour = self.viewer.layers[layer_name].data[0]
                contour_list.append((z_pos,contour))
                i = i+1
            except:
                break

        contour_list_sorted = sorted(contour_list, key=lambda tup: tup[0])





        print("contour_list_sorted[0][0] {}".format(contour_list_sorted[0][0]))
        
            
        new_shape = []

        if pos < contour_list_sorted[0][0]:
            for i in range(0, data_length):
                x = pos
                y = float(contour_list_sorted[0][1][i][1])
                z = float(contour_list_sorted[0][1][i][2])

                new_shape.append((x,y,z))

        elif pos > contour_list_sorted[len(contour_list_sorted)-1][0]:
            for i in range(0, data_length):
                x = pos
                y = float(contour_list_sorted[len(contour_list_sorted)-1][1][i][1])
                z = float(contour_list_sorted[len(contour_list_sorted)-1][1][i][2])

                new_shape.append((x,y,z))
            
        else:
            z_i = 0
            for i in range(0,len(contour_list_sorted)-1):
                z_level_start = contour_list_sorted[i][0]
                z_level_end = contour_list_sorted[i+1][0]
                if pos >= z_level_start and pos <= z_level_end:
                    z_i = i

            for i in range(0, data_length):
                x1 = float(contour_list_sorted[z_i][1][i][0])
                y1 = float(contour_list_sorted[z_i][1][i][1])
                z1 = float(contour_list_sorted[z_i][1][i][2])

                x2 = float(contour_list_sorted[z_i+1][1][i][0])
                y2 = float(contour_list_sorted[z_i+1][1][i][1])
                z2 = float(contour_list_sorted[z_i+1][1][i][2])

                weight2 = (pos - x1) / (x2 - x1)
                weight = 1 - weight2

                x = pos
                y = (weight * y1) + (weight2 * y2)
                z = (weight * z1) + (weight2 * z2)
                new_shape.append((x,y,z))


        new_shape = np.array(new_shape)
        output_resolution = float(self.pixel_size.text())
        shapes_layer = self.viewer.add_shapes(new_shape, shape_type='polygon', name = "Shapes", scale=(float(self.slice_spacing)/float(self.optical_slices), output_resolution, output_resolution),)
    




    def run_remove_outside(self):

        

        
        if self.aligned_1 is not None:
            aligned1_tmp = np.copy(self.aligned_1)
        if self.aligned_2 is not None:
            aligned2_tmp = np.copy(self.aligned_2)
        if self.aligned_3 is not None:
            aligned3_tmp = np.copy(self.aligned_3)
        if self.aligned_4 is not None:
            aligned4_tmp = np.copy(self.aligned_4)

        contour_list = []
        z_pos = self.viewer.layers["Shapes"].data[0][0][0] # z pos
        contour = self.viewer.layers["Shapes"].data[0]
        contour_list.append([z_pos,contour])


        i = 1
        while True:
            layer_name = "Shapes [{}]".format(i)
            try:
                z_pos = self.viewer.layers[layer_name].data[0][0][0] # z pos
                contour = self.viewer.layers[layer_name].data[0]
                contour_list.append((z_pos,contour))
                i = i+1
            except:
                break
        contour_list_sorted = sorted(contour_list, key=lambda tup: tup[0])
        

        

        output_resolution = float(self.pixel_size.text())
        data_length = len(self.viewer.layers['Shapes'].data[0])
        print("data_length {}".format(data_length))

        c = 0
        width = 0
        height = 0
        if self.aligned_1 is not None:
            c, width, height = self.aligned_1.shape
        if self.aligned_2 is not None:
            c, width, height = self.aligned_2.shape
        if self.aligned_3 is not None:
            c, width, height = self.aligned_3.shape
        if self.aligned_4 is not None:
            c, width, height = self.aligned_4.shape

        for z_level in range(0,c):

            z_i = -1
            for i in range(0,len(contour_list_sorted)-1):
                z_level_start = contour_list_sorted[i][0]
                z_level_end = contour_list_sorted[i+1][0]
                if z_level >= z_level_start and z_level <= z_level_end:
                    z_i = i

            if z_i == -1:
                mask = np.zeros((width, height), dtype=int)
            else:
                polygon_values = []
                for i in range(0, data_length):
                    x1 = float(contour_list_sorted[z_i][1][i][0])
                    y1 = float(contour_list_sorted[z_i][1][i][1])
                    z1 = float(contour_list_sorted[z_i][1][i][2])

                    x2 = float(contour_list_sorted[z_i+1][1][i][0])
                    y2 = float(contour_list_sorted[z_i+1][1][i][1])
                    z2 = float(contour_list_sorted[z_i+1][1][i][2])

                    weight2 = (z_level - x1) / (x2 - x1)
                    weight = 1 - weight2

                    y = (weight * y1) + (weight2 * y2)# + (self.image_translation[0]/15)
                    z = (weight * z1) + (weight2 * z2)#  + (self.image_translation[1]/15)
                    polygon_values.append((y, z))


                img = Image.new('L', (width, height), 0)
                ImageDraw.Draw(img).polygon(polygon_values, outline=1, fill=1)

                mask = np.array(img)
                mask = np.transpose(mask, (1,0))

                #self.viewer.add_image([mask], name='mask', scale=(output_resolution, output_resolution), blending='additive', colormap='bop purple')
                #aligned3_tmp[z_level,:,:] = aligned3_tmp[z_level,:,:] * mask

            
            if self.aligned_1 is not None:
                aligned1_tmp[z_level,:,:] = aligned1_tmp[z_level,:,:] * mask
            if self.aligned_2 is not None:
                aligned2_tmp[z_level,:,:] = aligned2_tmp[z_level,:,:] * mask
            if self.aligned_3 is not None:
                aligned3_tmp[z_level,:,:] = aligned3_tmp[z_level,:,:] * mask
            if self.aligned_4 is not None:
                aligned4_tmp[z_level,:,:] = aligned4_tmp[z_level,:,:] * mask

        
        #self.viewer.layers.remove('C3')
           


        if self.aligned_1 is not None:
            self.viewer.add_image([aligned1_tmp], name='C1_masked', scale=(15, output_resolution, output_resolution), blending='additive', colormap='bop purple', translate=(0, 0, 0))
        if self.aligned_2 is not None:
            self.viewer.add_image([aligned2_tmp], name='C2_masked', scale=(15, output_resolution, output_resolution), blending='additive', colormap='red', translate=(0, 0, 0))
        if self.aligned_3 is not None:
            self.viewer.add_image([aligned3_tmp], name='C3_masked', scale=(15, output_resolution, output_resolution), blending='additive', colormap='green', translate=(0, 0, 0))
        if self.aligned_4 is not None:
            self.viewer.add_image([aligned4_tmp], name='C4_masked', scale=(15, output_resolution, output_resolution), blending='additive', colormap='bop blue', translate=(0, 0, 0))


        #self.viewer.layers.remove('Shapes')


    def SelectFolder(self):
        file = str(QtWidgets.QFileDialog.getExistingDirectory())
        file = file + '/'
        # text_files = glob.glob(path + "/**/**/*.zarr", recursive = False)
        # print(text_files)
        # return

        self.search_folder.setText(file)
        self.comboBoxPath.clear()
        for f in os.scandir(self.search_folder.text()):
            if f.is_dir():
                if os.path.exists(f.path + "/mos.zarr"):
                    s = f.path
                    s = s.replace(self.search_folder.text(), '')
                    print(s)
                    self.comboBoxPath.addItem(s)

    def on_combobox_changed(self):
        if os.path.exists(self.search_folder.text() + str(self.comboBoxPath.currentText())+'/mos.zarr/.zmetadata'):
            print("metadata is available")
            ds = xr.open_zarr(self.search_folder.text(
            ) + str(self.comboBoxPath.currentText())+'/mos.zarr', consolidated=True)
            # print(ds.attrs["cube_reg"]['slice'])
            length = len(ds.attrs["multiscale"]['datasets'])
            
            #self.comboBoxResolution.clear()
            #for n in range(0, length):
            #    print("adding resolution {}".format(
            #        ds.attrs["multiscale"]['datasets'][n]['level']))
            #    self.comboBoxResolution.addItem(
            #        str(ds.attrs["multiscale"]['datasets'][n]['level']))
            #    self.comboBoxResolution.setCurrentText(
            #        str(ds.attrs["multiscale"]['datasets'][n]['level']))
            try:
                self.scroll.setValue(0)
                self.image_slice.setText("0")
                self.slice_names = ds.attrs['cube_reg']['slice']
                self.scroll.setRange(0, len(self.slice_names))
            except Exception:
                # print("not initialised")
                pass
        else:
            print("adding default resolution")
            #self.comboBoxResolution.clear()
            #self.comboBoxResolution.addItem(str(1))
            #self.comboBoxResolution.addItem(str(2))
            #self.comboBoxResolution.addItem(str(4))
            #self.comboBoxResolution.addItem(str(8))
            #self.comboBoxResolution.addItem(str(16))
            #self.comboBoxResolution.addItem(str(32))
            #self.comboBoxResolution.setCurrentText(str(32))

    def set_image_slice_value(self):
        value = self.scroll.value()
        self.image_slice.setText(str(value))

    def run_stardist(self):
        print("run_stardist")

        print(self.viewer.layers['Shapes'].data)
        # print()

        data_length = len(self.viewer.layers['Shapes'].data[0])
        print(data_length)
        minX = 100000000
        maxX = 0
        minY = 100000000
        maxY = 0
        for i in range(0, data_length):
            print(self.viewer.layers['Shapes'].data[0][i])
            if minX > self.viewer.layers['Shapes'].data[0][i][0]:
                minX = self.viewer.layers['Shapes'].data[0][i][0]
            if maxX < self.viewer.layers['Shapes'].data[0][i][0]:
                maxX = self.viewer.layers['Shapes'].data[0][i][0]
            if minY > self.viewer.layers['Shapes'].data[0][i][1]:
                minY = self.viewer.layers['Shapes'].data[0][i][1]
            if maxY < self.viewer.layers['Shapes'].data[0][i][1]:
                maxY = self.viewer.layers['Shapes'].data[0][i][1]

        minX2= minX / 2
        #maxX = maxX / 2
        minY2 = minY / 2
        #maxY = maxY / 2
                
        maxX2 = (minX + (maxX - minX)/2)/2
        maxY2 = (minY + (maxY - minY)/2)/2
                
        print("{} {} {} {}".format(minX, maxX, minY, maxY))
        # return
        # prints a list of available models
        StarDist2D.from_pretrained()

        # creates a pretrained model
        model = StarDist2D.from_pretrained('2D_versatile_fluo')

        # import image
        # image = skimage.io.imread('H:/GitRepositories/STPT/crop.tif')

        bscale = self.ds1.attrs['bscale']
        bzero = self.ds1.attrs['bzero']

        z = self.scroll.value()

        self.slice_names = self.ds1.attrs['cube_reg']['slice']

        self.scroll.setRange(0, len(self.slice_names)-1)
        if z >= len(self.slice_names):
            z = len(self.slice_names)-1
            self.scroll.setValue(z)

        slice_name = self.slice_names[z]

        # print(image.shape)
        # print(type(image))
        # return
        # normalize image
        #ds1 = xr.open_zarr(self.search_folder.text(
        #) + str(self.comboBoxPath.currentText())+'/mos.zarr', group='l.1')
        
        im12 = self.ds1[slice_name].sel(
            channel=2, type='mosaic', z=0).data * bscale + bzero
        im12 = np.array(im12[int(minX):int(maxX), int(minY):int(maxY)])
        im13 = self.ds1[slice_name].sel(
            channel=3, type='mosaic', z=0).data * bscale + bzero
        im13 = np.array(im13[int(minX):int(maxX), int(minY):int(maxY)])

        #im12 = np.array(self.aligned_2)
        #im13 = np.array(self.aligned_3)
        
        
        # img = np.array(im1[(im1.shape[0]/2)-1000:(im1.shape[0]/2)+1000,(im1.shape[1]/2)-1000:(im1.shape[1]/2)+1000,])

        img = im12 - ((0.128 * im13) + 17.371)

        print(type(img))
        print(img.shape)

        # img = np.array(im1[(im1.shape[0]/2)-1000:(im1.shape[0]/2)+1000,(im1.shape[1]/2)-1000:(im1.shape[1]/2)+1000,])
        # print(type(img))
        # print(img.shape)

        img = normalize(img, 1, 99.8, axis=(0, 1))
        # return
        # predict segmentation
        labels, details = model.predict_instances(img, n_tiles=(10, 10))
        print(type(details['points'].shape[0]))
        print(details['points'].shape)
        #print(details)
        #print(labels.shape)
        #print(type(labels))
        
        max_value = np.max(labels)
        print(f'max_value {max_value}')
        
        refined_labels = labels.copy()
        minimum_cc_sum = 100000
        for label in range(max_value+1):
            if np.sum(refined_labels[labels == label]) < minimum_cc_sum:
                refined_labels[labels == label] = 0
                
                
        refined_labels2 = labels.copy()
        minimum_cc_sum = 100000
        for label in range(max_value+1):
            if np.sum(refined_labels2[labels == label]) > minimum_cc_sum:
                refined_labels2[labels == label] = 0

        # with napari.gui_qt():
        self.viewer.add_image(
            img, name='img', blending='additive', scale=(1,1), translate=(int(minX), int(minY)))

        self.viewer.add_image(
            im12, name='im12', blending='additive', scale=(1,1), translate=(int(minX), int(minY)))
        self.viewer.add_image(
            im13, name='im13', blending='additive', scale=(1,1), translate=(int(minX), int(minY)))
        self.viewer.add_labels(
            labels, name='labels', blending='additive', scale=(1,1), translate=(int(minX), int(minY)))
        self.viewer.add_labels(
            refined_labels, name='refined_labels', blending='additive', scale=(1,1), translate=(int(minX), int(minY)))
        self.viewer.add_labels(
            refined_labels2, name='refined_labels2', blending='additive', scale=(1,1), translate=(int(minX), int(minY)))
        self.viewer.add_points(
            details['points'], translate=(int(minX), int(minY)))

    #def SaveMovie(self):

        #widget = QtWidgets.QWidget()
        #fname = QtWidgets.QFileDialog.getSaveFileName(
        #    widget, 'Save file as', 'c:\\', "Image files (*.gif *.mp4)")

        #if isinstance(fname, tuple):
        #    print(fname[0])
        #    filename, file_extension = os.path.splitext(fname[0])
        #    if file_extension == '.mp4':
        #        self.movie.make_movie(name='movie.mp4', resolution=300, fps=20)
        #    if file_extension == '.gif':
        #        self.movie.make_gif('S:/Tristan/gifmovie.gif')

    def SaveVolume(self):
        widget = QtWidgets.QWidget()
        if sys.platform == 'linux':
            fname, _ = QtWidgets.QFileDialog.getSaveFileName(
                widget, 'Save file as', '/home/', "Image files (*.tiff *.mha)")
        else:
            fname, _ = QtWidgets.QFileDialog.getSaveFileName(
                widget, 'Save file as', 'c:\\', "Image files (*.tiff *.mha)")

        if fname != "":
            print(fname)

            file_name = fname#"/home/tristan/test.tiff"

            output_resolution = float(self.pixel_size.text())
            spacing = (output_resolution, output_resolution,
                       15/self.optical_slices)

            # Not the most efficinet method
            #volume = np.zeros(self.shape)
            #volume = np.expand_dims(volume, axis=3)
            volume = []

            #print(f"volume {volume.shape}")
            if self.cb_C1.isChecked():
                #volume = np.concatenate((volume, np.expand_dims(self.aligned_1, axis=3)), axis=3)
                volume.append(self.aligned_1)
            #except Exception:
            #   pass
            if self.cb_C2.isChecked():
                #volume = np.concatenate((volume, np.expand_dims(self.aligned_2, axis=3)), axis=3)
                volume.append(self.aligned_2)
            #except Exception:
                pass
            if self.cb_C3.isChecked():
                #volume = np.concatenate((volume, np.expand_dims(self.aligned_3, axis=3)), axis=3)
                volume.append(self.aligned_3)
            #except Exception:
                pass
            if self.cb_C4.isChecked():
                #volume = np.concatenate((volume, np.expand_dims(self.aligned_4, axis=3)), axis=3)
                volume.append(self.aligned_4)
            #except Exception:
            #    pass

            volume = np.array(volume)
            print(f"volume {volume.shape}")
            volume = np.moveaxis(volume,0,3)
            print(f"volume {volume.shape}")

            #volume = np.delete(volume, 0, 3)
            

            #print(f"volume {volume.shape}")
            #100 774 851 4

            volume_itk = sitk.GetImageFromArray(volume)
            volume_itk.SetSpacing(spacing)
            caster = sitk.CastImageFilter()
            caster.SetOutputPixelType(sitk.sitkVectorFloat32)
            volume_itk = caster.Execute(volume_itk)
            sitk.WriteImage(volume_itk, file_name)

    def SaveSlice(self):
        widget = QtWidgets.QWidget()
        fname = QtWidgets.QFileDialog.getSaveFileName(
            widget, 'Save file as', 'c:\\', "Image files (*.tiff *.mha)")

        if isinstance(fname, tuple):
            print(fname[0])

            bscale = self.ds1.attrs['bscale']
            bzero = self.ds1.attrs['bzero']

            z = self.scroll.value()

            self.slice_names = self.ds1.attrs['cube_reg']['slice']

            self.scroll.setRange(0, len(self.slice_names)-1)
            if z >= len(self.slice_names):
                z = len(self.slice_names)-1
                self.scroll.setValue(z)

            slice_name = self.slice_names[z]



            output_resolution = float(self.pixel_size.text())
            spacing = (output_resolution, output_resolution)

            #im1 = (self.ds1[slice_name].sel(
            #    channel=3, type='mosaic', z=0).data * bscale + bzero).astype(dtype=np.float32)

                
            ds32 = xr.open_zarr(self.search_folder.text(
                ) + str(self.comboBoxPath.currentText())+'/mos.zarr', group='l.32')
            im32 = (ds32[slice_name].sel(
                channel=3, type='mosaic', z=0).data * bscale + bzero).astype(dtype=np.float32)

            image_itk = sitk.GetImageFromArray(im32)
            image_itk.SetSpacing(spacing)
            caster = sitk.CastImageFilter()
            caster.SetOutputPixelType(sitk.sitkVectorFloat32)
            image_itk = caster.Execute(image_itk)
            sitk.WriteImage(image_itk, fname[0])


    def Normalize_slices(self, volume, optical_sections):
        import statistics
        print("Normalize optical sections")

        #output_resolution = float(self.pixel_size.text())

        #self.aligned_3 = self.aligned_3[4::5,:,:]

        #with napari.gui_qt():
        #    self.viewer.layers.remove('C3')
        #    self.viewer.add_image([self.aligned_3], name='C3', scale=(float(self.slice_spacing)/float(self.optical_slices), output_resolution, output_resolution), blending='additive', colormap='green', contrast_limits=[0,30])#, translate=(0, int(start_y*output_resolution), int(start_z*output_resolution)))

        #return 


        #output_resolution = float(self.pixel_size.text())

        #if self.cb_C2.isChecked():


        slices, size_x, size_y = volume.shape
        for i in range(1,slices):
            if i%optical_sections != 0:
                values_x = []
                values_y = []
                for j in range(1000):
                    rand_x = random.randint(1,size_x-1)
                    rand_y = random.randint(1,size_y-1)
                    if(volume[i-1,rand_x,rand_y] > 0.0 and volume[i,rand_x,rand_y] > 0.0):
                        x = volume[i,rand_x,rand_y]
                        y = volume[i-1,rand_x,rand_y]
                        values_x.append(x)
                        values_y.append(y)

                slope, intercept, r, p, std_err = stats.linregress(values_x, values_y)

                volume[i,:,:] = volume[i,:,:] * slope + intercept

        return




        slices, size_x, size_y = volume.shape
        for i in range(1,slices):
            if i%5 != 0:
                values = []
                for j in range(10000):
                    rand_x = random.randint(1,size_x-1)
                    rand_y = random.randint(1,size_y-1)
                    if(volume[i-1,rand_x,rand_y] > 0.1 and volume[i,rand_x,rand_y] > 0.1 and volume[i-1,rand_x,rand_y] < 1 and volume[i,rand_x,rand_y] < 1):
                        division = volume[i-1,rand_x,rand_y] / volume[i,rand_x,rand_y]
                        values.append(division)

                mean = sum(values) / len(values)
                volume[i,:,:] = volume[i,:,:] * mean

        return

        if self.cb_C3.isChecked():  
            slices, size_x, size_y = self.aligned_3.shape
            for i in range(1,slices):
                if i%5 != 0:
                    values = []
                    for j in range(1000):
                        rand_x = random.randint(1,size_x-1)
                        rand_y = random.randint(1,size_y-1)
                        if(self.aligned_3[i-1,rand_x,rand_y] > 0.1 and self.aligned_3[i,rand_x,rand_y] > 0.1):
                            division = self.aligned_3[i-1,rand_x,rand_y] / self.aligned_3[i,rand_x,rand_y]
                            values.append(division)

                    mean = sum(values) / len(values)
                    self.aligned_3[i,:,:] = self.aligned_3[i,:,:] * mean

            with napari.gui_qt():
                self.viewer.layers.remove('C3')
                self.viewer.add_image([self.aligned_3], name='C3', scale=(float(self.slice_spacing)/float(self.optical_slices), output_resolution, output_resolution), blending='additive', colormap='green', contrast_limits=[0,30])#, translate=(0, int(start_y*output_resolution), int(start_z*output_resolution)))



        return



    def Normalize(self):

        output_resolution = float(self.pixel_size.text())
        norm_value = float(self.normalize_value.text())

        if self.cb_C1.isChecked():
            self.aligned_1[0::10,:,:] = self.aligned_1[0::10,:,:] * norm_value
            self.aligned_1[1::10,:,:] = self.aligned_1[1::10,:,:] * norm_value
            self.aligned_1[2::10,:,:] = self.aligned_1[2::10,:,:] * norm_value
            self.aligned_1[3::10,:,:] = self.aligned_1[3::10,:,:] * norm_value
            self.aligned_1[4::10,:,:] = self.aligned_1[4::10,:,:] * norm_value

            with napari.gui_qt():
                self.viewer.layers.remove('C1')
                self.viewer.add_image([self.aligned_1], name='C1', scale=(float(self.slice_spacing)/float(self.optical_slices), output_resolution, output_resolution), blending='additive', colormap='bop purple', contrast_limits=[0,30])#, translate=(0, int(start_y*output_resolution), int(start_z*output_resolution)))


        if self.cb_C2.isChecked():  
            self.aligned_2[0::10,:,:] = self.aligned_2[0::10,:,:] * norm_value
            self.aligned_2[1::10,:,:] = self.aligned_2[1::10,:,:] * norm_value
            self.aligned_2[2::10,:,:] = self.aligned_2[2::10,:,:] * norm_value
            self.aligned_2[3::10,:,:] = self.aligned_2[3::10,:,:] * norm_value
            self.aligned_2[4::10,:,:] = self.aligned_2[4::10,:,:] * norm_value

            with napari.gui_qt():
                self.viewer.layers.remove('C2')
                self.viewer.add_image([self.aligned_2], name='C2', scale=(float(self.slice_spacing)/float(self.optical_slices), output_resolution, output_resolution), blending='additive', colormap='red', contrast_limits=[0,30])#, translate=(0, int(start_y*output_resolution), int(start_z*output_resolution)))

        if self.cb_C3.isChecked():

            slices, size_x, size_y = self.aligned_3.shape
            for i in range(1,slices):
                if i%2 != 0:
                    values_x = []
                    values_y = []
                    for j in range(1000):
                        rand_x = random.randint(1,size_x-1)
                        rand_y = random.randint(1,size_y-1)
                        if(self.aligned_3[i-1,rand_x,rand_y] > 0.0 and self.aligned_3[i,rand_x,rand_y] > 0.0):
                            x = self.aligned_3[i,rand_x,rand_y]
                            y = self.aligned_3[i-1,rand_x,rand_y]
                            values_x.append(x)
                            values_y.append(y)

                    slope, intercept, r, p, std_err = stats.linregress(values_x, values_y)

                    print(slope)

                    self.aligned_3[i,:,:] = self.aligned_3[i,:,:] * slope + intercept

            with napari.gui_qt():
                self.viewer.layers.remove('C3')
                self.viewer.add_image([self.aligned_3], name='C3', scale=(float(self.slice_spacing)/float(self.optical_slices), output_resolution, output_resolution), blending='additive', colormap='green', contrast_limits=[0,30])#, translate=(0, int(start_y*output_resolution), int(start_z*output_resolution)))

        if self.cb_C4.isChecked():
            self.aligned_4[0::10,:,:] = self.aligned_4[0::10,:,:] * norm_value
            self.aligned_4[1::10,:,:] = self.aligned_4[1::10,:,:] * norm_value
            self.aligned_4[2::10,:,:] = self.aligned_4[2::10,:,:] * norm_value
            self.aligned_4[3::10,:,:] = self.aligned_4[3::10,:,:] * norm_value
            self.aligned_4[4::10,:,:] = self.aligned_4[4::10,:,:] * norm_value

            with napari.gui_qt():
                self.viewer.layers.remove('C4')
                self.viewer.add_image([self.aligned_4], name='C4', scale=(float(self.slice_spacing)/float(self.optical_slices), output_resolution, output_resolution), blending='additive', colormap='bop blue', contrast_limits=[0,30])#, translate=(0, int(start_y*output_resolution), int(start_z*output_resolution)))

            

    def main(self):

        with napari.gui_qt():

            self.viewer = napari.Viewer()

            widget = QtWidgets.QWidget()
            vbox = QtWidgets.QVBoxLayout()

            hbox = QtWidgets.QHBoxLayout()
            hbox.addWidget(QtWidgets.QLabel("Data folder:"))
            
            print(sys.platform)
            if sys.platform == 'linux':
                if self.old_method:
                    self.search_folder = QtWidgets.QLineEdit('/data/meds1_c/storage/processed/stpt')
                else:
                    self.search_folder = QtWidgets.QLineEdit('/data/meds1_a/eglez/stpt_out')
            else:
                self.search_folder = QtWidgets.QLineEdit('N:/stpt/')
            
            self.search_folder.setMaximumWidth(250)
            hbox.addWidget(self.search_folder)
            bSelectFolder = QtWidgets.QPushButton('Set folder')
            bSelectFolder.setCheckable(True)
            bSelectFolder.clicked.connect(self.SelectFolder)
            hbox.addWidget(bSelectFolder)
            vbox.addLayout(hbox)

            hbox = QtWidgets.QHBoxLayout()
            self.comboBoxPath = QtWidgets.QComboBox()

            #text_file = open("sample.txt", "wt")

            self.comboBoxPath.clear()
            if os.path.exists(self.search_folder.text()):
                for f in os.scandir(self.search_folder.text()):
                    if f.is_dir():
                        if os.path.exists(f.path + "/mos.zarr"):
                            s = f.path
                            s = s.replace(self.search_folder.text(), '')
                            print(s)
                            # print(type(s))
                            # text_file.write(s)
                            # text_file.write('\n')
                            self.comboBoxPath.addItem(s)

            # text_file.close()

            self.comboBoxPath.setMaximumWidth(300)
            self.comboBoxPath.currentIndexChanged.connect(
                self.on_combobox_changed)
            hbox.addWidget(self.comboBoxPath)
            hbox.addStretch(1)
            vbox.addLayout(hbox)

            hbox = QtWidgets.QHBoxLayout()
            #hbox.addWidget(QtWidgets.QLabel("Load level:"))
            #self.comboBoxResolution = QtWidgets.QComboBox()
            #self.comboBoxResolution.setMaximumWidth(100)
            #self.on_combobox_changed()
            #hbox.addWidget(self.comboBoxResolution)
            #hbox.addStretch(1)

            hbox.addWidget(QtWidgets.QLabel("Slice spacing:"))
            self.m_slice_spacing = QtWidgets.QLineEdit("15")
            self.m_slice_spacing.setMaximumWidth(50)
            hbox.addWidget(self.m_slice_spacing)

            hbox.addWidget(QtWidgets.QLabel("Output pixel size:"))
            self.pixel_size = QtWidgets.QLineEdit("15")
            self.pixel_size.setMaximumWidth(50)
            hbox.addWidget(self.pixel_size)
            hbox.addStretch(1)
            vbox.addLayout(hbox)


            hbox = QtWidgets.QHBoxLayout()
            self.cb_correct_brightness_optical_section = QtWidgets.QCheckBox('Normalize brightness optical sections')
            self.cb_correct_brightness_optical_section.setChecked(True)
            hbox.addWidget(self.cb_correct_brightness_optical_section)
            hbox.addStretch(1)
            vbox.addLayout(hbox)
            

            hbox = QtWidgets.QHBoxLayout()
            hbox.addWidget(QtWidgets.QLabel(
                "Slices range:"))
            self.start_slice = QtWidgets.QLineEdit("")
            self.start_slice.setMaximumWidth(50)
            hbox.addWidget(self.start_slice)
            hbox.addWidget(QtWidgets.QLabel(
                "to"))
            self.end_slice = QtWidgets.QLineEdit("")
            self.end_slice.setMaximumWidth(50)
            hbox.addWidget(self.end_slice)
            hbox.addStretch(1)
            vbox.addLayout(hbox)

            # hbox = QtWidgets.QHBoxLayout()
            # hbox.addWidget(QtWidgets.QLabel(
            #     "Maximum number of optical slices:"))
            # self.nr_optical_slices = QtWidgets.QLineEdit("1")
            # self.nr_optical_slices.setMaximumWidth(50)
            # hbox.addWidget(self.nr_optical_slices)
            # hbox.addStretch(1)
            # vbox.addLayout(hbox)

            hbox = QtWidgets.QHBoxLayout()
            hbox.addWidget(QtWidgets.QLabel("Load channels:"))
            self.cb_C1 = QtWidgets.QCheckBox('1')
            self.cb_C1.setChecked(True)
            hbox.addWidget(self.cb_C1)
            self.cb_C2 = QtWidgets.QCheckBox('2')
            self.cb_C2.setChecked(True)
            hbox.addWidget(self.cb_C2)
            self.cb_C3 = QtWidgets.QCheckBox('3')
            self.cb_C3.setChecked(True)
            hbox.addWidget(self.cb_C3)
            self.cb_C4 = QtWidgets.QCheckBox('4')
            self.cb_C4.setChecked(True)
            hbox.addWidget(self.cb_C4)
            hbox.addStretch(1)
            vbox.addLayout(hbox)

            hbox = QtWidgets.QHBoxLayout()
            bLoad3D = QtWidgets.QPushButton('Load volume')
            bLoad3D.setCheckable(True)
            bLoad3D.clicked.connect(self.Load)
            hbox.addWidget(bLoad3D)
            # hbox.addStretch(1)
            # vbox.addLayout(hbox)

            # hbox = QtWidgets.QHBoxLayout()
            bLoad3D = QtWidgets.QPushButton('Reload in shape')
            bLoad3D.setCheckable(True)
            bLoad3D.clicked.connect(self.LoadInRegion)
            hbox.addWidget(bLoad3D)
            bCrop = QtWidgets.QPushButton('Crop to shape')
            bCrop.setCheckable(True)
            bCrop.clicked.connect(self.CropToRegion)
            hbox.addWidget(bCrop)
            hbox.addStretch(1)
            vbox.addLayout(hbox)

            hbox = QtWidgets.QHBoxLayout()
            hbox.addWidget(QtWidgets.QLabel("Tissue threshold value:"))
            self.thresholdN = QtWidgets.QLineEdit("0.3")
            hbox.addWidget(self.thresholdN)
            hbox.addStretch(1)
            vbox.addLayout(hbox)

            hbox = QtWidgets.QHBoxLayout()
            hbox.addWidget(QtWidgets.QLabel("Number of regions to retain:"))
            self.spinN = QtWidgets.QSpinBox()
            self.spinN.setValue(1)
            hbox.addWidget(self.spinN)
            hbox.addStretch(1)
            vbox.addLayout(hbox)

            hbox = QtWidgets.QHBoxLayout()
            hbox.addWidget(QtWidgets.QLabel("Minimal size:"))
            self.maxSizeN = QtWidgets.QLineEdit("5000")
            hbox.addWidget(self.maxSizeN)
            hbox.addStretch(1)
            vbox.addLayout(hbox)

            hbox = QtWidgets.QHBoxLayout()
            bKeepN = QtWidgets.QPushButton('Show only large regions')
            bKeepN.setCheckable(True)
            bKeepN.clicked.connect(self.Keep_n_Regions)
            hbox.addWidget(bKeepN)
            bRemoveN = QtWidgets.QPushButton('Remove small regions')
            bRemoveN.setCheckable(True)
            bRemoveN.clicked.connect(self.Remove_Small_Regions)
            hbox.addWidget(bRemoveN)
            hbox.addStretch(1)
            vbox.addLayout(hbox)

            bLoad2D = QtWidgets.QPushButton('Load slice')
            bLoad2D.setCheckable(True)
            bLoad2D.clicked.connect(self.Load2D)
            hbox = QtWidgets.QHBoxLayout()
            hbox.addWidget(bLoad2D)

            hbox.addWidget(QtWidgets.QLabel("Slice:"))
            self.scroll = QtWidgets.QScrollBar()
            self.scroll.setOrientation(QtCore.Qt.Horizontal)
            self.scroll.setRange(0, 100)
            self.scroll.setMinimumWidth(150)
            self.scroll.valueChanged.connect(self.set_image_slice_value)
            hbox.addWidget(self.scroll)

            self.image_slice = QtWidgets.QLineEdit("0")
            self.image_slice.setMaximumWidth(30)
            hbox.addWidget(self.image_slice)

            hbox.addStretch(1)
            vbox.addLayout(hbox)


            hbox = QtWidgets.QHBoxLayout()

            bSaveSlice = QtWidgets.QPushButton('Save slice')
            bSaveSlice.setCheckable(True)
            bSaveSlice.clicked.connect(self.SaveSlice)
            hbox.addWidget(bSaveSlice)

            bSaveVolume = QtWidgets.QPushButton('Save volume')
            bSaveVolume.setCheckable(True)
            bSaveVolume.clicked.connect(self.SaveVolume)
            hbox.addWidget(bSaveVolume)

            # bNormalize = QtWidgets.QPushButton('Normalize')
            # bNormalize.setCheckable(True)
            # bNormalize.clicked.connect(self.Normalize)
            # hbox.addWidget(bNormalize)

            # self.normalize_value = QtWidgets.QLineEdit("0")
            # self.normalize_value.setMaximumWidth(30)
            # hbox.addWidget(self.normalize_value)

            #bSaveMovie = QtWidgets.QPushButton('Save movie')
            #bSaveMovie.setCheckable(True)
            #bSaveMovie.clicked.connect(self.SaveMovie)
            #hbox.addWidget(bSaveMovie)
            hbox.addStretch(1)
            vbox.addLayout(hbox)

            #bStardist = QtWidgets.QPushButton('Run Stardist in shape')
            #bStardist.setCheckable(True)
            #bStardist.clicked.connect(self.run_stardist)
            #hbox = QtWidgets.QHBoxLayout()
            #hbox.addWidget(bStardist)
            #hbox.addStretch(1)
            #vbox.addLayout(hbox)

            #vbox.addStretch(1)
            
            bAddPolygon = QtWidgets.QPushButton('Add')
            bAddPolygon.setCheckable(True)
            bAddPolygon.clicked.connect(self.add_polygon)
            bRemoveOutside = QtWidgets.QPushButton('Remove outside interpolated region')
            bRemoveOutside.setCheckable(True)
            bRemoveOutside.clicked.connect(self.run_remove_outside)
            hbox = QtWidgets.QHBoxLayout()
            hbox.addWidget(bAddPolygon)
            hbox.addWidget(bRemoveOutside)
            hbox.addStretch(1)
            vbox.addLayout(hbox)

            vbox.addStretch(1)



            
            widget.setLayout(vbox)
            self.viewer.window.add_dock_widget(widget, area="right")

            animation_widget = AnimationWidget(self.viewer)
            self.viewer.window.add_dock_widget(animation_widget, area='right')

            #@viewer.bind_key('z')
            #def print_z_slice(viewer):
            #    print(viewer.dims.point[0])

            # inputfile = ''
            # try:
            #     opts, args = getopt.getopt(
            #         argv, "hi:c:", ["ifile=", "channel="])
            # except getopt.GetoptError:'
            #     sys.exit(2)
            # for opt, arg in opts:
            #     if opt == '-h':
            #         sys.exit()
            #     elif opt in ("-i", "--ifile"):
            #         inputfile = arg
            #     elif opt in ("-c", "--channel"):
            #         channel = channel


# if __name__ == "__main__":
#    NapariSTPT().main(sys.argv[1:])
