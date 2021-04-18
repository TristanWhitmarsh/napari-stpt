"""
napari-STPT reads zarr files and displays them
"""

import os
import sys
import napari
import numpy as np
import xarray as xr
import glob
from qtpy import QtCore, QtGui, QtWidgets
import SimpleITK as sitk
import getopt
#import cv2 as cv
#from skimage import filters

#from skimage.segmentation import watershed
#from skimage.feature import peak_local_max
from scipy import ndimage
#from skimage import measure
#from skimage.measure import regionprops
import math
import cv2


class NapariSTPT:

    def __init__(self):
        self.scroll = None
        self.image_slice = None

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
        self.thresholdN = None
        self.slice_names = None

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
        global ds1
        global ds
        global bscale
        global bzero
        global align_x
        global align_y
        global spacing

        print(self.search_folder.text() +
              str(self.comboBoxPath.currentText()) + '/mos.zarr')

        try:
            ds1 = xr.open_zarr(self.search_folder.text(
            ) + str(self.comboBoxPath.currentText())+'/mos.zarr', consolidated=True)
            print("consolidated")
        except:
            print("none-consolidated")
            ds1 = xr.open_zarr(self.search_folder.text() +
                               str(self.comboBoxPath.currentText())+'/mos.zarr')
        #ds2 = xr.open_zarr(str(self.comboBoxPath.currentText())+'/mos.zarr', group='l.2')
        #ds4 = xr.open_zarr(str(self.comboBoxPath.currentText())+'/mos.zarr', group='l.4')
        #ds8 = xr.open_zarr(str(self.comboBoxPath.currentText())+'/mos.zarr', group='l.8')
        #ds32 = xr.open_zarr(str(self.comboBoxPath.currentText())+'/mos.zarr', group='l.32')

        self.slice_names = ds1.attrs['cube_reg']['slice']

        bscale = ds1.attrs['bscale']
        bzero = ds1.attrs['bzero']

        # print(ds1.attrs.keys())
        # print(ds1['S001'].attrs.keys())
        # print(ds1['S001'].attrs)
        # print(ds1.attrs["multiscale"]['datasets'][1]['level'])
        # print(ds1.attrs["multiscale"]['datasets'][1]['path'])

        align_x = ds1.attrs['cube_reg']['abs_dx']
        align_y = ds1.attrs['cube_reg']['abs_dy']

        #print("size: {}".format(ds1.dims))

        resolution = int(self.comboBoxResolution.currentText())
        print("loading at resolution {}".format(resolution))
        index = self.comboBoxResolution.currentIndex()
        print("group: " + ds1.attrs["multiscale"]['datasets'][index]['path'])
        gr = ds1.attrs["multiscale"]['datasets'][index]['path']
        ds = xr.open_zarr(self.search_folder.text(
        ) + str(self.comboBoxPath.currentText())+'/mos.zarr', group=gr)

        output_resolution = float(self.pixel_size.text())

        with napari.gui_qt():
            if self.cb_C1.isChecked():
                print("loading C1")
                volume_1 = (ds.sel(channel=1, type='mosaic',
                            z=0).to_array().data * bscale + bzero)
                print("cropping C1")
                volume_1 = self.CropOriginal(volume_1)
                print("aligning C1")
                self.aligned_1 = self.Align(
                    volume_1, resolution, output_resolution)
                self.viewer.add_image([self.aligned_1], name='C1', scale=(
                    15, output_resolution, output_resolution), blending='additive', colormap='bop purple')

            if self.cb_C2.isChecked():
                print("loading C2")
                volume_2 = (ds.sel(channel=2, type='mosaic',
                            z=0).to_array().data * bscale + bzero)
                print("cropping C2")
                volume_2 = self.CropOriginal(volume_2)
                print("aligning C2")
                self.aligned_2 = self.Align(
                    volume_2, resolution, output_resolution)
                self.viewer.add_image([self.aligned_2], name='C2', scale=(
                    15, output_resolution, output_resolution), blending='additive', colormap='red')

            if self.cb_C3.isChecked():
                print("loading C3")
                volume_3 = (ds.sel(channel=3, type='mosaic',
                            z=0).to_array().data * bscale + bzero)
                print("cropping C3")
                volume_3 = self.CropOriginal(volume_3)
                print("aligning C3")
                self.aligned_3 = self.Align(
                    volume_3, resolution, output_resolution)
                self.viewer.add_image([self.aligned_3], name='C3', scale=(
                    15, output_resolution, output_resolution), blending='additive', colormap='green')

            if self.cb_C4.isChecked():
                print("loading C4")
                volume_4 = (ds.sel(channel=4, type='mosaic',
                            z=0).to_array().data * bscale + bzero)
                print("cropping C4")
                volume_4 = self.CropOriginal(volume_4)
                print("aligning C4")
                self.aligned_4 = self.Align(
                    volume_4, resolution, output_resolution)
                self.viewer.add_image([self.aligned_4], name='C4', scale=(
                    15, output_resolution, output_resolution), blending='additive', colormap='bop blue')

            if any(i.name == 'C1' for i in self.viewer.layers):
                print("is there")
            else:
                print("is not there")

    def Load2D(self, text):

        print(self.search_folder.text() +
              str(self.comboBoxPath.currentText()) + '/mos.zarr')

        ds1 = xr.open_zarr(self.search_folder.text() +
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

        bscale = ds1.attrs['bscale']
        bzero = ds1.attrs['bzero']

        z = self.scroll.value()

        self.slice_names = ds1.attrs['cube_reg']['slice']

        self.scroll.setRange(0, len(self.slice_names)-1)
        if z >= len(self.slice_names):
            z = len(self.slice_names)-1
            self.scroll.setValue(z)

        slice_name = self.slice_names[z]

        #z = 0
        #slice_name = 'S00'+str(z+1)
        # if(z+1 > 9):
        #    slice_name = 'S0'+str(z+1)
        # if(z+1 > 99):
        #    slice_name = 'S'+str(z+1)

        if self.cb_C1.isChecked():
            im1 = ds1[slice_name].sel(
                channel=1, type='mosaic', z=0).data * bscale + bzero
            im2 = ds2[slice_name].sel(
                channel=1, type='mosaic', z=0).data * bscale + bzero
            im4 = ds4[slice_name].sel(
                channel=1, type='mosaic', z=0).data * bscale + bzero
            im8 = ds8[slice_name].sel(
                channel=1, type='mosaic', z=0).data * bscale + bzero
            im16 = ds16[slice_name].sel(
                channel=1, type='mosaic', z=0).data * bscale + bzero
            im32 = ds32[slice_name].sel(
                channel=1, type='mosaic', z=0).data * bscale + bzero

            with napari.gui_qt():
                self.viewer.add_image([im1, im2, im4, im8, im16, im32], multiscale=True,
                                      name='C1', blending='additive', colormap='bop purple')

        if self.cb_C2.isChecked():
            im1 = ds1[slice_name].sel(
                channel=2, type='mosaic', z=0).data * bscale + bzero
            im2 = ds2[slice_name].sel(
                channel=2, type='mosaic', z=0).data * bscale + bzero
            im4 = ds4[slice_name].sel(
                channel=2, type='mosaic', z=0).data * bscale + bzero
            im8 = ds8[slice_name].sel(
                channel=2, type='mosaic', z=0).data * bscale + bzero
            im16 = ds16[slice_name].sel(
                channel=2, type='mosaic', z=0).data * bscale + bzero
            im32 = ds32[slice_name].sel(
                channel=2, type='mosaic', z=0).data * bscale + bzero

            with napari.gui_qt():
                self.viewer.add_image([im1, im2, im4, im8, im16, im32],
                                      multiscale=True, name='C2', blending='additive', colormap='red')

        if self.cb_C3.isChecked():
            im1 = ds1[slice_name].sel(
                channel=3, type='mosaic', z=0).data * bscale + bzero
            im2 = ds2[slice_name].sel(
                channel=3, type='mosaic', z=0).data * bscale + bzero
            im4 = ds4[slice_name].sel(
                channel=3, type='mosaic', z=0).data * bscale + bzero
            im8 = ds8[slice_name].sel(
                channel=3, type='mosaic', z=0).data * bscale + bzero
            im16 = ds16[slice_name].sel(
                channel=3, type='mosaic', z=0).data * bscale + bzero
            im32 = ds32[slice_name].sel(
                channel=3, type='mosaic', z=0).data * bscale + bzero

            with napari.gui_qt():
                self.viewer.add_image([im1, im2, im4, im8, im16, im32], multiscale=True,
                                      name='C3', blending='additive', colormap='green')
                #self.viewer.add_image(float_img, multiscale=False, name='C3', blending='additive', colormap='green')

        if self.cb_C4.isChecked():
            im1 = ds1[slice_name].sel(
                channel=4, type='mosaic', z=0).data * bscale + bzero
            im2 = ds2[slice_name].sel(
                channel=4, type='mosaic', z=0).data * bscale + bzero
            im4 = ds4[slice_name].sel(
                channel=4, type='mosaic', z=0).data * bscale + bzero
            im8 = ds8[slice_name].sel(
                channel=4, type='mosaic', z=0).data * bscale + bzero
            im16 = ds16[slice_name].sel(
                channel=4, type='mosaic', z=0).data * bscale + bzero
            im32 = ds32[slice_name].sel(
                channel=4, type='mosaic', z=0).data * bscale + bzero

            with napari.gui_qt():
                self.viewer.add_image([im1, im2, im4, im8, im16, im32], multiscale=True,
                                      name='C4', blending='additive', colormap='bop blue')

    def Align(self, volume, resolution, output_resolution):

        spacing = (ds1['S001'].attrs['scale'])
        #size_multiplier = (resolution*0.1*spacing[0])/7.5
        # print(resolution*0.1*spacing[0])
        size_multiplier = (resolution*0.1*spacing[0])/output_resolution
        size = (volume.shape[0], int(
            size_multiplier*volume.shape[1]), int(size_multiplier*volume.shape[2]))
        aligned = np.zeros(size, dtype=float)
        size2D = (
            int(size_multiplier*volume.shape[2]), int(size_multiplier*volume.shape[1]))

        global align_x
        global align_y

        z_size = volume.shape[0]
        for z in range(0, z_size):

            fixed = sitk.GetImageFromArray(volume[z, :, :])
            fixed.SetOrigin((0, 0))

            #slice_name = 'S00'+str(z+1)
            # if(z+1 > 9):
            #    slice_name = 'S0'+str(z+1)
            # if(z+1 > 99):
            #    slice_name = 'S'+str(z+1)
            slice_name = self.slice_names[z]

            spacing = (ds1[slice_name].attrs['scale'])
            # print(spacing)
            # fixed.SetSpacing([1.6*spacing[0],1.6*spacing[1]])
            fixed.SetSpacing([resolution*0.1*spacing[1],
                             resolution*0.1*spacing[0]])

            transform = sitk.Euler2DTransform()

            alignY = 0
            if not np.isnan(align_y[z]):
                alignY = -align_y[z]*0.1*spacing[1]

            alignX = 0
            if not np.isnan(align_x[z]):
                alignX = -align_x[z]*0.1*spacing[0]

            transform.SetTranslation([alignY, alignX])

            # transform.SetTranslation([-align_y[z]*0.1*spacing[1],-align_x[z]*0.1*spacing[0]])
            # transform.SetTranslation([0.5*-align_y[z],0.5*-align_x[z]])
            # transform.SetTranslation([0,0])
            print('{}/{}, translate_x: {}, translate_y: {}'.format(z,
                  z_size, alignX, alignY))

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

        return aligned

    def Keep_n_Regions(self):

        output_resolution = float(self.pixel_size.text())
        threshold = float(self.thresholdN.text())

        if self.aligned_1 is not None and self.cb_C1.isChecked():

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

                max_size = sizes_sorted[len(sizes_sorted)-1-keep_n]
                for i in range(1, nb_components):
                    if sizes[i] < max_size:
                        volume_z[output == i] = threshold

                self.aligned_1[z, :, :] = volume_z

            self.viewer.layers['C1'].visible = False
            self.viewer.layers['C1'].visible = True

        if self.aligned_2 is not None and self.cb_C2.isChecked():

            threholded = self.aligned_2 > threshold
            threholded = ndimage.binary_fill_holes(threholded)
            threholded = threholded.astype(np.uint8)

            #regions = threholded.copy()

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

                #regions[z,:,:] = output

                max_size = sizes_sorted[len(sizes_sorted)-1-keep_n]
                for i in range(1, nb_components):
                    if sizes[i] < max_size:
                        volume_z[output == i] = threshold

                self.aligned_2[z, :, :] = volume_z

            self.viewer.layers['C2'].visible = False
            self.viewer.layers['C2'].visible = True

            # with napari.gui_qt():
            #    self.viewer.add_image(regions, name='C2 regions', scale=(15,output_resolution,output_resolution), blending='additive', colormap='blue')

        if self.aligned_3 is not None and self.cb_C3.isChecked():

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

                max_size = sizes_sorted[len(sizes_sorted)-1-keep_n]
                for i in range(1, nb_components):
                    if sizes[i] < max_size:
                        volume_z[output == i] = threshold

                self.aligned_3[z, :, :] = volume_z

            self.viewer.layers['C3'].visible = False
            self.viewer.layers['C3'].visible = True

        if self.aligned_4 is not None and self.cb_C4.isChecked():

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

                max_size = sizes_sorted[len(sizes_sorted)-1-keep_n]
                for i in range(1, nb_components):
                    if sizes[i] < max_size:
                        volume_z[output == i] = threshold

                self.aligned_4[z, :, :] = volume_z

            self.viewer.layers['C4'].visible = False
            self.viewer.layers['C4'].visible = True

        #self.viewer.layers['C2'].contrast_limits=(0, 3)
        # with napari.gui_qt():
        #    if self.cb_C2.isChecked():
        #        self.viewer.add_image([self.aligned_2], name='C2 processed', scale=(15,output_resolution,output_resolution), blending='additive', colormap='red')

    def Crop(self):

        output_resolution = float(self.pixel_size.text())
        threshold = float(self.thresholdN.text())

        maxX_final = 0
        minX_final = sys.maxsize
        maxY_final = 0
        minY_final = sys.maxsize
        maxZ_final = 0
        minZ_final = sys.maxsize

        if self.aligned_1 is not None and self.cb_C1.isChecked():

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

        if self.aligned_2 is not None and self.cb_C2.isChecked():

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

        if self.aligned_3 is not None and self.cb_C3.isChecked():

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

        if self.aligned_4 is not None and self.cb_C4.isChecked():

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

        if self.aligned_1 is not None and self.cb_C1.isChecked():
            self.aligned_1 = self.aligned_1[minX_final:maxX_final,
                                            minY_final:maxY_final, minZ_final:maxZ_final]
            self.viewer.layers['C1'].visible = False
            self.viewer.layers['C1'].visible = True
        if self.aligned_2 is not None and self.cb_C2.isChecked():
            self.aligned_2 = self.aligned_2[minX_final:maxX_final,
                                            minY_final:maxY_final, minZ_final:maxZ_final]
            self.viewer.layers.remove('C2')
            self.viewer.add_image([self.aligned_2], name='C2', scale=(
                15, output_resolution, output_resolution), blending='additive', colormap='red')
        if self.aligned_3 is not None and self.cb_C3.isChecked():
            self.aligned_3 = self.aligned_3[minX_final:maxX_final,
                                            minY_final:maxY_final, minZ_final:maxZ_final]
            self.viewer.layers['C3'].visible = False
            self.viewer.layers['C3'].visible = True
        if self.aligned_4 is not None and self.cb_C4.isChecked():
            self.aligned_4 = self.aligned_4[minX_final:maxX_final,
                                            minY_final:maxY_final, minZ_final:maxZ_final]
            self.viewer.layers['C4'].visible = False
            self.viewer.layers['C4'].visible = True

    def SelectFolder(self):
        file = str(QtWidgets.QFileDialog.getExistingDirectory())
        file = file + '/'
        #text_files = glob.glob(path + "/**/**/*.zarr", recursive = False)
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
            self.comboBoxResolution.clear()
            # for n in reversed(range(0, length)):
            for n in range(0, length):
                print("adding resolution {}".format(
                    ds.attrs["multiscale"]['datasets'][n]['level']))
                self.comboBoxResolution.addItem(
                    str(ds.attrs["multiscale"]['datasets'][n]['level']))
                self.comboBoxResolution.setCurrentText(
                    str(ds.attrs["multiscale"]['datasets'][n]['level']))
            try:
                self.scroll.setValue(0)
                self.image_slice.setText("0")
                self.slice_names = ds.attrs['cube_reg']['slice']
                self.scroll.setRange(0, len(self.slice_names))
            except:
                #print("not initialised")
                pass
        else:
            print("adding default resolution")
            self.comboBoxResolution.clear()
            self.comboBoxResolution.addItem(str(1))
            self.comboBoxResolution.addItem(str(2))
            self.comboBoxResolution.addItem(str(4))
            self.comboBoxResolution.addItem(str(8))
            self.comboBoxResolution.addItem(str(16))
            self.comboBoxResolution.addItem(str(32))
            self.comboBoxResolution.setCurrentText(str(32))

    def set_image_slice_value(self):
        value = self.scroll.value()
        self.image_slice.setText(str(value))

    def main(self, argv):
        with napari.gui_qt():

            loaded16_2 = 0
            loaded16_3 = 0
            loaded16_4 = 0

            self.viewer = napari.Viewer()

            widget = QtWidgets.QWidget()
            vbox = QtWidgets.QVBoxLayout()

            hbox = QtWidgets.QHBoxLayout()
            hbox.addWidget(QtWidgets.QLabel("Data folder:"))
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

            self.comboBoxPath.clear()
            if os.path.exists(self.search_folder.text()):
                for f in os.scandir(self.search_folder.text()):
                    if f.is_dir():
                        if os.path.exists(f.path + "/mos.zarr"):
                            s = f.path
                            s = s.replace(self.search_folder.text(), '')
                            print(s)
                            self.comboBoxPath.addItem(s)

            self.comboBoxPath.setMaximumWidth(300)
            self.comboBoxPath.currentIndexChanged.connect(
                self.on_combobox_changed)
            hbox.addWidget(self.comboBoxPath)
            vbox.addLayout(hbox)

            hbox = QtWidgets.QHBoxLayout()
            hbox.addWidget(QtWidgets.QLabel("Load level:"))
            self.comboBoxResolution = QtWidgets.QComboBox()
            self.comboBoxResolution.setMaximumWidth(100)
            self.on_combobox_changed()
            hbox.addWidget(self.comboBoxResolution)
            hbox.addStretch(1)

            hbox.addWidget(QtWidgets.QLabel("Output pixel size:"))
            self.pixel_size = QtWidgets.QLineEdit("7.5")
            self.pixel_size.setMaximumWidth(50)
            hbox.addWidget(self.pixel_size)
            hbox.addStretch(1)
            vbox.addLayout(hbox)

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
            bLoad3D = QtWidgets.QPushButton('Load image 3D')
            bLoad3D.setCheckable(True)
            bLoad3D.clicked.connect(self.Load)
            hbox.addWidget(bLoad3D)
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
            hbox.addWidget(QtWidgets.QLabel("Tissue threshold value:"))
            self.thresholdN = QtWidgets.QLineEdit("0.1")
            hbox.addWidget(self.thresholdN)
            hbox.addStretch(1)
            vbox.addLayout(hbox)

            hbox = QtWidgets.QHBoxLayout()
            bKeepN = QtWidgets.QPushButton('Show only large regions')
            bKeepN.setCheckable(True)
            bKeepN.clicked.connect(self.Keep_n_Regions)
            hbox.addWidget(bKeepN)
            bCrop = QtWidgets.QPushButton('Crop')
            bCrop.setCheckable(True)
            bCrop.clicked.connect(self.Crop)
            hbox.addWidget(bCrop)
            hbox.addStretch(1)
            vbox.addLayout(hbox)

            bLoad2D = QtWidgets.QPushButton('Load image 2D')
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

            vbox.addStretch(1)
            widget.setLayout(vbox)
            self.viewer.window.add_dock_widget(widget, area="right")

            inputfile = ''
            try:
                opts, args = getopt.getopt(
                    argv, "hi:c:", ["ifile=", "channel="])
            except getopt.GetoptError:
                # print 'napari-stpt -i <inputfile> -c <channel>'
                sys.exit(2)
            for opt, arg in opts:
                if opt == '-h':
                    # print 'napari-stpt'
                    sys.exit()
                elif opt in ("-i", "--ifile"):
                    inputfile = arg
                elif opt in ("-c", "--channel"):
                    channel = channel
            # print 'Input file is "', inputfile
            # print 'channel is "', channel


if __name__ == "__main__":
    NapariSTPT().main(sys.argv[1:])
