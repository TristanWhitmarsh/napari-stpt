"""
napari-STPT reads zarr files and displays them
"""

import os
import sys
import re
import napari
from napari._qt.widgets import qt_dims_slider

try:
    from .custom_qt_dims_slider import CustomQtDimSliderWidget
except ImportError:
    from custom_qt_dims_slider import CustomQtDimSliderWidget

# Apply monkey patch
napari._qt.widgets.qt_dims_slider.QtDimSliderWidget._update_range = CustomQtDimSliderWidget._update_range
napari._qt.widgets.qt_dims_slider.QtDimSliderWidget._update_slice_labels = CustomQtDimSliderWidget._update_slice_labels

import numpy as np
import xarray as xr
import vispy.color
import colorsys

# if 'PyQt5' in sys.modules:
#     print("Using PyQt5")
#     from qtpy import QtCore, QtWidgets
#     from qtpy.QtCore import Qt, QSortFilterProxyModel
#     from qtpy.QtGui import QColor, QPixmap, QIcon, QStandardItemModel, QStandardItem, QPainter, QFont
#     from qtpy.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QFrame, QScrollArea,
#                                 QComboBox, QCheckBox, QListWidget, QListWidgetItem, QMessageBox, QFileDialog, QMainWindow)
# else:
#print("Using PySide2")
from PySide2 import QtCore, QtWidgets, QtGui
from PySide2.QtCore import Qt, QSortFilterProxyModel
from PySide2.QtGui import QColor, QPixmap, QIcon, QStandardItemModel, QStandardItem, QPainter, QFont
from PySide2.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QFrame, QScrollArea,
                               QComboBox, QCheckBox, QListWidget, QListWidgetItem, QMessageBox, QFileDialog, QMainWindow)
 
import SimpleITK as sitk
from scipy import ndimage, stats
import cv2
#from stardist.models import StarDist2D
from napari_animation import AnimationWidget
from PIL import Image, ImageDraw
import random
# import skimage.io
# from stardist.models import StarDist3D
#from csbdeep.utils import normalize
#from naparimovie import Movie
import math
import gc
import json
import tifffile
import pandas as pd
from tifffile import imwrite
from tifffile import TiffWriter
from skimage.transform import resize

try:
    from .my_widgets import ExtendedComboBox
    from .my_widgets import CheckableItemWidget
    from .my_widgets import CustomComboBox
    from .my_widgets import LegendWidget
    from .utilities import parse_channel_input
    from .data import Data 
except ImportError:
    from my_widgets import ExtendedComboBox
    from my_widgets import CheckableItemWidget
    from my_widgets import CustomComboBox
    from my_widgets import LegendWidget
    from utilities import parse_channel_input
    from data import Data
    

        
        
       

class NapariSTPT:

    def __init__(self):
        self.legend_widget = None
        self.scroll = None
        self.scroll_overall_brightness = None
        self.overall_brightness = 1
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
        self.current_output_resolution = None
        self.output_resolution = None

        self.cb_R_Axio = None
        # self.cb_R_Old = None
        self.cb_R_STPT = None
        self.spacing = [0,0,0]

        self.cb_perspective = None
        self.cb_isometric = None
        
        self.origin_x = None
        self.origin_y = None
        self.crop_start_x = None
        self.crop_start_y = None
        self.crop_end_x = None
        self.crop_end_y = None

        self.spinN = None
        self.maxSizeN = None
        self.thresholdN = None
        self.slice_names = None

        self.align_x = None
        self.align_y = None
        # data.ds1 = None
        # data.ds2 = None
        # data.ds4 = None
        # data.ds8 = None
        # data.ds16 = None
        # data.ds32 = None
        
        self.m_volume_1 = None
        self.m_volume_1_multiplier = None
        self.m_volume_2 = None
        self.m_volume_2_multiplier = None
        self.m_volume_new = None

        #self.movie = None
        self.optical_slices = None
        self.nr_optical_slices = None
        self.optical_slices_available = None
        self.cb_correct_brightness_optical_section = None
        self.shape = None
        self.spacing_loaded = None

        self.crop_start_ratio_x = None
        self.crop_size_ratio_x = None
        self.crop_start_ratio_y = None
        self.crop_size_ratio_y = None

        self.crop = False
        self.old_method = False

        self.start_slice = None
        self.end_slice = None

        self.image_translation = None
        self.m_slice_spacing = None
        self.slice_spacing = None
        self.shapeText = None
        
        self.default_contrast_limits = None
        self.channels_start_at_0 = None
        self.channels_start = None
        
        self.bscale = None
        self.bzero = None
        
        self.layerC1 = None
        self.layerC2 = None
        self.layerC3 = None
        self.layerC4 = None
        
        self.loaded_2D = None
        self.loaded_3D = None
        
        self.bLoad3D = None
        self.bLoad2D = None
        
        self.number_of_sections = None
        
        self.image_folder = None
        self.axio = False
        self.selected_slices = None
        self.selected_channels = None
        self.value_range = None
        
        self.data = Data()
        
        self.ignore_gui_call = False;

    def Make3DShape(self):

        output_resolution = float(self.pixel_size.text())
        
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

        #print("crop to: {} {} {} {}".format(minX, maxX, minY, maxY))

        
        size_x = self.shape[0] * float(self.slice_spacing)/float(self.optical_slices) / output_resolution
        size_y = maxX
        size_z = maxY

        line_locations = []
        line_locations.append([[0, minX, minY], [size_x, minX, minY]])
        line_locations.append([[size_x, minX, minY], [size_x, size_y, minY]])
        line_locations.append([[size_x, size_y, minY], [0, size_y, minY]])
        line_locations.append([[0, size_y, minY], [0, minX, minY]])

        line_locations.append([[0, minX, size_z], [size_x, minX, size_z]])
        line_locations.append([[size_x, minX, size_z], [size_x, size_y, size_z]])
        line_locations.append([[size_x, size_y, size_z], [0, size_y, size_z]])
        line_locations.append([[0, size_y, size_z], [0, minX, size_z]])
        
        line_locations.append([[0, minX, minY], [0, minX, size_z]])
        line_locations.append([[size_x, minX, minY], [size_x, minX, size_z]])
        line_locations.append([[size_x, size_y, minY], [size_x, size_y, size_z]])
        line_locations.append([[0, size_y, minY], [0, size_y, size_z]])

        width = np.mean([size_x, size_y, size_z]) / 500
        
        self.viewer.add_shapes(np.asarray(line_locations), name = 'crop box',shape_type='line', edge_color = "white", edge_width = width)



    def MakeBoundingBox(self):
        #with napari.gui_qt():
        
        size_x = (self.shape[0]-1)# * float(self.slice_spacing)/float(self.optical_slices) / output_resolution) - 1
        size_y = (self.shape[1]-1)# * self.output_resolution
        size_z = (self.shape[2]-1)# * self.output_resolution

        line_locations = []
        line_locations.append([[0, 0, 0], [size_x, 0, 0]])
        line_locations.append([[size_x, 0, 0], [size_x, size_y, 0]])
        line_locations.append([[size_x, size_y, 0], [0, size_y, 0]])
        line_locations.append([[0, size_y, 0], [0, 0, 0]])

        line_locations.append([[0, 0, size_z], [size_x, 0, size_z]])
        line_locations.append([[size_x, 0, size_z], [size_x, size_y, size_z]])
        line_locations.append([[size_x, size_y, size_z], [0, size_y, size_z]])
        line_locations.append([[0, size_y, size_z], [0, 0, size_z]])
        
        line_locations.append([[0, 0, 0], [0, 0, size_z]])
        line_locations.append([[size_x, 0, 0], [size_x, 0, size_z]])
        line_locations.append([[size_x, size_y, 0], [size_x, size_y, size_z]])
        line_locations.append([[0, size_y, 0], [0, size_y, size_z]])

        width = np.mean([size_x, size_y, size_z]) / 500
        
        #output_resolution = float(self.pixel_size.text())
        scale_x = float(self.slice_spacing)/float(self.optical_slices)# / output_resolution
        self.viewer.add_shapes(np.asarray(line_locations), name = 'bounding box',shape_type='line', scale=(float(self.slice_spacing)/float(self.optical_slices), self.output_resolution, self.output_resolution), edge_color = "white", edge_width = width)

        
    def parse_channel_input(self, input_string):
        channels = []
        # Remove spaces and split input by commas
        input_string = input_string.replace(" ", "")
        parts = input_string.split(',')

        for part in parts:
            if '-' in part:
                # Handle ranges
                start, end = part.split('-')
                try:
                    start = int(start)
                    end = int(end)
                    channels.extend(range(start, end + 1))
                except ValueError:
                    #print(f"Invalid range: {part}")
                    pass
            else:
                # Single channel
                try:
                    channels.append(int(part))
                except ValueError:
                    #print(f"Invalid channel: {part}")
                    pass

        return channels
        
        
        
        

    def LoadInRegion(self, text):
        
        zoom = self.viewer.window.qt_viewer.camera.zoom
        angles = self.viewer.window.qt_viewer.camera.angles
        
        if any(i.name == 'bounding box' for i in self.viewer.layers):
            self.viewer.layers.remove('bounding box')
        
        if any(i.name == 'crop box' for i in self.viewer.layers):
            self.viewer.layers.remove('crop box')

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


        print("crop to: {} {} {} {}".format(minX, maxX, minY, maxY))
        print("self.origin_x: {}".format(self.origin_x))
        print("self.origin_y: {}".format(self.origin_y))
        print("self.current_output_resolution: {}".format(self.current_output_resolution))
        
        self.crop_start_x = self.origin_x + (minX * self.current_output_resolution)
        self.crop_start_y = self.origin_y + (minY * self.current_output_resolution)
        self.crop_end_x = self.origin_x + (maxX * self.current_output_resolution)
        self.crop_end_y = self.origin_y + (maxY * self.current_output_resolution)

        self.origin_x = self.crop_start_x
        self.origin_y = self.crop_start_y
        
        
        for layer in self.viewer.layers:
            if isinstance(layer, napari.layers.Image):
                data = layer.data  # Get the data of the layer
                self.crop_start_ratio_x = self.aligned_1.shape[1]/minX
                self.crop_size_ratio_x = self.aligned_1.shape[1]/(maxX-minX)
                self.crop_start_ratio_y = self.aligned_1.shape[2]/minY
                self.crop_size_ratio_y = self.aligned_1.shape[2]/(maxY-minY)
                break
                
        self.viewer.layers.remove('Shapes')
        
        # for layer in self.viewer.layers[:]:
        #     self.viewer.layers.remove(layer)

        print(f"crop_start_x {self.crop_start_x}, crop_start_y {self.crop_start_y}, crop_size_x {self.crop_end_x}, crop_size_y {self.crop_end_y}")

        self.crop = True
        self.Load(text)
        self.crop = False

        # self.MakeBoundingBox()
        
        self.viewer.window.qt_viewer.camera.center = (self.shape[0]/2, self.shape[1]/2, self.shape[2]/2 )
        self.viewer.window.qt_viewer.camera.angles = angles
        self.viewer.window.qt_viewer.camera.zoom = zoom


    def CropToRegion(self):

        if any(i.name == 'bounding box' for i in self.viewer.layers):
            self.viewer.layers.remove('bounding box')
        
        if any(i.name == 'crop box' for i in self.viewer.layers):
            self.viewer.layers.remove('crop box')

        output_resolution = float(self.pixel_size.text())
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



        self.origin_x = self.origin_x + (output_resolution * minX)
        self.origin_y = self.origin_y + (output_resolution * minY)

        print("crop to: {} {} {} {}".format(minX, maxX, minY, maxY))
        
        
        
        layers_copy = self.viewer.layers.copy()
        
        for layer in layers_copy:
            if isinstance(layer, napari.layers.Image):
                desired_layer = layer
                layer_name = layer.name
                print("layer_name: {}".format(layer_name))

                data = layer.data  # Get the data of the layer
                layer_colormap = layer.colormap.name  # Get the name of the colormap
                print("layer_colormap: {}".format(layer_colormap))
                layer_blending = layer.blending  # Get the blending method
                print("layer_blending: {}".format(layer_blending))
                layer_contrast_limits = layer.contrast_limits
                print("layer_contrast_limits: {}".format(layer_contrast_limits))
                layer_scale = layer.scale
                print("layer_scale: {}".format(layer_scale))

                data = data[:, int(minX):int(maxX), int(minY):int(maxY)]
                
                self.shape = data.shape

                self.viewer.layers.remove(layer_name)

                self.image_translation = (int(minX*output_resolution), int(minY*output_resolution))

                self.viewer.add_image([data], name=layer_name, scale=layer_scale, 
                                  blending=layer_blending, colormap=layer_colormap, contrast_limits=layer_contrast_limits)

        
            
        self.viewer.layers.remove('Shapes')
        self.MakeBoundingBox()
        
        self.viewer.window.qt_viewer.camera.center = (self.shape[0]/2, self.shape[1]/2, self.shape[2]/2 )


          
    def set_image_slice_value(self):
        
        if not self.loaded_2D:
            return
        if self.ignore_gui_call:
            return
        
        self.image_slice.setText(str(self.scroll.value()))  
        #print(f"self.optical_slices_available {self.optical_slices_available}")
        optical_slice = (self.scroll.value()-1) % self.optical_slices_available
        #print(f"optical_slice {optical_slice}")
        z = math.floor((self.scroll.value()-1) / self.optical_slices_available)
        #print(f"z {z}")
        
        channel_names = self.data.ds1.coords['channel'].values.tolist()
        #print(f"channel_names: {channel_names}")
        
        
        if self.old_method:
            #Old method
            self.bscale = self.data.ds1.attrs['bscale']
            self.bzero = self.data.ds1.attrs['bzero']

            self.slice_names = data.ds1.attrs['cube_reg']['slice']

            # self.scroll.setRange(1, len(self.slice_names))
            # if z >= len(self.slice_names):
            #     z = len(self.slice_names)-1
            #     self.scroll.setValue(z)

            slice_name = self.slice_names[z]
        else:
            try:
                self.bscale = self.data.ds1['S001'].attrs['bscale']
                self.bzero = self.data.ds1['S001'].attrs['bzero']
            except:
                self.bscale = 1
                self.bzero = 0
            slice_name = f"S{(z+1):03d}"
        
        #print(self.selected_channels)
        for chn in range(50):
            if chn in self.selected_channels:
                #print("loading")

                        
                try:
                    im1 = (self.data.ds1[slice_name].sel(type='mosaic', z=optical_slice).data[chn] * self.bscale + self.bzero).squeeze()
                    im2 = (self.data.ds2[slice_name].sel(type='mosaic', z=optical_slice).data[chn] * self.bscale + self.bzero).squeeze()
                    im4 = (self.data.ds4[slice_name].sel(type='mosaic', z=optical_slice).data[chn] * self.bscale + self.bzero).squeeze()
                    im8 = (self.data.ds8[slice_name].sel(type='mosaic', z=optical_slice).data[chn] * self.bscale + self.bzero).squeeze()
                    im16 = (self.data.ds16[slice_name].sel(type='mosaic', z=optical_slice).data[chn] * self.bscale + self.bzero).squeeze()
                    im32 = (self.data.ds32[slice_name].sel(type='mosaic', z=optical_slice).data[chn] * self.bscale + self.bzero).squeeze()
                except:
                    try:
                        im1 = (self.data.ds1[slice_name].sel(z=optical_slice).data[chn] * self.bscale + self.bzero).squeeze()
                        im2 = (self.data.ds2[slice_name].sel(z=optical_slice).data[chn] * self.bscale + self.bzero).squeeze()
                        im4 = (self.data.ds4[slice_name].sel(z=optical_slice).data[chn] * self.bscale + self.bzero).squeeze()
                        im8 = (self.data.ds8[slice_name].sel(z=optical_slice).data[chn] * self.bscale + self.bzero).squeeze()
                        im16 = (self.data.ds16[slice_name].sel(z=optical_slice).data[chn] * self.bscale + self.bzero).squeeze()
                        im32 = (self.data.ds32[slice_name].sel(z=optical_slice).data[chn] * self.bscale + self.bzero).squeeze()
                    except:
                        try:
                            im1 = (self.data.ds1[slice_name].data[chn] * self.bscale + self.bzero).squeeze()
                            im2 = (self.data.ds2[slice_name].data[chn] * self.bscale + self.bzero).squeeze()
                            im4 = (self.data.ds4[slice_name].data[chn] * self.bscale + self.bzero).squeeze()
                            im8 = (self.data.ds8[slice_name].data[chn] * self.bscale + self.bzero).squeeze()
                            im16 = (self.data.ds16[slice_name].data[chn] * self.bscale + self.bzero).squeeze()
                            im32 = (self.data.ds32[slice_name].data[chn] * self.bscale + self.bzero).squeeze()
                        except:
                            #print("skipping this channel since it can't be read")
                            continue
                        

            
                channel_name = channel_names[chn]
                #self.viewer.layers[channel_name].data = [im1, im2, im4, im8, im16, im32]
                
                # Check if the layer exists before attempting to update
                if str(channel_name) in [str(layer.name) for layer in self.viewer.layers]:
                    self.viewer.layers[str(channel_name)].data = [im1, im2, im4, im8, im16, im32]


        
        


    def Remove_Regions(self, use_size):

        # output_resolution = float(self.pixel_size.text())
        threshold = float(self.thresholdN.text())
        
        for i in self.viewer.layers:
            name = i.name
            
            if isinstance(i.data, np.ndarray):
                # Check if the layer contains numpy data
                volume_data = i.data

                threholded = volume_data > threshold
                threholded = ndimage.binary_fill_holes(threholded)
                threholded = threholded.astype(np.uint8)

                keep_n = self.spinN.value()

                for z in range(0, threholded.shape[0]):
                    print('{}/{}'.format(z+1, threholded.shape[0]))

                    threholded_z = threholded[z, :, :]
                    volume_z = volume_data[z, :, :]

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

                    volume_data[z, :, :] = volume_z

                self.viewer.layers[name].visible = False
                self.viewer.layers[name].visible = True

    def Remove_Small_Regions(self):
        self.Remove_Regions(True)

    def Keep_n_Regions(self):
        self.Remove_Regions(False)
    
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

        output_resolution = float(self.pixel_size.text())
        
        pos = self.viewer.dims.point[0] / ((float(self.slice_spacing)/float(self.optical_slices)) / output_resolution)
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

        new_shape = []

        if pos < contour_list_sorted[0][0]:
            for i in range(0, data_length):
                x = pos
                y = float(contour_list_sorted[0][1][i][1])
                z = float(contour_list_sorted[0][1][i][2])

                print(f"x {x}, y {y}, z {z}")

                new_shape.append((x,y,z))

        elif pos > contour_list_sorted[len(contour_list_sorted)-1][0]:
            for i in range(0, data_length):
                x = pos
                y = float(contour_list_sorted[len(contour_list_sorted)-1][1][i][1])
                z = float(contour_list_sorted[len(contour_list_sorted)-1][1][i][2])

                print(f"x {x}, y {y}, z {z}")

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
        shapes_layer = self.viewer.add_shapes(new_shape, shape_type='polygon', name = "Shapes", scale=(float(self.slice_spacing)/float(self.optical_slices) / output_resolution, 1, 1),)

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
                self.viewer.layers[layer_name].visible = False
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
            
            if self.aligned_1 is not None:
                aligned1_tmp[z_level,:,:] = aligned1_tmp[z_level,:,:] * mask
            if self.aligned_2 is not None:
                aligned2_tmp[z_level,:,:] = aligned2_tmp[z_level,:,:] * mask
            if self.aligned_3 is not None:
                aligned3_tmp[z_level,:,:] = aligned3_tmp[z_level,:,:] * mask
            if self.aligned_4 is not None:
                aligned4_tmp[z_level,:,:] = aligned4_tmp[z_level,:,:] * mask

        
           

        if self.aligned_1 is not None:
            self.viewer.add_image([aligned1_tmp], name='C1_masked', scale=(float(self.slice_spacing)/float(self.optical_slices) / output_resolution, 1, 1), 
                                  blending='additive', colormap='bop purple', contrast_limits=self.default_contrast_limits)
            self.viewer.layers['C1'].visible = False
        if self.aligned_2 is not None:
            self.viewer.add_image([aligned2_tmp], name='C2_masked', scale=(float(self.slice_spacing)/float(self.optical_slices) / output_resolution, 1, 1), 
                                  blending='additive', colormap='red', contrast_limits=self.default_contrast_limits)
            self.viewer.layers['C2'].visible = False
        if self.aligned_3 is not None:
            self.viewer.add_image([aligned3_tmp], name='C3_masked', scale=(float(self.slice_spacing)/float(self.optical_slices) / output_resolution, 1, 1), 
                                  blending='additive', colormap='green', contrast_limits=self.default_contrast_limits)
            self.viewer.layers['C3'].visible = False
        if self.aligned_4 is not None:
            self.viewer.add_image([aligned4_tmp], name='C4_masked', scale=(float(self.slice_spacing)/float(self.optical_slices) / output_resolution, 1, 1), 
                                  blending='additive', colormap='bop blue', contrast_limits=self.default_contrast_limits)
            self.viewer.layers['C4'].visible = False

        self.viewer.layers['Shapes'].visible = False





    def SaveVolume(self):
        widget = QtWidgets.QWidget()
        if sys.platform == 'linux':
            fname, _ = QtWidgets.QFileDialog.getSaveFileName(
                widget, 'Save file as', '/home/', "Image files (*.tiff *.mha *.nii)")
        else:
            fname, _ = QtWidgets.QFileDialog.getSaveFileName(
                widget, 'Save file as', 'c:\\', "Image files (*.tiff *.mha *.nii)")

        if fname != "":
            print(fname)

            file_name = fname#"/home/tristan/test.tiff"

            output_resolution = float(self.pixel_size.text())
            output_resolution_z = float(self.slice_spacing)/float(self.optical_slices)
            spacing = (output_resolution, output_resolution, output_resolution_z)
            print(f"spacing {spacing}")
            
            volume = []
            
            
            for layer in self.viewer.layers:
                if isinstance(layer, napari.layers.Image):
                    volume.append(layer.data)
            

            name, extension = os.path.splitext(file_name)
            print(f"extension {extension}")
            
            if extension == ".mha" or extension == ".nii" :
                
                volume = np.array(volume)
                print(f"volume {volume.shape}")
                volume = np.moveaxis(volume,0,3)
                print(f"spacing {spacing}")
                print(f"volume {volume.shape}")

                volume_itk = sitk.GetImageFromArray(volume)
                caster = sitk.CastImageFilter()
                caster.SetOutputPixelType(sitk.sitkVectorFloat32)
                volume_itk = caster.Execute(volume_itk)
                volume_itk.SetSpacing(spacing)
                
                writer = sitk.ImageFileWriter()
                writer.SetFileName(file_name)
                writer.UseCompressionOff()
                #writer.SetCompressionLevel(0)
                writer.Execute(volume_itk)
        
            else:

                volume = np.array(volume)
                print(f"volume {volume.shape}")
                volume = np.moveaxis(volume,0,1)
                print(f"volume {volume.shape}")
                print(volume.shape)
                
                # tifffile.imsave(file_name, volume.astype('float32'))
                # tifffile.imwrite(file_name, volume.astype('float32'), compress=9, imagej=True, metadata={'spacing': output_resolution_z, 'unit': 'um', 'axes': 'ZCYX'})
                # tifffile.imwrite(file_name, volume.astype('float32'), compression="zlib", imagej=True, metadata={'spacing': output_resolution_z, 'unit': 'um', 'axes': 'ZCYX'})
                tifffile.imwrite(file_name, volume.astype('float32'), compression="zlib", compressionargs={'level':5}, imagej=True, resolution=(1/output_resolution, 1/output_resolution), metadata={'spacing': output_resolution_z, 'unit': 'um', 'axes': 'ZCYX'})

            
    def LoadVolume(self):
        widget = QtWidgets.QWidget()
        if sys.platform == 'linux':
            fname, _ = QtWidgets.QFileDialog.getOpenFileName(
                widget, 'Load file', '/home/', "Image files (*.tiff *.mha *.nii)")
        else:
            fname, _ = QtWidgets.QFileDialog.getOpenFileName(
                widget, 'Load file', 'c:\\', "Image files (*.tiff *.mha *.nii)")

        if fname != "":
            print(fname)

            file_name = fname
            
            reader = sitk.ImageFileReader()
            reader.SetFileName(file_name)
            image = reader.Execute()
            
            print(image.GetSpacing())
            np_image = sitk.GetArrayFromImage(image)
            
            self.viewer.add_image([np_image], name='Volume', scale=(float(image.GetSpacing()[2]), float(image.GetSpacing()[1]), float(image.GetSpacing()[0])), 
                      blending='additive', colormap='gray', contrast_limits=self.default_contrast_limits)
    
    def LoadMask(self):
        widget = QtWidgets.QWidget()
        if sys.platform == 'linux':
            fname, _ = QtWidgets.QFileDialog.getOpenFileName(
                widget, 'Load file', '/home/', "Image files (*.tiff *.mha *.nii)")
        else:
            fname, _ = QtWidgets.QFileDialog.getOpenFileName(
                widget, 'Load file', 'c:\\', "Image files (*.tiff *.mha *.nii)")

        if fname != "":
            print(fname)

            file_name = fname
            
            reader = sitk.ImageFileReader()
            reader.SetFileName(file_name)
            image = reader.Execute()
            
            print(image.GetSpacing())
            np_image = sitk.GetArrayFromImage(image)
            np_image = np_image.astype(np.uint8)
            
            self.viewer.add_labels([np_image], name='Mask', scale=(float(image.GetSpacing()[2]), float(image.GetSpacing()[1]), float(image.GetSpacing()[0])))

            
    def SaveSlice(self):
        widget = QtWidgets.QWidget()
        if sys.platform == 'linux':
            fname, _ = QtWidgets.QFileDialog.getSaveFileName(
                widget, 'Save file as', '/home/', "Image files (*.tiff)")
        else:
            fname, _ = QtWidgets.QFileDialog.getSaveFileName(
                widget, 'Save file as', 'c:\\', "Image files (*.tiff)")

        if fname != "":
            print(fname)

            file_name = fname#"/home/tristan/test.tiff"

            output_resolution = float(self.pixel_size.text())
            spacing = (output_resolution, output_resolution, 1)

            volume = []
            # z = self.scroll.value()
            z = self.viewer.dims.current_step[0]
            
            for layer in self.viewer.layers:
                if isinstance(layer, napari.layers.Image):
                    volume.append(layer.data[z,:,:])
            
            volume = np.array(volume)
            
            if False:
                volume = np.expand_dims(volume, axis=1)
                volume = np.moveaxis(volume,0,3)

                volume_itk = sitk.GetImageFromArray(volume)
                volume_itk.SetSpacing(spacing)
                caster = sitk.CastImageFilter()
                caster.SetOutputPixelType(sitk.sitkVectorFloat32)
                volume_itk = caster.Execute(volume_itk)
                sitk.WriteImage(volume_itk, file_name)  
            
            # volume = np.moveaxis(volume,0,2)
            # print(volume.shape)
            tifffile.imsave(file_name, volume)
            
          
        
    
    def SaveSlice2D(self):
        widget = QtWidgets.QWidget()
        if sys.platform == 'linux':
            file_name, _ = QtWidgets.QFileDialog.getSaveFileName(
                widget, 'Save file as', '/home/', "Image files (*.tiff)")
        else:
            file_name, _ = QtWidgets.QFileDialog.getSaveFileName(
                widget, 'Save file as', 'c:\\', "Image files (*.tiff)")

        if file_name != "":
#             print(file_name)
            
            
#             layer_data = []
#             layer_names = []
#             for layer in self.viewer.layers:
#                 if hasattr(layer, 'data'):
#                     data = np.asarray(layer.data[0])
#                     print(data.shape)
#                     layer_data.append(data)
#                     layer_names.append(layer.name)

#             # Create a single array where the first dimension represents the layers
#             if len(layer_data) > 0:
#                 # Create an array with an extra dimension for layers
#                 stacked_data = np.stack(layer_data, axis=0)

#                 # Generate pyramid levels
#                 pyramid_levels = [stacked_data]
#                 for level in range(1, 4):  # Create 3 additional pyramid levels (2x, 4x, 8x downsampling)
#                     downsampled = resize(
#                         stacked_data, 
#                         (stacked_data.shape[0], stacked_data.shape[1] // (2 ** level), stacked_data.shape[2] // (2 ** level)),
#                         anti_aliasing=True
#                     )
#                     pyramid_levels.append(downsampled.astype(stacked_data.dtype))

#                 # Prepare OME metadata
#                 metadata = {
#                     'axes': 'CYX',  # Channel, Y, X dimensions
#                     'Channel': [{'Name': name} for name in layer_names]
#                 }

#                 # Write pyramidal OME-TIFF
#                 with TiffWriter(file_name, ome=True, bigtiff=True) as tif:
#                     options = dict(tile=(256, 256), compression='deflate')
#                     tif.write(pyramid_levels[0], subifds=len(pyramid_levels) - 1, metadata=metadata, **options)
#                     for level in pyramid_levels[1:]:
#                         tif.write(level, subfiletype=1, **options)

#                 print(f"Saved {len(layer_data)} layers to {file_name}")
#             else:
#                 print("No layers to save.")

            
            import numpy as np
            from tifffile import imwrite

            def crop_to_nonzero(data):
                """Crop the 2D array to the minimum bounding box of non-zero elements."""
                rows = np.any(data, axis=1)
                cols = np.any(data, axis=0)
                ymin, ymax = np.where(rows)[0][[0, -1]]
                xmin, xmax = np.where(cols)[0][[0, -1]]
                return data[ymin:ymax+1, xmin:xmax+1]

            layer_data = []
            layer_names = []
            for layer in self.viewer.layers:
                if hasattr(layer, 'data'):
                    data = np.asarray(layer.data[0])
                    # Crop data to non-zero region
                    cropped_data = crop_to_nonzero(data)
                    if cropped_data.size > 0:
                        print(cropped_data.shape)
                        layer_data.append(cropped_data)
                        layer_names.append(layer.name)

            # Create a single array where the first dimension represents the layers
            if len(layer_data) > 0:
                # Stack all cropped layers
                stacked_data = np.stack(layer_data, axis=0)

                # Prepare OME metadata
                metadata = {
                    'axes': 'CYX',  # Channel, Y, X dimensions
                    'Channel': [{'Name': name} for name in layer_names]
                }

                # Write data to OME-TIFF
                imwrite(file_name, stacked_data, photometric='minisblack', metadata=metadata, compression='zlib')
                print(f"Saved {len(layer_data)} layers to {file_name}")
            else:
                print("No layers to save.")




#             output_resolution = float(self.pixel_size.text())
#             spacing = (output_resolution, output_resolution, 1)

#             volume = []
#             z = self.viewer.dims.current_step[0]
            
#             for layer in self.viewer.layers:
#                 if isinstance(layer, napari.layers.Image):
#                     volume.append(layer.data[z,:,:])
            
#             volume = np.array(volume) 
            
#             tifffile.imsave(file_name, volume)
        
        
    def Normalize_slices(self, volume, optical_sections):
        import statistics
        print("Normalize optical sections")

        slices, size_x, size_y = volume.shape
        for i in range(1,slices):
            if i%optical_sections != 0:
                values_x = []
                values_y = []
                for j in range(10000):
                    rand_x = random.randint(1,size_x-1)
                    rand_y = random.randint(1,size_y-1)
                    if(volume[i-1,rand_x,rand_y] > 0 and volume[i,rand_x,rand_y] > 0):
                        x = volume[i,rand_x,rand_y]
                        y = volume[i-1,rand_x,rand_y]
                        values_x.append(x)
                        values_y.append(y)

                if len(values_x) > 3:
                    slope, intercept, r, p, std_err = stats.linregress(values_x, values_y)
                    volume[i,:,:] = volume[i,:,:] * slope + intercept

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

            self.viewer.layers.remove('C1')
            self.viewer.add_image([self.aligned_1], name='C1', scale=(float(self.slice_spacing)/float(self.optical_slices), output_resolution, output_resolution), 
                                  blending='additive', colormap='bop purple', contrast_limits=self.default_contrast_limits)


        if self.cb_C2.isChecked():  
            self.aligned_2[0::10,:,:] = self.aligned_2[0::10,:,:] * norm_value
            self.aligned_2[1::10,:,:] = self.aligned_2[1::10,:,:] * norm_value
            self.aligned_2[2::10,:,:] = self.aligned_2[2::10,:,:] * norm_value
            self.aligned_2[3::10,:,:] = self.aligned_2[3::10,:,:] * norm_value
            self.aligned_2[4::10,:,:] = self.aligned_2[4::10,:,:] * norm_value

            # with napari.gui_qt():
            self.viewer.layers.remove('C2')
            self.viewer.add_image([self.aligned_2], name='C2', scale=(float(self.slice_spacing)/float(self.optical_slices), output_resolution, output_resolution), 
                                  blending='additive', colormap='red', contrast_limits=self.default_contrast_limits)

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

            self.viewer.layers.remove('C3')
            self.viewer.add_image([self.aligned_3], name='C3', scale=(float(self.slice_spacing)/float(self.optical_slices), output_resolution, output_resolution), 
                                  blending='additive', colormap='green', contrast_limits=self.default_contrast_limits)

        if self.cb_C4.isChecked():
            self.aligned_4[0::10,:,:] = self.aligned_4[0::10,:,:] * norm_value
            self.aligned_4[1::10,:,:] = self.aligned_4[1::10,:,:] * norm_value
            self.aligned_4[2::10,:,:] = self.aligned_4[2::10,:,:] * norm_value
            self.aligned_4[3::10,:,:] = self.aligned_4[3::10,:,:] * norm_value
            self.aligned_4[4::10,:,:] = self.aligned_4[4::10,:,:] * norm_value

            self.viewer.layers.remove('C4')
            self.viewer.add_image([self.aligned_4], name='C4', scale=(float(self.slice_spacing)/float(self.optical_slices), output_resolution, output_resolution), 
                                  blending='additive', colormap='bop blue', contrast_limits=self.default_contrast_limits)


    
    def SelectFolder(self):
        self.old_method = False
        file_list = []
        
        folder = str(QtWidgets.QFileDialog.getExistingDirectory())
        print(folder)
        folder = folder
        print(folder)
        self.image_folder = folder
        print(self.image_folder)
        self.search_folder = QtWidgets.QLineEdit(folder)
        if os.path.exists(self.search_folder.text()):
            for f in os.scandir(self.search_folder.text()):
                if f.is_dir():
                    if os.path.exists(f.path + "/mos.zarr"):
                        s = f.path
                        s = s.replace(self.search_folder.text(), '')
                        file_list.append(s)
                    elif os.path.exists(f.path + "/mos"):
                        s = f.path
                        s = s.replace(self.search_folder.text(), '')
                        file_list.append(s)
                            
            print(file_list)

            file_list.sort()

            self.comboBoxPath.clear()
            for i in file_list:
                self.comboBoxPath.addItem(i)
                
                
            if self.cb_R_Axio.isChecked():
                self.old_method = False
                self.axio = True
                self.image_folders = [
                    folder,
                    '/storage/imaxt.processed.2022/axio/',
                    '/storage/processed.2022/axio/',
                ]
                print("setting axio method")
            else:
                self.old_method = False
                self.axio = False
                self.image_folders = [
                    folder,
                    '/storage/processed.2022/stpt/',
                    '/storage/processed/stpt/',
                    '/storage/imaxt.processed.2022/stpt/',
                    '/data/meds1_c/storage/processed0/stpt/',
                    '/storage/imaxt/imaxt_zarr_imc/'
                ]
                print("setting STPT method")   
            
            print(self.image_folders)
        
    def MethodChanged(self):

        if self.cb_R_Axio.isChecked():
            self.old_method = False
            self.axio = True
            self.image_folders = {
                '/storage/imaxt.processed.2022/axio/',
                '/storage/processed.2022/axio/',
            }
            print("setting AXIO method")
        elif self.cb_R_IMC.isChecked():
            self.old_method = False
            self.axio = False
            self.image_folders = {
                '/storage/imaxt/imaxt_zarr_imc/',
                #'/storage/imaxt/msa51/serial sections 144 CM/',
                #'/storage/imaxt/msa51/Lung_serial_sections/zarr/imc/',
                #'/storage/imaxt/msa51/imc-serial-sections/output_zarr/'
            }
            print("setting IMC method")
        elif self.cb_R_Merged.isChecked():
            self.old_method = False
            self.axio = False
            self.image_folders = {
                '/storage/imaxt/imaxt_zarr_merged/'
            }
            print("setting IMC method")
        else:
            self.old_method = False
            self.axio = False
            self.image_folders = {
                '/storage/processed.2022/stpt/',
                '/storage/processed/stpt/',
                '/storage/imaxt.processed.2022/stpt/',
                '/data/meds1_c/storage/processed0/stpt/',
            }
            print("setting STPT method")

        file_list = []
        for folder in self.image_folders:
            print(sys.platform)
            if sys.platform == 'linux':
                self.search_folder = QtWidgets.QLineEdit(folder)
                if os.path.exists(self.search_folder.text()):
                    for f in os.scandir(self.search_folder.text()):
                        if f.is_dir():
                            if os.path.exists(f.path + "/mos.zarr"):
                                s = f.path
                                s = s.replace(self.search_folder.text(), '')
                                file_list.append(s)
                            elif os.path.exists(f.path + "/mos"):
                                s = f.path
                                s = s.replace(self.search_folder.text(), '')
                                file_list.append(s)
                            elif os.path.exists(f.path + "/"):
                                s = f.path
                                s = s.replace(self.search_folder.text(), '')
                                file_list.append(s)
                                
        else:
            self.search_folder = QtWidgets.QLineEdit('N:/stpt/')
        
        print(file_list)

        file_list.sort()

        self.comboBoxPath.clear()
        for i in file_list:
            self.comboBoxPath.addItem(i)

            
    def extract_um_number(self, filename):
        # This regex looks for a pattern where a number is followed by 'x', another number, and 'um'
        # It captures the number right before 'um'
        match = re.search(r'\d+x(\d+)um', filename)
        if match:
            return int(match.group(1))  # Convert the captured group to an integer and return it
        else:
            return 15  # Return 15 as default if the pattern is not found

        
    def on_combobox_changed(self):
        
        for folder in self.image_folders:
            
            self.image_folder = folder
    
            if folder == '/data/meds1_c/storage/processed0/stpt/':
                self.old_method = True
            else:
                self.old_method = False
                
            file_name = folder + str(self.comboBoxPath.currentText()) + '/mos'
            if os.path.exists(file_name):
                self.image_folder = folder
                self.channels_start = 0
                break
                
            file_name = folder + str(self.comboBoxPath.currentText()) + '/mos.zarr'
            if os.path.exists(file_name):
                self.image_folder = folder
                self.channels_start = 1
                break
            
            file_name = folder + str(self.comboBoxPath.currentText()) + '/'
            if os.path.exists(file_name):
                self.image_folder = folder
                self.channels_start = 1
                break
                
                     
        
        slice_spacing = self.extract_um_number(str(self.comboBoxPath.currentText()))
        self.m_slice_spacing.setText(str(slice_spacing))
        
        if os.path.exists(file_name + '/.zmetadata'):
            print("metadata is available")
            try:
                # Extract the channels list
                channels = data["metadata"][".zattrs"]["meta"][0]["acquisitions"][0]["acquisition_rois"][0]["acquisitions"][0]["channels"]

                # Extract the target values
                targets = [channel["target"] for channel in channels]
                
                print(targets)
                pass
                
            except Exception:
                print("none-consolidated")
                ds = xr.open_zarr(file_name)
                try:
                    def remove_trailing_zeros(number):
                        return str(number).rstrip('0').rstrip('.') if '.' in str(number) else str(number)


                    slice_thickness = remove_trailing_zeros(ds['S001'].attrs['thickness'])
                    self.m_slice_spacing.setText(str(slice_thickness))
                    print(slice_thickness)
                except:
                    pass

                channel_names = []
                try:
                    channel_names = ds.coords['channel'].values.tolist()
                    print(f"channel_names: {channel_names}")
                except:
                    pass

                self.comboBox.clearItems()
                for i in range(len(channel_names)):
                    self.comboBox.addItem(f"{channel_names[i]}")


                selected_channels = parse_channel_input(self.selected_slices.text())
                for i in range(50):
                    if i in selected_channels:
                        self.comboBox.checkItem(i, True)
                    else:
                        self.comboBox.checkItem(i, False)

                self.ignore_gui_call = True;
                self.scroll.setValue(0)
                self.ignore_gui_call = False;
                self.image_slice.setText("0")
                #self.slice_names = ds.attrs['cube_reg']['slice']
                #self.scroll.setRange(0, len(self.slice_names))

                pass
        else:
            print("not changing number of slices")
            
            
    def SetPerspective(self):
        if(self.cb_perspective.isChecked()):
            self.viewer.window.qt_viewer.camera._3D_camera.fov = 45
        else:
            self.viewer.window.qt_viewer.camera._3D_camera.fov = 0
            
    def MergeLayers(self):
        f1 = float(self.m_volume_1_multiplier.text())
        f2 = float(self.m_volume_2_multiplier.text())
        C1 = self.viewer.layers[self.m_volume_1.text()].data
        C2 = self.viewer.layers[self.m_volume_2.text()].data
        C_new = (f1*C1)+(f2*C2)
        self.viewer.add_image(C_new, name=self.m_volume_new.text(), scale=self.viewer.layers[self.m_volume_1.text()].scale, blending='additive', colormap='gray')
        


    def SetShapeText(self):
        print(self.viewer.layers.selection)
        
        shapes = self.viewer.add_points(properties={'box_label': "test"})
        
        shapes.text = 'box_label'
        shapes.text.color = 'white'
        shapes.text.size = 30
        shapes.anchor = 'center'
        shapes.face_color = "transparent"
        shapes.edge_color = "transparent"
        shapes.current_face_color = "transparent"
        shapes.current_edge_color = "transparent"
        shapes.edge_width = 0
        shapes.current_edge_width = 0
        shapes.out_of_slice_display = True

        def on_data(event):
            shapes.text = 'box_label'
            shapes.text.color = 'white'
            shapes.text.size = 30
            shapes.anchor = 'center'
            shapes.face_color = "transparent"
            shapes.edge_color = "transparent"
            shapes.current_face_color = "transparent"
            shapes.current_edge_color = "transparent"
            shapes.edge_width = 0
            shapes.current_edge_width = 0
            shapes.out_of_slice_display = True
            
        shapes.events.set_data.connect(on_data)
        
    def set_overall_brightness(self):
        
        if self.selected_channels is None:
            return
        
        number_of_channels = len(self.selected_channels)
        
        if self.cb_R_IMC.isChecked():
            self.overall_brightness = 0.1 * number_of_channels * (1.0001 - (float(self.scroll_overall_brightness.value()) / 1000))
        else:
            self.overall_brightness = 0.5 * number_of_channels * (1.0001 - (float(self.scroll_overall_brightness.value()) / 1000))
        
        contrast_limits = self.default_contrast_limits
        contrast_limits = [self.value_range[0],self.value_range[1]*self.overall_brightness]
        for layer in self.viewer.layers:
            layer.contrast_limits = contrast_limits

    def OnLoad2D(self):
        
        self.loaded_2D = True
        self.loaded_3D = False
        
        self.bLoad2D.setChecked(True)
        self.bLoad3D.setChecked(False)
        
        # Remove all layers (images)
        self.viewer.layers.select_all()  # Select all layers
        self.viewer.layers.remove_selected()  # Remove selected layers (which are all layers)
        
        
        self.selected_channels = self.comboBox.getCheckedItems()        
        
        
        self.optical_slices_available, self.value_range, self.number_of_sections, channel_names, colors = self.data.Load2D(self.viewer, self.image_folder, self.comboBoxPath.currentText(), self.selected_channels, self.default_contrast_limits, self.thresholdN, self.channels_start, self.axio, self.old_method, self.overall_brightness, self.scroll, self.scroll_overall_brightness)

        number_of_channels = len(self.selected_channels)
        
        print(f"number_of_channels {number_of_channels}")
        
        if self.cb_R_IMC.isChecked():
            self.overall_brightness = 0.1 * number_of_channels * (1.0001 - (float(self.scroll_overall_brightness.value()) / 1000))
        else:
            self.overall_brightness = 0.5 * number_of_channels * (1.0001 - (float(self.scroll_overall_brightness.value()) / 1000))
        
        contrast_limits = self.default_contrast_limits
        contrast_limits = [self.value_range[0],self.value_range[1]*self.overall_brightness]
        for layer in self.viewer.layers:
            layer.contrast_limits = contrast_limits
            
        self.ignore_gui_call = True
        self.scroll.setRange(1, self.number_of_sections * self.optical_slices_available);
        self.ignore_gui_call = False
        
        self.viewer.scale_bar.unit = "um"
        
        #print(channel_names)
        #print(colors)
        #channel_names = ["Channel 1", "Channel 2"]  # Define your channel names
        #colors = [(128, 0, 128), (128, 0, 128)]
        
        self.legend_widget.populate_legend(channel_names, colors)
            

    def onSelectedSlicesTextChanged(self):
        selected_channels = parse_channel_input(self.selected_slices.text())
        for i in range(50):
            if i in selected_channels:
                self.comboBox.checkItem(i, True)
            else:
                self.comboBox.checkItem(i, False)


    def Load3D(self, text):

        self.origin_x = 0
        self.origin_y = 0

        self.crop = False
        self.loaded_2D = False
        self.loaded_3D = True
        
        self.bLoad2D.setChecked(False)
        self.bLoad3D.setChecked(True)
        
        # Remove all layers (images)
        self.viewer.layers.select_all()  # Select all layers
        self.viewer.layers.remove_selected()  # Remove selected layers (which are all layers)
        
        self.selected_channels = self.comboBox.getCheckedItems()  
        
        self.optical_slices_available, self.value_range, self.shape, self.slice_spacing, self.optical_slices, self.output_resolution = self.data.Load3D(self.viewer, self.image_folder, self.comboBoxPath.currentText(), self.selected_channels, self.default_contrast_limits, self.thresholdN, self.channels_start, self.axio, self.old_method, self.overall_brightness, self.scroll, self.scroll_overall_brightness, self.pixel_size.text(), self.m_slice_spacing.text(), self.start_slice.text(), self.end_slice.text(), self.crop, self.crop_start_x, self.crop_end_x, self.crop_start_y, self.crop_end_y)
        
        self.MakeBoundingBox()
        
        
        number_of_channels = len(self.selected_channels)
        
        if self.cb_R_IMC.isChecked():
            self.overall_brightness = 0.1 * number_of_channels * (1.0001 - (float(self.scroll_overall_brightness.value()) / 1000))
        else:
            self.overall_brightness = 0.5 * number_of_channels * (1.0001 - (float(self.scroll_overall_brightness.value()) / 1000))
        
        contrast_limits = self.default_contrast_limits
        contrast_limits = [self.value_range[0],self.value_range[1]*self.overall_brightness]
        for layer in self.viewer.layers:
            layer.contrast_limits = contrast_limits
            
        self.viewer.scale_bar.unit = "um"
            
            
    def on_visibility_changed(self):
        colors = []
        
        channel_names_used = []
        
        for layer in self.viewer.layers:
            if isinstance(layer, napari.layers.Image):
                far_end_color = layer.colormap.map([1.0])[0]
                rgb_color = far_end_color[:3]
                color_map_tuple = tuple(int(c * 255) for c in rgb_color)
                colors.append(color_map_tuple)
                channel_names_used.append(layer.name)
        
        self.legend_widget.populate_legend(channel_names_used, colors)
      
            
            
            
        
    def main(self):

        self.viewer = napari.Viewer()
        self.viewer.theme = 'dark' 

        widget = QtWidgets.QWidget()
        layout = QtWidgets.QVBoxLayout()
        vbox = QtWidgets.QVBoxLayout()

        cb_group1 = QtWidgets.QButtonGroup()
        hbox = QtWidgets.QHBoxLayout()
        hbox.addStretch(1)
        
        self.cb_R_Axio = QtWidgets.QRadioButton('AXIO')
        self.cb_R_Axio.setChecked(False)
        self.cb_R_Axio.toggled.connect(self.MethodChanged)
        hbox.addWidget(self.cb_R_Axio)
        
        self.cb_R_STPT = QtWidgets.QRadioButton('STPT')
        self.cb_R_STPT.setChecked(True)
        self.cb_R_STPT.toggled.connect(self.MethodChanged)
        hbox.addWidget(self.cb_R_STPT)
        
        self.cb_R_IMC= QtWidgets.QRadioButton('IMC')
        self.cb_R_IMC.setChecked(False)
        self.cb_R_IMC.toggled.connect(self.MethodChanged)
        hbox.addWidget(self.cb_R_IMC)
        
        self.cb_R_Merged= QtWidgets.QRadioButton('Merged')
        self.cb_R_Merged.setChecked(False)
        self.cb_R_Merged.toggled.connect(self.MethodChanged)
        hbox.addWidget(self.cb_R_Merged)
        
        
        cb_group1.addButton(self.cb_R_Axio)
        cb_group1.addButton(self.cb_R_STPT)
        cb_group1.addButton(self.cb_R_IMC)
        cb_group1.addButton(self.cb_R_Merged)
        
        bSelectFolder = QtWidgets.QPushButton('Select folder')
        bSelectFolder.clicked.connect(self.SelectFolder)
        hbox.addWidget(bSelectFolder)
        
        hbox.addStretch(1)
        vbox.addLayout(hbox)
        
        
        
        
        hbox = QtWidgets.QHBoxLayout()
        self.comboBoxPath = ExtendedComboBox()

        self.image_folders = {
            '/storage/processed.2022/stpt/',
            '/storage/processed/stpt/',
            '/storage/imaxt.processed.2022/stpt/',
            '/data/meds1_c/storage/processed0/stpt/',
            #'/storage/imaxt/imaxt_zarr_imc/'
            '/storage/imaxt/eglez/processed/lightsheet/',
        }
        file_list = []
        for folder in self.image_folders:
            self.search_folder = QtWidgets.QLineEdit(folder)
            if os.path.exists(self.search_folder.text()):
                for f in os.scandir(self.search_folder.text()):
                    if f.is_dir():
                        if os.path.exists(f.path + "/mos.zarr"):
                            s = f.path
                            s = s.replace(self.search_folder.text(), '')
                            file_list.append(s)
                        elif os.path.exists(f.path + "/mos"):
                            s = f.path
                            s = s.replace(self.search_folder.text(), '')
                            file_list.append(s)

        
        file_list.sort()
        
        self.axio = False

        self.comboBoxPath.clear()
        for i in file_list:
            self.comboBoxPath.addItem(i)

        self.comboBoxPath.setMinimumWidth(300)
        self.comboBoxPath.currentIndexChanged.connect(
            self.on_combobox_changed)
        hbox.addWidget(self.comboBoxPath)
        self.comboBoxPath.setMaximumWidth(300)
        vbox.addLayout(hbox)
        
        layout.addLayout(vbox)
        
        vbox = QtWidgets.QVBoxLayout()

        hbox = QtWidgets.QHBoxLayout()

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

        # Create the label and set its size policy
        label = QtWidgets.QLabel("Brightness:")
        label.setSizePolicy(QtWidgets.QSizePolicy.Preferred, QtWidgets.QSizePolicy.Fixed)
        hbox.addWidget(label)

        # Create the scrollbar
        self.scroll_overall_brightness = QtWidgets.QScrollBar()
        self.scroll_overall_brightness.setOrientation(QtCore.Qt.Horizontal)
        self.scroll_overall_brightness.setRange(0, 1000)
        self.scroll_overall_brightness.setValue(800)
        self.scroll_overall_brightness.setMinimumWidth(150)
        self.scroll_overall_brightness.valueChanged.connect(self.set_overall_brightness)

        # Set the size policy for the scrollbar to expanding
        self.scroll_overall_brightness.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)

        # Add the scrollbar to the layout
        hbox.addWidget(self.scroll_overall_brightness)

        vbox.addLayout(hbox)
        

        hbox = QtWidgets.QHBoxLayout()
        hbox.addWidget(QtWidgets.QLabel("Slices range:"))
        self.start_slice = QtWidgets.QLineEdit("")
        self.start_slice.setMaximumWidth(50)
        hbox.addWidget(self.start_slice)
        hbox.addWidget(QtWidgets.QLabel("to"))
        self.end_slice = QtWidgets.QLineEdit("")
        self.end_slice.setMaximumWidth(50)
        hbox.addWidget(self.end_slice)
        hbox.addStretch(1)
        vbox.addLayout(hbox)
        
        
        hbox = QtWidgets.QHBoxLayout()
        hbox.addWidget(QtWidgets.QLabel("Channels:"))
        self.selected_slices = QtWidgets.QLineEdit("0-2,3,4")
        #self.selected_slices.setMinimumWidth(50)
        self.selected_slices.setMaximumWidth(100)
        self.selected_slices.textChanged.connect(self.onSelectedSlicesTextChanged)
        hbox.addWidget(self.selected_slices)
        #hbox.addStretch(1)
        
        
        
        self.comboBox = CustomComboBox()
        # for i in range(100):
        #     self.comboBox.addItem(f"Channel {i+1}")
        self.comboBox.setSizePolicy(QtWidgets.QSizePolicy.Expanding, QtWidgets.QSizePolicy.Fixed)
        hbox.addWidget(self.comboBox)
        
        vbox.addLayout(hbox)
        

        hbox = QtWidgets.QHBoxLayout()
        self.bLoad3D = QtWidgets.QPushButton('Load volume')
        self.bLoad3D.setCheckable(True)
        self.bLoad3D.clicked.connect(self.Load3D)
        hbox.addWidget(self.bLoad3D)

        bReload3D = QtWidgets.QPushButton('Reload in shape')
        #bReload3D.setCheckable(True)
        bReload3D.clicked.connect(self.LoadInRegion)
        hbox.addWidget(bReload3D)
        bCrop = QtWidgets.QPushButton('Crop to shape')
        #bCrop.setCheckable(True)
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
        self.maxSizeN = QtWidgets.QLineEdit("1000")
        hbox.addWidget(self.maxSizeN)
        hbox.addStretch(1)
        vbox.addLayout(hbox)

        hbox = QtWidgets.QHBoxLayout()
        bKeepN = QtWidgets.QPushButton('Show only large regions')
        #bKeepN.setCheckable(True)
        bKeepN.clicked.connect(self.Keep_n_Regions)
        hbox.addWidget(bKeepN)
        bRemoveN = QtWidgets.QPushButton('Remove small regions')
        #bRemoveN.setCheckable(True)
        bRemoveN.clicked.connect(self.Remove_Small_Regions)
        hbox.addWidget(bRemoveN)
        hbox.addStretch(1)
        vbox.addLayout(hbox)
        
        
        
        
        hbox = QtWidgets.QHBoxLayout()

        bSaveSlice = QtWidgets.QPushButton('Save slice')
        bSaveSlice.clicked.connect(self.SaveSlice)
        hbox.addWidget(bSaveSlice)

        bSaveVolume = QtWidgets.QPushButton('Save volume')
        bSaveVolume.clicked.connect(self.SaveVolume)
        hbox.addWidget(bSaveVolume)

        bLoadVolume = QtWidgets.QPushButton('Load volume')
        bLoadVolume.clicked.connect(self.LoadVolume)
        hbox.addWidget(bLoadVolume)

        bLoadMask = QtWidgets.QPushButton('Load mask')
        bLoadMask.clicked.connect(self.LoadMask)
        hbox.addWidget(bLoadMask)

        hbox.addStretch(1)
        vbox.addLayout(hbox)
        
        bAddPolygon = QtWidgets.QPushButton('Add')
        #bAddPolygon.setCheckable(True)
        bAddPolygon.clicked.connect(self.add_polygon)
        bRemoveOutside = QtWidgets.QPushButton('Remove outside interpolated region')
        #bRemoveOutside.setCheckable(True)
        bRemoveOutside.clicked.connect(self.run_remove_outside)
        hbox = QtWidgets.QHBoxLayout()
        hbox.addWidget(bAddPolygon)
        hbox.addWidget(bRemoveOutside)
        hbox.addStretch(1)
        vbox.addLayout(hbox)


        bAdd3DShape = QtWidgets.QPushButton('2D to 3D shape')
        #bAdd3DShape.setCheckable(True)
        bAdd3DShape.clicked.connect(self.Make3DShape)
        hbox = QtWidgets.QHBoxLayout()
        hbox.addWidget(bAdd3DShape)

        hbox.addStretch(1)
        vbox.addLayout(hbox)
        
        hbox = QtWidgets.QHBoxLayout()
        cb_group2 = QtWidgets.QButtonGroup()
        self.cb_perspective = QtWidgets.QRadioButton('Perspective')
        self.cb_perspective.setChecked(True)
        self.cb_perspective.toggled.connect(self.SetPerspective)
        hbox.addWidget(self.cb_perspective)
        self.cb_isometric = QtWidgets.QRadioButton('Isometric')
        self.cb_isometric.setChecked(False)
        self.cb_isometric.toggled.connect(self.SetPerspective)
        hbox.addWidget(self.cb_isometric)
        cb_group2.addButton(self.cb_perspective)
        cb_group2.addButton(self.cb_isometric)
        hbox.addStretch(1)

        vbox.addLayout(hbox)
        vbox.addStretch(1)
        
        
        
        
        
        groupbox = QtWidgets.QGroupBox("3D")
        groupbox.setCheckable(False)
        groupbox.setLayout(vbox)
        layout.addWidget(groupbox)
        #widget.setLayout(layout)
        
        
        
        
        
        

        vbox = QtWidgets.QVBoxLayout()

        
        hbox = QtWidgets.QHBoxLayout()
        self.bLoad2D = QtWidgets.QPushButton('Load slice')
        self.bLoad2D.clicked.connect(self.OnLoad2D)
        hbox.addWidget(self.bLoad2D)
        
        
        bSaveSlice = QtWidgets.QPushButton('Save slice')
        bSaveSlice.clicked.connect(self.SaveSlice2D)
        hbox.addWidget(bSaveSlice)
        

        hbox.addStretch(1)
        vbox.addLayout(hbox)
        

        hbox = QtWidgets.QHBoxLayout()
        #hbox.addWidget(QtWidgets.QLabel("Slice:"))
        self.scroll = QtWidgets.QScrollBar()
        self.scroll.setOrientation(QtCore.Qt.Horizontal)
        self.scroll.setRange(1, 100)
        self.scroll.setMinimumWidth(150)
        self.scroll.valueChanged.connect(self.set_image_slice_value)
        hbox.addWidget(self.scroll)

        self.image_slice = QtWidgets.QLineEdit("0")
        self.image_slice.setMaximumWidth(30)
        hbox.addWidget(self.image_slice)

        #hbox.addStretch(1)
        vbox.addLayout(hbox)
        
        
        
        groupbox = QtWidgets.QGroupBox("2D")
        groupbox.setCheckable(False)
        groupbox.setLayout(vbox)
        layout.addWidget(groupbox)
        
        
        widget.setLayout(layout)
        
        
        


        #widget.setLayout(vbox)
        
        dw1 = self.viewer.window.add_dock_widget(widget, area="right")
        dw1.setWindowTitle('Main')
        
        self.comboBoxPath.setMaximumWidth(1000000)
        
        
        
        
        # Math widget
        widget_merge = QtWidgets.QWidget()
        
        vbox = QtWidgets.QVBoxLayout()
        
        hbox = QtWidgets.QHBoxLayout()
        hbox.addWidget(QtWidgets.QLabel("Layer 1:"))
        self.m_volume_1 = QtWidgets.QLineEdit("C2")
        self.m_volume_1.setMaximumWidth(75)
        hbox.addWidget(self.m_volume_1)
        hbox.addWidget(QtWidgets.QLabel("Multiplier:"))
        self.m_volume_1_multiplier = QtWidgets.QLineEdit("-5")
        self.m_volume_1_multiplier.setMaximumWidth(50)
        hbox.addWidget(self.m_volume_1_multiplier)
        hbox.addStretch(1)
        vbox.addLayout(hbox)
        
        hbox = QtWidgets.QHBoxLayout()
        hbox.addWidget(QtWidgets.QLabel("Layer 2:"))
        self.m_volume_2 = QtWidgets.QLineEdit("C3")
        self.m_volume_2.setMaximumWidth(75)
        hbox.addWidget(self.m_volume_2)
        hbox.addWidget(QtWidgets.QLabel("Multiplier:"))
        self.m_volume_2_multiplier = QtWidgets.QLineEdit("1")
        self.m_volume_2_multiplier.setMaximumWidth(50)
        hbox.addWidget(self.m_volume_2_multiplier)
        hbox.addStretch(1)
        vbox.addLayout(hbox)
        
        hbox = QtWidgets.QHBoxLayout()
        bMerge = QtWidgets.QPushButton('Merge layers')
        #bMerge.setCheckable(False)
        bMerge.clicked.connect(self.MergeLayers)
        hbox.addWidget(bMerge)
        hbox.addWidget(QtWidgets.QLabel("Output layer:"))
        self.m_volume_new = QtWidgets.QLineEdit("C_new")
        self.m_volume_new.setMaximumWidth(75)
        hbox.addWidget(self.m_volume_new)
        
        hbox.addStretch(1)
        vbox.addLayout(hbox)
        
        
        vbox.addStretch(1)
        widget_merge.setLayout(vbox)
        
        dw2 = self.viewer.window.add_dock_widget(widget_merge, area='right')
        dw2.setWindowTitle('Processing')
        
        
        # Create and add the legend widget
        self.legend_widget = LegendWidget()
        self.legend_widget.setAutoFillBackground(True)
        p = self.legend_widget.palette()
        p.setColor(self.legend_widget.backgroundRole(), Qt.red)
        self.legend_widget.setPalette(p)
        
        dw4 = self.viewer.window.add_dock_widget(self.legend_widget, area='right')
        dw4.setWindowTitle('Legend')
        
        # Connect the visibilityChanged signal to the on_visibility_changed function
        dw4.visibilityChanged.connect(self.on_visibility_changed)
        
        # Animation widget
        animation_widget = AnimationWidget(self.viewer)
        
        dw3 = self.viewer.window.add_dock_widget(animation_widget, area='right')
        dw3.setWindowTitle('Animation')
        
        
        self.viewer.window._qt_window.tabifyDockWidget(dw1, dw2)
        self.viewer.window._qt_window.tabifyDockWidget(dw1, dw3)
        self.viewer.window._qt_window.tabifyDockWidget(dw1, dw4)
        
        
        
        
        


        
        
        

        napari.run()