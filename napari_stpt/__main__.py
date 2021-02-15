"""
napari-STPT reads zarr files and displays them
"""

#from skimage import data, measure
#from skimage.transform import pyramid_gaussian
#from skimage.filters import threshold_otsu, threshold_local
#from skimage.segmentation import clear_border
#from skimage.measure import label
#from skimage.morphology import closing, square, remove_small_objects, binary_opening, skeletonize

#from scipy import ndimage, stats

import os
import sys
#import math
#import cv2
import dask
import napari
#import trimesh
#import statistics
import numpy as np
import xarray as xr
#from PIL import Image
#import openmesh as om
#import raster_geometry as rg
#import matplotlib.pyplot as plt
#from matplotlib.backends.backend_qt5agg import FigureCanvas
#from matplotlib.figure import Figure

#from PyQt5.QtWidgets import QFileDialog, QWidget, QPushButton, QApplication, QHBoxLayout, QGroupBox, QVBoxLayout, QGridLayout, QSpacerItem, QLineEdit, QLabel, QComboBox, QCheckBox

from qtpy import QtCore, QtGui, QtWidgets
#from PySide2.QtCore import Qt

import SimpleITK as sitk

#import itk
#from distutils.version import StrictVersion as VS
#if VS(itk.Version.GetITKVersion()) < VS("5.0.0"):
#    print("ITK 5.0.0 or newer is required.")
#    sys.exit(1)

def Load():

    global ds1
    global ds
    global bscale
    global bzero
    global comboBoxPath
    global align_x
    global align_y
    global viewer
    global volume_1
    global volume_2
    global volume_3
    global volume_4
    global comboBoxResolution

    print('N:/stpt/' + str(comboBoxPath.currentText())+'/mos.zarr')

    ds1 = xr.open_zarr('N:/stpt/' + str(comboBoxPath.currentText())+'/mos.zarr')
    #ds2 = xr.open_zarr(str(comboBoxPath.currentText())+'/mos.zarr', group='l.2')
    #ds4 = xr.open_zarr(str(comboBoxPath.currentText())+'/mos.zarr', group='l.4')
    #ds8 = xr.open_zarr(str(comboBoxPath.currentText())+'/mos.zarr', group='l.8')
    #ds32 = xr.open_zarr(str(comboBoxPath.currentText())+'/mos.zarr', group='l.32')

    
    bscale = ds1.attrs['bscale']
    bzero = ds1.attrs['bzero']

    #print(ds1.attrs.keys())
    #print(ds1['S001'].attrs.keys())
    #print(ds1['S001'].attrs)
    align_x = ds1.attrs['cube_reg']['abs_dx']
    align_y = ds1.attrs['cube_reg']['abs_dy']

    print("size: {}".format(ds1.dims))
    if str(comboBoxResolution.currentText()) == '1':
        print("loading at resolution 1")
        ds = ds1
        resolution = 1
    elif str(comboBoxResolution.currentText()) == '2':
        print("loading at resolution 2")
        ds = xr.open_zarr('N:/stpt/' + str(comboBoxPath.currentText())+'/mos.zarr', group='l.2')
        resolution = 2
    elif str(comboBoxResolution.currentText()) == '4':
        print("loading at resolution 4")
        ds = xr.open_zarr('N:/stpt/' + str(comboBoxPath.currentText())+'/mos.zarr', group='l.4')
        resolution = 4
    elif str(comboBoxResolution.currentText()) == '8':
        print("loading at resolution 8")
        ds = xr.open_zarr('N:/stpt/' + str(comboBoxPath.currentText())+'/mos.zarr', group='l.8')
        resolution = 8
    elif str(comboBoxResolution.currentText()) == '16':
        print("loading at resolution 16")
        ds = xr.open_zarr('N:/stpt/' + str(comboBoxPath.currentText())+'/mos.zarr', group='l.16')
        resolution = 16
    else:
        print("loading at resolution 32")
        ds = xr.open_zarr('N:/stpt/' + str(comboBoxPath.currentText())+'/mos.zarr', group='l.32')
        resolution = 32

    global pixel_size
    output_resolution = float(pixel_size.text())

    with napari.gui_qt():
        if cb_C1.isChecked():
            volume_1 = (ds.sel(channel=1, type='mosaic', z=0).to_array().data * bscale + bzero)
            aligned_1 = Align(volume_1, resolution, output_resolution)
            viewer.add_image([aligned_1], name='C1', scale=(15,output_resolution,output_resolution), blending='additive', colormap='red')
        
        if cb_C2.isChecked():
            volume_2 = (ds.sel(channel=2, type='mosaic', z=0).to_array().data * bscale + bzero)
            aligned_2 = Align(volume_2, resolution, output_resolution)
            viewer.add_image([aligned_2], name='C2', scale=(15,output_resolution,output_resolution), blending='additive', colormap='green')

        if cb_C3.isChecked():
            volume_3 = (ds.sel(channel=3, type='mosaic', z=0).to_array().data * bscale + bzero)
            aligned_3 = Align(volume_3, resolution, output_resolution)
            viewer.add_image([aligned_3], name='C3', scale=(15,output_resolution,output_resolution), blending='additive', colormap='blue')

        if cb_C4.isChecked():
            volume_4 = (ds.sel(channel=4, type='mosaic', z=0).to_array().data * bscale + bzero)
            aligned_4 = Align(volume_4, resolution, output_resolution)
            viewer.add_image([aligned_4], name='C4', scale=(15,output_resolution,output_resolution), blending='additive', colormap='bop purple')



    
def Align(volume, resolution, output_resolution):

    spacing = (ds1['S001'].attrs['scale'])
    #size_multiplier = (resolution*0.1*spacing[0])/7.5
    #print(resolution*0.1*spacing[0])
    size_multiplier = (resolution*0.1*spacing[0])/output_resolution
    size = (volume.shape[0], int(size_multiplier*volume.shape[1]), int(size_multiplier*volume.shape[2]))
    aligned = np.zeros(size, dtype=float)
    size2D = ( int(size_multiplier*volume.shape[2]), int(size_multiplier*volume.shape[1]))

    global align_x
    global align_y

    z_size = volume.shape[0]
    for z in range(0, z_size):
        print('{}/{}'.format(z,z_size))
        fixed = sitk.GetImageFromArray(volume[z,:,:])
        fixed.SetOrigin((0, 0))

        slice_name = 'S00'+str(z+1)
        if(z+1 > 9):
            slice_name = 'S0'+str(z+1)
        if(z+1 > 99):
            slice_name = 'S'+str(z+1)
        spacing = (ds1[slice_name].attrs['scale'])
        #print(spacing)
        #fixed.SetSpacing([1.6*spacing[0],1.6*spacing[1]])
        fixed.SetSpacing([resolution*0.1*spacing[1],resolution*0.1*spacing[0]])

        transform = sitk.Euler2DTransform()

        alignY = 0
        if not np.isnan(align_y[z]):
            alignY = -align_y[z]*0.1*spacing[1]
            
        alignX = 0
        if not np.isnan(align_x[z]):
            alignX = -align_x[z]*0.1*spacing[0] 

        transform.SetTranslation([alignY,alignX])
        
        #transform.SetTranslation([-align_y[z]*0.1*spacing[1],-align_x[z]*0.1*spacing[0]])
        #transform.SetTranslation([0.5*-align_y[z],0.5*-align_x[z]])
        #transform.SetTranslation([0,0])
        print(align_y[z])

        resampler = sitk.ResampleImageFilter()
        #resampler.SetReferenceImage(fixed)
        
        resampler.SetSize(size2D)
        resampler.SetOutputSpacing([output_resolution,output_resolution])
        resampler.SetOutputOrigin((0, 0))
        resampler.SetInterpolator(sitk.sitkLinear)
        resampler.SetDefaultPixelValue(0)
        resampler.SetTransform(transform)

        out = resampler.Execute(fixed)

        np_out = sitk.GetArrayFromImage(out)
        aligned[z,:,:] = np_out

    return aligned

def SelectFolder():
    global search_folder
    global comboBoxPath
    file = str(QFileDialog.getExistingDirectory())
    search_folder.setText(file)
    comboBoxPath.clear()
    for f in os.scandir(search_folder.text()):
         if f.is_dir():
            s = f.path
            s = s.replace(search_folder.text(), '')
            print(s)
            comboBoxPath.addItem(s)


def main():
    with napari.gui_qt():

        global viewer
        global comboBoxPath
        global comboBoxResolution
        global cb_C1
        global cb_C2
        global cb_C3
        global cb_C4

        loaded16_2 = 0
        loaded16_3 = 0
        loaded16_4 = 0

        viewer = napari.Viewer()

        widget = QtWidgets.QWidget()
        vbox = QtWidgets.QVBoxLayout()


        hbox = QtWidgets.QHBoxLayout()
        hbox.addWidget(QtWidgets.QLabel("Data folder:"))
        global search_folder
        search_folder = QtWidgets.QLineEdit('N:/stpt/')
        search_folder.setMaximumWidth(250)
        hbox.addWidget(search_folder)
        bSelectFolder = QtWidgets.QPushButton('Set folder')
        bSelectFolder.setCheckable(True)
        bSelectFolder.clicked.connect(SelectFolder)
        hbox.addWidget(bSelectFolder)
        #hbox.addStretch(1)
        vbox.addLayout(hbox)




        hbox = QtWidgets.QHBoxLayout()
        global comboBoxPath
        comboBoxPath = QtWidgets.QComboBox()    
        #for f in os.scandir('K:/STPT/'):
        #     if f.is_dir():
        #        print(f.path)
        #        comboBoxPath.addItem(f.path)
        if os.path.exists(search_folder.text()):
            for f in os.scandir(search_folder.text()):
                if f.is_dir():
                    s = f.path
                    s = s.replace(search_folder.text(), '')
                    print(s)
                    comboBoxPath.addItem(s)
        comboBoxPath.removeItem(0)
        comboBoxPath.setMaximumWidth(300)
        hbox.addWidget(comboBoxPath)
        hbox.addStretch(1)
        vbox.addLayout(hbox)

        
        hbox = QtWidgets.QHBoxLayout()
        hbox.addWidget(QtWidgets.QLabel("Resolution:"))
        comboBoxResolution = QtWidgets.QComboBox() 
        comboBoxResolution.addItem(str(1))
        comboBoxResolution.addItem(str(2))
        comboBoxResolution.addItem(str(4))
        comboBoxResolution.addItem(str(8))
        comboBoxResolution.addItem(str(16))
        comboBoxResolution.addItem(str(32))
        comboBoxResolution.setCurrentText(str(16))
        comboBoxResolution.setMaximumWidth(100)
        hbox.addWidget(comboBoxResolution)
        hbox.addStretch(1)

        hbox.addWidget(QtWidgets.QLabel("Output pixel size:"))
        global pixel_size
        pixel_size = QtWidgets.QLineEdit("7.5")
        pixel_size.setMaximumWidth(50)
        hbox.addWidget(pixel_size)
        hbox.addStretch(1)
        vbox.addLayout(hbox)

            
        hbox = QtWidgets.QHBoxLayout()
        hbox.addWidget(QtWidgets.QLabel("Load channels:"))
        cb_C1 = QtWidgets.QCheckBox('1')
        cb_C1.setChecked(True)
        hbox.addWidget(cb_C1)
        cb_C2 = QtWidgets.QCheckBox('2')
        cb_C2.setChecked(True)
        hbox.addWidget(cb_C2)
        cb_C3 = QtWidgets.QCheckBox('3')
        cb_C3.setChecked(True)
        hbox.addWidget(cb_C3)
        cb_C4 = QtWidgets.QCheckBox('4')
        cb_C4.setChecked(True)
        hbox.addWidget(cb_C4)
        hbox.addStretch(1)
        vbox.addLayout(hbox)


        bLoad = QtWidgets.QPushButton('Load image')
        bLoad.setCheckable(True)
        bLoad.clicked.connect(Load)
        hbox = QtWidgets.QHBoxLayout()
        hbox.addWidget(bLoad)
        hbox.addStretch(1)
        vbox.addLayout(hbox)
        
        vbox.addStretch(1)
        widget.setLayout(vbox)
        viewer.window.add_dock_widget(widget, area="right" )



if __name__ == "__main__":
    main()
