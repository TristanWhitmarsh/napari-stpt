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

class NapariSTPT(object):

    def Load(self, text):
        global ds1
        global ds
        global bscale
        global bzero
        global comboBoxPath
        global align_x
        global align_y
        global viewer
        global comboBoxResolution
        global search_folder
        global pixel_size
        global spacing

        print(search_folder.text() + str(comboBoxPath.currentText()) + '/mos.zarr')
        

        try:
            ds1 = xr.open_zarr(search_folder.text() + str(comboBoxPath.currentText())+'/mos.zarr', consolidated = True)
            print("consolidated = True") 
        except:
            print("An exception occurred") 
            ds1 = xr.open_zarr(search_folder.text() + str(comboBoxPath.currentText())+'/mos.zarr')
        #ds2 = xr.open_zarr(str(comboBoxPath.currentText())+'/mos.zarr', group='l.2')
        #ds4 = xr.open_zarr(str(comboBoxPath.currentText())+'/mos.zarr', group='l.4')
        #ds8 = xr.open_zarr(str(comboBoxPath.currentText())+'/mos.zarr', group='l.8')
        #ds32 = xr.open_zarr(str(comboBoxPath.currentText())+'/mos.zarr', group='l.32')

        
        bscale = ds1.attrs['bscale']
        bzero = ds1.attrs['bzero']
        
        #print(ds1.attrs.keys())
        #print(ds1['S001'].attrs.keys())
        #print(ds1['S001'].attrs)
        #print(ds1.attrs["multiscale"]['datasets'][1]['level'])
        #print(ds1.attrs["multiscale"]['datasets'][1]['path'])
        
        align_x = ds1.attrs['cube_reg']['abs_dx']
        align_y = ds1.attrs['cube_reg']['abs_dy']

        print("size: {}".format(ds1.dims))

        resolution = int(comboBoxResolution.currentText())
        print("loading at resolution {}".format(resolution))
        index = comboBoxResolution.currentIndex()
        print("group " + ds1.attrs["multiscale"]['datasets'][index]['path'])
        gr = ds1.attrs["multiscale"]['datasets'][index]['path']
        ds = xr.open_zarr(search_folder.text() + str(comboBoxPath.currentText())+'/mos.zarr', group=gr)


        output_resolution = float(pixel_size.text())

        with napari.gui_qt():
            if cb_C1.isChecked():
                volume_1 = (ds.sel(channel=1, type='mosaic', z=0).to_array().data * bscale + bzero)
                aligned_1 = self.Align(volume_1, resolution, output_resolution)
                viewer.add_image([aligned_1], name='C1', scale=(15,output_resolution,output_resolution), blending='additive', colormap='red')
            
            if cb_C2.isChecked():
                volume_2 = (ds.sel(channel=2, type='mosaic', z=0).to_array().data * bscale + bzero)
                aligned_2 = self.Align(volume_2, resolution, output_resolution)
                viewer.add_image([aligned_2], name='C2', scale=(15,output_resolution,output_resolution), blending='additive', colormap='green')

            if cb_C3.isChecked():
                volume_3 = (ds.sel(channel=3, type='mosaic', z=0).to_array().data * bscale + bzero)
                aligned_3 = self.Align(volume_3, resolution, output_resolution)
                viewer.add_image([aligned_3], name='C3', scale=(15,output_resolution,output_resolution), blending='additive', colormap='blue')

            if cb_C4.isChecked():
                volume_4 = (ds.sel(channel=4, type='mosaic', z=0).to_array().data * bscale + bzero)
                aligned_4 = self.Align(volume_4, resolution, output_resolution)
                viewer.add_image([aligned_4], name='C4', scale=(15,output_resolution,output_resolution), blending='additive', colormap='bop purple')



        
    def Align(self, volume, resolution, output_resolution):

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
            print('{}/{}, translate_x: {}, translate_y: {}'.format(z,z_size,alignX,alignY))

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

    def SelectFolder(self):
        global search_folder
        global comboBoxPath
        file = str(QtWidgets.QFileDialog.getExistingDirectory())
        file = file + '/'
        #text_files = glob.glob(path + "/**/**/*.zarr", recursive = False)
        #print(text_files)
        #return

        search_folder.setText(file)
        comboBoxPath.clear()
        for f in os.scandir(search_folder.text()):
            if f.is_dir():
                if os.path.exists(f.path + "/mos.zarr"):
                    s = f.path
                    s = s.replace(search_folder.text(), '')
                    print(s)
                    comboBoxPath.addItem(s)


    def on_combobox_changed(self):
        if os.path.exists(search_folder.text() + str(comboBoxPath.currentText())+'/mos.zarr/.zmetadata'):
            print("metadata is available")
            ds = xr.open_zarr(search_folder.text() + str(comboBoxPath.currentText())+'/mos.zarr', consolidated = True)
            print(ds.attrs["multiscale"]['datasets'])
            length = len(ds.attrs["multiscale"]['datasets'])
            comboBoxResolution.clear()
            for n in range(0, length):
                comboBoxResolution.addItem(str(ds.attrs["multiscale"]['datasets'][n]['level']))
        else:
            comboBoxResolution.addItem(str(1))
            comboBoxResolution.addItem(str(2))
            comboBoxResolution.addItem(str(4))
            comboBoxResolution.addItem(str(8))
            comboBoxResolution.addItem(str(16))
            comboBoxResolution.addItem(str(32))
            comboBoxResolution.setCurrentText(str(16))


    def main(self, argv):
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
            bSelectFolder.clicked.connect(self.SelectFolder)
            hbox.addWidget(bSelectFolder)
            vbox.addLayout(hbox)

            hbox = QtWidgets.QHBoxLayout()
            global comboBoxPath
            comboBoxPath = QtWidgets.QComboBox()    
            
            comboBoxPath.clear()
            if os.path.exists(search_folder.text()):
                for f in os.scandir(search_folder.text()):
                    if f.is_dir():
                        if os.path.exists(f.path + "/mos.zarr"):
                            s = f.path
                            s = s.replace(search_folder.text(), '')
                            print(s)
                            comboBoxPath.addItem(s)

            comboBoxPath.setMaximumWidth(300)
            comboBoxPath.currentIndexChanged.connect(self.on_combobox_changed)
            hbox.addWidget(comboBoxPath)
            vbox.addLayout(hbox)

            
            hbox = QtWidgets.QHBoxLayout()
            hbox.addWidget(QtWidgets.QLabel("Load level:"))
            comboBoxResolution = QtWidgets.QComboBox() 
            comboBoxResolution.setMaximumWidth(100)
            self.on_combobox_changed()
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
            bLoad.clicked.connect(self.Load)
            hbox = QtWidgets.QHBoxLayout()
            hbox.addWidget(bLoad)
            hbox.addStretch(1)
            vbox.addLayout(hbox)
            
            vbox.addStretch(1)
            widget.setLayout(vbox)
            viewer.window.add_dock_widget(widget, area="right" )

            inputfile = ''
            try:
                opts, args = getopt.getopt(argv,"hi:c:",["ifile=","channel="])
            except getopt.GetoptError:
                #print 'napari-stpt -i <inputfile> -c <channel>'
                sys.exit(2)
            for opt, arg in opts:
                if opt == '-h':
                    #print 'napari-stpt'
                    sys.exit()
                elif opt in ("-i", "--ifile"):
                    inputfile = arg
                elif opt in ("-c", "--channel"):
                    channel = channel
            #print 'Input file is "', inputfile
            #print 'channel is "', channel

if __name__ == "__main__":
    NapariSTPT().main(sys.argv[1:])
