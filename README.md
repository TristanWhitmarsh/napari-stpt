# napari-stpt
A napari viewer which reads stpt images as zarr files

## Installation

You can install the Napari STPT viewer from [PyPI](https://pypi.org/project/napari-stpt/):

    pip install napari-stpt

The viewer is supported on Python 3.7 and above.

Type the following to install napari movie:

    pip install git+https://github.com/guiwitz/naparimovie.git@master#egg=naparimovie


## How to run

The napari viewer can be run from command line by typing:

    napari-stpt
    
## How to use

    - Select the folder containing the folders with STPT Zarr files (windows server: N:/stpt/; linux: /data/meds1_c/storage/processed/stpt)
    - Select the image file
    - The "Load level" is the level in which the raw data is loaded. Level 1 corresponds to an approximate 0.56um pixel size, level 2: 2x0.56um, level 4: 4x0.56um etc.
    - The "Output pixel size" is the resolution to which the images are reformatted after applying the translations. It would be best to make sure the raw data pixel size is smaller than this size.
    - The "Maximum number of optical slices" can be set in case the optical slices go beyond the slice thickness of 15um.
    - Select which channels to load.
    - Press "Load image 3D".
    - To crop the volume draw a shap using the add rectangles button of napari.
    - Pressing "Crop" will just crop the colume to this region
    - Pressing "Reload image 3D in region" will reload the slices. In this case you can set a different Load level or Output pixel size.
    - Press "Save volume" to save the multi-channel volume to a tiff.

    - At any time you can us the keys defined in https://github.com/guiwitz/naparimovie to set key frames. Press "Save movie" to save a gif (save as mp4 seems to get stuck in a loop)
