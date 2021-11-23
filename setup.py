"""Setup script for realpython-napari-stpt"""

import os.path
from setuptools import setup

# The directory containing this file
HERE = os.path.abspath(os.path.dirname(__file__))

# The text of the README file
with open(os.path.join(HERE, "README.md")) as fid:
    README = fid.read()

# this grabs the requirements from requirements.txt
#REQUIREMENTS = [i.strip() for i in open("requirements.txt").readlines()]

# This call to setup() does all the work
setup(
    name="napari-stpt",
    version="0.0.8.23",
    description="napari viewer which can read stpt images as zarr files",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/TristanWhitmarsh/napari-stpt",
    author="Tristan Whitmarsh",
    author_email="tw401@cam.ac.uk",
    license="GNU",
    classifiers=[
        "License :: OSI Approved :: GNU Affero General Public License v3",
        "Programming Language :: Python :: 3",
    ],
    packages=["napari_stpt"],
    include_package_data=True,
    install_requires=[
        'napari',
        'numpy',
        'xarray',
        'qtpy',
        'pyqt5',
        'SimpleITK',
        'opencv-python',
        'zarr',
        'stardist',
        'tensorflow'
    ],
    entry_points={"console_scripts": ["napari-stpt=napari_stpt.__main__:main"]},
)
