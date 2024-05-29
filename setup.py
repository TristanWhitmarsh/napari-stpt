"""Setup script for realpython-napari-stpt"""

import os.path
from setuptools import setup, find_packages

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
    version="0.1.2.6",
    description="napari viewer which can read stpt, axio and imc images as zarr files",
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
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        'opencv-python>=4.5.1.48',
        'napari==0.4.19',
        'numpy==1.24.3',
        'xarray==2023.4.2',
        'zarr==2.14.2',
        'PySide2==5.13.2',
        'SimpleITK==2.2.1',
        'napari-animation==0.0.7',
        'tifffile==2023.4.12'
    ],
    entry_points={"console_scripts": ["napari-stpt=napari_stpt.__main__:main"]},
)
