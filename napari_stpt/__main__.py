"""
napari-STPT reads zarr files and displays them
"""
import sys
import napari_stpt

if __name__ == "__main__":
    napari_stpt.NapariSTPT().main(sys.argv[1:])
