"""
napari-STPT reads zarr files and displays them
"""
# import sys

try:
    from napari_stpt.napari_stpt import NapariSTPT
except:
    from napari_stpt import NapariSTPT

def main():
    napari = NapariSTPT()
    napari.main()

if __name__ == "__main__":
    main()