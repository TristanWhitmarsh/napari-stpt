# custom_qt_dims_slider.py

from napari._qt.widgets.qt_dims_slider import QtDimSliderWidget
from qtpy.QtWidgets import QLineEdit
from qtpy.QtCore import Qt
from qtpy.QtWidgets import QWidget

class CustomQtDimSliderWidget(QtDimSliderWidget):
    def __init__(self, dims, axis, parent=None):
        super().__init__(dims, axis, parent)
        
    # def _update_slice_labels(self):
    #     """Update slice labels to match current dimension slider position."""
    #     self.curslice_label.setText(str(self.dims.current_step[self.axis]+1))
    #     self.curslice_label.setAlignment(Qt.AlignmentFlag.AlignRight)

    def _update_slice_labels(self):
        self.curslice_label.setText(str(self.dims.current_step[self.axis] + 1))
        self.curslice_label.setAlignment(Qt.AlignmentFlag.AlignRight)
        nsteps = self.dims.nsteps[self.axis] - 1
        self.totslice_label.setText(str(nsteps + 1))
        #self.axis_label.setText("Slide #")

        for labl in self.findChildren(QWidget, 'slice_label'):
            labl.setFixedWidth(8 * 4 + 2)
        
        
    def _update_range(self):
        """Updates range for slider."""
        # super()._update_range()
        displayed_sliders = self.qt_dims._displayed_sliders

        nsteps = self.dims.nsteps[self.axis] - 1
        if nsteps == 0:
            displayed_sliders[self.axis] = False
            self.qt_dims.last_used = 0
            self.hide()
        else:
            if (
                not displayed_sliders[self.axis]
                and self.axis not in self.dims.displayed
            ):
                displayed_sliders[self.axis] = True
                self.last_used = self.axis
                self.show()
            self.slider.setMinimum(0)
            self.slider.setMaximum(nsteps)
            self.slider.setSingleStep(1)
            self.slider.setPageStep(1)
            self.slider.setValue(self.dims.current_step[self.axis])
            self.totslice_label.setText(str(nsteps+1))
            self.totslice_label.adjustSize()
            self.totslice_label.setAlignment(Qt.AlignmentFlag.AlignLeft)
            self.totslice_label.adjustSize()
            self._update_slice_labels()