    
import sys
import re

# if 'PyQt5' in sys.modules:
#     print("Using PyQt5")
#     from qtpy import QtCore, QtWidgets
#     from qtpy.QtWidgets import QComboBox, QApplication, QCompleter, QMessageBox
#     from qtpy.QtCore import QSortFilterProxyModel, Qt
#     from qtpy.QtGui import QColor, QPixmap, QIcon, QStandardItemModel, QStandardItem, QPainter, QFont
#     from qtpy.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QFrame, QScrollArea,
#                                 QComboBox, QCheckBox, QListWidget, QListWidgetItem, QMessageBox, QFileDialog, QMainWindow)
# else:
print("Using PySide2")
from PySide2 import QtCore, QtWidgets, QtGui
from PySide2.QtCore import Qt, QSortFilterProxyModel
from PySide2.QtGui import QColor, QPixmap, QIcon, QStandardItemModel, QStandardItem, QPainter, QFont
from PySide2.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QFrame, QScrollArea,
                               QComboBox, QCheckBox, QListWidget, QListWidgetItem, QMessageBox, QFileDialog, QMainWindow, QCompleter)


class ExtendedComboBox(QComboBox):
    def __init__(self, parent=None):
        super(ExtendedComboBox, self).__init__(parent)

        self.setFocusPolicy(Qt.StrongFocus)
        self.setEditable(True)

        # add a filter model to filter matching items
        self.pFilterModel = QSortFilterProxyModel(self)
        self.pFilterModel.setFilterCaseSensitivity(Qt.CaseInsensitive)
        self.pFilterModel.setSourceModel(self.model())

        # add a completer, which uses the filter model
        self.completer = QCompleter(self.pFilterModel, self)
        # always show all (filtered) completions
        self.completer.setCompletionMode(QCompleter.UnfilteredPopupCompletion)
        self.setCompleter(self.completer)

        # connect signals
        self.lineEdit().textEdited[str].connect(self.pFilterModel.setFilterFixedString)
        self.completer.activated.connect(self.on_completer_activated)


    # on selection of an item from the completer, select the corresponding item from combobox
    def on_completer_activated(self, text):
        if text:
            index = self.findText(text)
            self.setCurrentIndex(index)
            self.activated[str].emit(self.itemText(index))


    # on model change, update the models of the filter and completer as well
    def setModel(self, model):
        super(ExtendedComboBox, self).setModel(model)
        self.pFilterModel.setSourceModel(model)
        self.completer.setModel(self.pFilterModel)


    # on model column change, update the model column of the filter and completer as well
    def setModelColumn(self, column):
        self.completer.setCompletionColumn(column)
        self.pFilterModel.setFilterKeyColumn(column)
        super(ExtendedComboBox, self).setModelColumn(column)

        
        
class CheckableItemWidget(QtWidgets.QWidget):
    def __init__(self, text, parent=None):
        super().__init__(parent)

        self.checkBox = QtWidgets.QCheckBox()
        self.label = QtWidgets.QLabel(text)

        layout = QtWidgets.QHBoxLayout(self)
        layout.addWidget(self.checkBox)
        layout.addWidget(self.label)
        layout.addStretch()

        # Adjust these margins and spacing as needed
        layout.setContentsMargins(0, 0, 0, 0)  # Reduce the margins
        layout.setSpacing(0)  # Reduce the spacing between the checkbox and label

        self.setLayout(layout)

    def isChecked(self):
        return self.checkBox.isChecked()

    

class CustomComboBox(QtWidgets.QComboBox):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._listWidget = QtWidgets.QListWidget()
        self.setModel(self._listWidget.model())
        self.setView(self._listWidget)

    def addItem(self, text):
        item = QtWidgets.QListWidgetItem(self._listWidget)
        widget = CheckableItemWidget(text)
        item.setSizeHint(widget.sizeHint())
        self._listWidget.setItemWidget(item, widget)

    def itemChecked(self, index):
        item = self._listWidget.item(index)
        widget = self._listWidget.itemWidget(item)
        return widget.isChecked()
    
    def checkItem(self, index, check=True):
        if 0 <= index < self._listWidget.count():
            item = self._listWidget.item(index)
            widget = self._listWidget.itemWidget(item)
            if widget:
                widget.checkBox.setChecked(check)
    
    def clearItems(self):
        self._listWidget.clear()

    def getCheckedItems(self):
        checkedItems = []
        for index in range(self._listWidget.count()):
            item = self._listWidget.item(index)
            check = self._listWidget.itemWidget(item)
            if check.isChecked():
                #checkedItems.append(item.text())
                checkedItems.append(index)
        return checkedItems
    
    
    
    
    


class LegendWidget(QWidget):
    def __init__(self):
        super().__init__()
        self.setAutoFillBackground(True)
        p = self.palette()
        p.setColor(self.backgroundRole(), Qt.black)  # Set the widget's background to black
        self.setPalette(p)

        self.layout = QVBoxLayout(self)  # Set the main layout for the widget

        # Create a scroll area to contain the frame
        self.scrollArea = QScrollArea(self)
        self.scrollArea.setWidgetResizable(True)  # Allow the contained widget to resize
        self.layout.addWidget(self.scrollArea)  # Add the scroll area to the widget's main layout

        # Frame setup
        self.frame = QFrame()
        self.frame.setAutoFillBackground(False)
        self.frameLayout = QVBoxLayout(self.frame)  # Create a layout for the frame
        self.frame.setLayout(self.frameLayout)  # Set the layout on the frame

        self.scrollArea.setWidget(self.frame)  # Set the frame as the scroll area's widget

        self.channel_names = None
        self.colors = None

        # Ensure there's a stretch at the end of the frameLayout to push items to the top
        self.frameLayout.addStretch(1)

    def clear_layout(self, layout):
        """Clears a layout."""
        while layout.count():
            item = layout.takeAt(0)
            if item.layout() is not None:
                # Clear nested layout first
                self.clear_layout(item.layout())
                item.layout().deleteLater()
            elif item.widget() is not None:
                item.widget().deleteLater()

    def populate_legend(self, channel_names, colors):
        # First, clear the existing content in the layout
        self.clear_layout(self.frameLayout)
        
        self.channel_names = channel_names
        self.colors = colors

        # Removing the stretch before repopulating
        while self.frameLayout.count() > 0:
            self.frameLayout.takeAt(0)

        for name, color in zip(channel_names, colors):
            channel_layout = QHBoxLayout()
            color_label = QLabel()
            pixmap = QPixmap(16, 16)
            pixmap.fill(QColor(*color))
            color_label.setPixmap(pixmap)

            name_str = str(name) if isinstance(name, int) else name
            
            name_str2 = name_str
            # name_str2 = re.sub(r'\d+\s*', '', name_str, count=1)
            # if name_str2 == "" or name_str2 == " ":
            #     name_str2 = name_str
            
            text_label = QLabel(name_str2)
            text_label.setStyleSheet("QLabel { color : white; }")

            channel_layout.addWidget(color_label)
            channel_layout.addWidget(text_label)
            channel_layout.addStretch(1)

            self.frameLayout.addLayout(channel_layout)

        # Add a final stretch to ensure all items are aligned at the top
        self.frameLayout.addStretch(1)

        self.add_save_button()  # Consider whether you need to re-add the save button every time

    def add_save_button(self):
        # Check if the save button already exists
        if hasattr(self, 'save_button'):
            return  # The button already exists, so don't add it again

        self.save_button = QPushButton("Save to Image")
        self.save_button.clicked.connect(self.save_legend_to_image)
        self.layout.addWidget(self.save_button)  # Add the save button to the widget's main layout



    def save_legend_to_image(self):
        filename, _ = QFileDialog.getSaveFileName(self, "Save Image", "", "BMP Files (*.bmp);;All Files (*)")
        if filename:
            self.render_legend_to_image(filename)

    def render_legend_to_image(self, filename):
        width = 200
        height_per_item = 20
        canvas_height = len(self.channel_names) * height_per_item

        pixmap = QPixmap(width, canvas_height)
        pixmap.fill(Qt.black)

        painter = QPainter(pixmap)
        painter.setFont(QFont("Arial", 10))
        text_offset_x = 30
        box_size = 15

        for index, (name, color) in enumerate(zip(self.channel_names, self.colors)):
            painter.setBrush(QColor(*color))
            painter.setPen(Qt.NoPen)

            y_pos = index * height_per_item
            painter.drawRect(5, y_pos + 2, box_size, box_size)
            painter.setPen(Qt.white)
            
            name_str = str(name) if isinstance(name, int) else name
            
            name_str2 = re.sub(r'\d+\s*', '', name_str, count=1)
            if name_str2 == "" or name_str2 == " ":
                name_str2 = name_str
                
            painter.drawText(text_offset_x, y_pos + box_size, name_str2)

        painter.end()
        pixmap.save(filename, "BMP")
        print(f"Legend saved as {filename}")



        