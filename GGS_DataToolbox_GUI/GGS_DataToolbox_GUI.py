"""
Author: Salvatore Fricano
Repository: https://github.com/STF33/Glider-Guidance-System
"""

# =========================
# IMPORTS
# =========================

import sys
import os
import shutil
import json
from PyQt6.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, 
                             QPushButton, QFrame, QCheckBox, QLineEdit, QTextEdit)
from PyQt6.QtCore import Qt

from GGS_DataToolbox_Main import GGS_DataToolbox_Main

# =========================

config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "config", "config.json")

### CLASS:
class GUI_FileBox(QFrame):
    
    '''
    A drag-and-drop file box for selecting data files.
    
    Args:
      parent (QWidget, optional): The parent widget.
      
    Returns:
      None
    '''

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setFrameStyle(QFrame.Shape.Box)
        self.setStyleSheet("border: 2px dashed #aaa; border-radius: 10px; background-color: #f0f0f0;")
        self.setMinimumSize(250, 200)
        self.setAcceptDrops(True)
        self.file_list = []
        self.label = QLabel("Drag and Drop Data Files Here", self)
        self.label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.label.setStyleSheet("font-size: 18px; color: black; font-weight: bold;")
        layout = QVBoxLayout()
        layout.addWidget(self.label)
        self.setLayout(layout)

    def dragEnterEvent(self, event):
        
        '''
        Handle the drag enter event to accept file URLs.
        
        Args:
          event (QDragEnterEvent): The drag enter event.
          
        Returns:
          None
        '''

        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event):
        
        '''
        Handle the drop event by collecting the dropped file paths and updating the label.
        
        Args:
          event (QDropEvent): The drop event.
          
        Returns:
          None
        '''

        new_files = [url.toLocalFile() for url in event.mimeData().urls()]
        self.file_list.extend(new_files)
        self.label.setText("\n".join(self.file_list))

    def GUI_FileBox_Copy(self, script_directory):
        
        '''
        Copy the files from the file list to the DBD_Files directory.
        
        Args:
          script_directory (str): The base directory where the DBD_Files folder exists.
          
        Returns:
          None
        '''

        dbd_root = os.path.join(script_directory, "DBD_Files")
        os.makedirs(dbd_root, exist_ok=True)
        for file in self.file_list:
            shutil.copy(file, dbd_root)
        print("Files copied to DBD_Files.")


### CLASS:
class CacheFileBox(GUI_FileBox):
    
    '''
    A drag-and-drop file box for selecting cache files.
    
    Args:
      parent (QWidget, optional): The parent widget.
      
    Returns:
      None
    '''

    def __init__(self, parent=None):
        
        super().__init__(parent)
        self.setMinimumSize(200, 150)
        self.label.setText("Drag and Drop Cache Files Here")
    
    def GUI_FileBox_Copy(self, script_directory):
        
        '''
        Copy the files from the file list to the cache folder.
        
        Args:
          script_directory (str): The base directory where the cache folder will be created.
          
        Returns:
          None
        '''

        cache_dir = os.path.join(script_directory, "cache")
        os.makedirs(cache_dir, exist_ok=True)
        for file in self.file_list:
            shutil.copy(file, cache_dir)
        print("Cache files copied to cache.")


### CLASS:
class GUI_MainWindow(QWidget):
    
    '''
    Main interface for the Data Toolbox GUI.
    
    Args:
      None
      
    Returns:
      None
    '''

    def __init__(self):
        
        super().__init__()
        self.setWindowTitle("Glider Guidance System: Data Toolbox GUI")
        self.setGeometry(100, 100, 600, 600)
        self.script_directory = os.path.dirname(os.path.abspath(__file__))
        self.execution_order = []
        self.tasks = {}
        
        self.config_map = {
            "Decompression": ("DECOMPRESSION", "run_decompression"),
            "Ascii Conversion": ("CONVERSION", "run_conversion"),
            "Dataframe": ("DATA", "run_dataframe"),
            "Plot": ("PRODUCTS", "run_plot"),
            "Excel": ("PRODUCTS", "run_excel"),
            "Data Sorter": ("PRODUCTS", "run_data_sorter"),
            "Cleanup": ("ADVANCED", "run_data_cleanup")
        }
        
        main_layout = QVBoxLayout()
        
        self.title_label = QLabel("Glider Guidance System: Data Toolbox GUI")
        self.title_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.title_label.setStyleSheet("font-size: 18px; font-weight: bold; padding: 2px; margin-bottom: 2px;")
        main_layout.addWidget(self.title_label, 0, Qt.AlignmentFlag.AlignTop)
        
        content_layout = QHBoxLayout()
        content_layout.setStretch(0, 1)
        content_layout.setStretch(1, 1)
        
        left_layout = QVBoxLayout()
        self.dataFileBox = GUI_FileBox(self)
        self.cacheFileBox = CacheFileBox(self)
        left_layout.addWidget(self.dataFileBox)
        left_layout.addWidget(self.cacheFileBox)
        content_layout.addLayout(left_layout, 1)
        
        right_section = QVBoxLayout()
        right_section.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.glider_unit = QLineEdit()
        self.glider_version = QLineEdit()
        self.glider_type = QLineEdit()
        self.sensor_list = QTextEdit()
        self.options = {}
        for option in self.config_map.keys():
            self.options[option] = QCheckBox(f"Enable: {option}")
        
        self.GUI_Main_CreateSection(right_section, "GLIDER", ["Unit:", self.glider_unit, "Version:", self.glider_version, "Type:", self.glider_type])
        self.GUI_Main_CreateSection(right_section, "SENSORS", ["List:", self.sensor_list])
        self.GUI_Main_CreateSection(right_section, "DECOMPRESSION", [self.options["Decompression"]])
        self.GUI_Main_CreateSection(right_section, "CONVERSION", [self.options["Ascii Conversion"]])
        self.GUI_Main_CreateSection(right_section, "DATA", [self.options["Dataframe"]])
        self.GUI_Main_CreateSection(right_section, "PRODUCTS", [self.options["Plot"], self.options["Excel"], self.options["Data Sorter"]])
        self.GUI_Main_CreateSection(right_section, "ADVANCED", [self.options["Cleanup"]])
        
        self.runButton = QPushButton("Run")
        self.runButton.clicked.connect(self.GUI_Main_RunFunction)
        right_section.addWidget(self.runButton, alignment=Qt.AlignmentFlag.AlignCenter)
        content_layout.addLayout(right_section, 1)
        main_layout.addLayout(content_layout, 1)
        self.setLayout(main_layout)
        
        self.GUI_Main_ConfigLoad()

    def GUI_Main_CreateSection(self, layout, title, widgets):
        
        '''
        Create a labeled section in the GUI and add the provided widgets.

        Args:
          layout (QLayout): The parent layout to add the section to.
          title (str): The title of the section.
          widgets (list): A list of widgets or strings to add to the section.
          
        Returns:
          None
        '''

        section_label = QLabel(title)
        section_label.setStyleSheet("font-size: 16px; font-weight: bold;")
        layout.addWidget(section_label, alignment=Qt.AlignmentFlag.AlignCenter)
        for widget in widgets:
            if isinstance(widget, str):
                label = QLabel(widget)
                layout.addWidget(label)
            else:
                layout.addWidget(widget)

    def GUI_Main_ConfigLoad(self):
        
        '''
        Load the configuration from the JSON file and update the GUI fields.

        Args:
          None
          
        Returns:
          None
        '''

        if os.path.exists(config_path):
            with open(config_path, "r") as f:
                config = json.load(f)
                self.glider_unit.setText(config["GLIDER"].get("glider_unit", ""))
                self.glider_version.setText(config["GLIDER"].get("glider_version", ""))
                self.glider_type.setText(config["GLIDER"].get("glider_type", ""))
                self.sensor_list.setText(json.dumps(config["SENSORS"].get("sensor_list", [])))
                
                for option, (section, key) in self.config_map.items():
                    if section in config and key in config[section]:
                        self.options[option].setChecked(config[section][key])

    def GUI_Main_ConfigSave(self):
        
        '''
        Save the current GUI configuration to the JSON config file.

        Args:
          None
          
        Returns:
          None
        '''

        sensor_text = self.sensor_list.toPlainText().strip()
        if sensor_text == "":
            sensor_list_value = []
        else:
            try:
                sensor_list_value = json.loads(sensor_text)
            except Exception as e:
                print(f"Error decoding sensor_list, defaulting to empty list: {e}")
                sensor_list_value = []
                
        config = {
            "GLIDER": {
                "glider_unit": self.glider_unit.text(),
                "glider_version": self.glider_version.text(),
                "glider_type": self.glider_type.text()
            },
            "SENSORS": {
                "sensor_list": sensor_list_value
            }
        }
        for option, (section, key) in self.config_map.items():
            if section not in config:
                config[section] = {}
            config[section][key] = self.options[option].isChecked()
        os.makedirs(os.path.dirname(config_path), exist_ok=True)
        with open(config_path, "w") as f:
            json.dump(config, f, indent=4)

    def GUI_Main_RunFunction(self):
        
        '''
        Save the configuration, copy the dragged files to their designated folders, and run the main Data Toolbox function.

        Args:
          None
          
        Returns:
          None
        '''

        self.GUI_Main_ConfigSave()
        self.dataFileBox.GUI_FileBox_Copy(self.script_directory)
        self.cacheFileBox.GUI_FileBox_Copy(self.script_directory)
        GGS_DataToolbox_Main(config_name="config")

### MAIN:
if __name__ == "__main__":
    
    '''
    Create the application, instantiate the GUI_MainWindow, and start the event loop.

    Args:
      None
      
    Returns:
      None
    '''

    app = QApplication(sys.argv)
    window = GUI_MainWindow()
    window.show()
    sys.exit(app.exec())
