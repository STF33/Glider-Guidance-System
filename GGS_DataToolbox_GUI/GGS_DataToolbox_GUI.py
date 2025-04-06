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
from PyQt6.QtWidgets import (QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QFrame, QCheckBox, QLineEdit, QTextEdit, QTabWidget, QFormLayout, QGroupBox, QMessageBox, QMenuBar, QStatusBar)
from PyQt6.QtCore import Qt

from GGS_DataToolbox_Main import GGS_DataToolbox_Main
from X_DataToolbox_Data import *
from X_DataToolbox_Products import *
from X_DataToolbox_Advanced import *

# =========================

config_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), "config", "config.json")

class FileDropBox(QFrame):

    ''' 
    Display a drag-and-drop file box.
    
    Arguments:
    - title (str): Tooltip title for the file box.
    - placeholder (str): Text to display initially.
    - min_size (tuple): Minimum width and height (int, int).
    - parent (QWidget, optional): Parent widget.
      
    Returns:
    - None
    '''

    def __init__(self, title, placeholder, min_size=(250, 150), parent=None):
        super().__init__(parent)
        self.setFrameStyle(QFrame.Shape.Box)
        self.setStyleSheet("border: 2px dashed #aaa; border-radius: 10px; background-color: #f0f0f0;")
        self.setMinimumSize(*min_size)
        self.setAcceptDrops(True)
        self.file_list = []
        self.label = QLabel(placeholder, self)
        self.label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.label.setStyleSheet("font-size: 16px; color: #333;")
        layout = QVBoxLayout()
        layout.addWidget(self.label)
        self.setLayout(layout)
        self.setToolTip(title)

    def dragEnterEvent(self, event):

        ''' 
        Handle drag enter events.
        
        Arguments:
        - event (QDragEnterEvent): The drag event.
        
        Returns:
        - None
        '''

        if event.mimeData().hasUrls():
            event.acceptProposedAction()

    def dropEvent(self, event):

        ''' 
        Handle drop events and update file list.
        
        Arguments:
        - event (QDropEvent): The drop event.
        
        Returns:
        - None
        '''

        new_files = [url.toLocalFile() for url in event.mimeData().urls()]
        self.file_list.extend(new_files)
        self.label.setText("\n".join(self.file_list))

    def copy_files(self, dest_folder):

        ''' 
        Copy files in the file list to the destination folder.
        
        Arguments:
        - dest_folder (str): Destination folder path.
        
        Returns:
        - None
        '''

        os.makedirs(dest_folder, exist_ok=True)
        for file in self.file_list:
            shutil.copy(file, dest_folder)
        print(f"Files copied to {dest_folder}")

class MainGUI(QMainWindow):

    ''' 
    Main GUI window with a tabbed interface.
    The "Home" tab allows configuration and file drop input.
    Additional product tabs (Plot, Energy Evaluation) are added dynamically.
    
    Arguments:
    - None
      
    Returns:
    - None
    '''

    def __init__(self):
        super().__init__()
        self.setWindowTitle("Glider Guidance System: Data Toolbox")
        self.resize(900, 700)
        self.script_directory = os.path.dirname(os.path.abspath(__file__))
        self.statusBar = QStatusBar()
        self.setStatusBar(self.statusBar)
        self.product_tabs = {}
        self.initUI()
        self.load_configuration()

    def initUI(self):

        ''' 
        Initialize the main GUI user interface.
        
        Arguments:
        - None
        
        Returns:
        - None
        '''

        menu_bar = self.menuBar()
        file_menu = menu_bar.addMenu("File")
        exit_action = file_menu.addAction("Exit")
        exit_action.triggered.connect(self.close)
        help_menu = menu_bar.addMenu("Help")
        about_action = help_menu.addAction("About")
        about_action.triggered.connect(self.show_about)
        
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        main_layout = QVBoxLayout(central_widget)
        self.tabs = QTabWidget()
        main_layout.addWidget(self.tabs)
        self.home_tab = QWidget()
        self.init_home_tab()
        self.tabs.addTab(self.home_tab, "Home")

    def init_home_tab(self):

        ''' 
        Initialize the Home tab with configuration, file inputs, and processing options.
        
        Arguments:
        - None
        
        Returns:
        - None
        '''

        layout = QVBoxLayout(self.home_tab)
        glider_group = QGroupBox("Glider Information")
        glider_layout = QFormLayout()
        self.glider_unit = QLineEdit()
        self.glider_version = QLineEdit()
        self.glider_type = QLineEdit()
        glider_layout.addRow("Unit:", self.glider_unit)
        glider_layout.addRow("Version:", self.glider_version)
        glider_layout.addRow("Type:", self.glider_type)
        glider_group.setLayout(glider_layout)
        layout.addWidget(glider_group)
        
        sensors_group = QGroupBox("Sensors")
        sensor_layout = QVBoxLayout()
        self.sensor_list = QTextEdit()
        self.sensor_list.setToolTip("Enter sensor list as a JSON array, e.g. [\"sensor1\", \"sensor2\"]")
        sensor_layout.addWidget(QLabel("Sensor List:"))
        sensor_layout.addWidget(self.sensor_list)
        sensors_group.setLayout(sensor_layout)
        layout.addWidget(sensors_group)
        
        tasks_group = QGroupBox("Processing Options")
        tasks_layout = QVBoxLayout()
        self.config_map = {
            "Decompression": ("DECOMPRESSION", "run_decompression"),
            "Ascii Conversion": ("CONVERSION", "run_conversion"),
            "Dataframe": ("DATA", "run_dataframe"),
            "Data Filter": ("DATA", "run_data_filter"),
            "Logfile Search": ("DATA", "run_logfile_search"),
            "Data Sorter": ("DATA", "run_data_sorter"),
            "Plot": ("PRODUCTS", "run_plot"),
            "Excel": ("PRODUCTS", "run_excel"),
            "Energy Evaluation": ("PRODUCTS", "run_energy_evaluation"),
            "Cleanup": ("ADVANCED", "run_data_cleanup")
        }
        self.options = {}
        for option in self.config_map.keys():
            cb = QCheckBox(f"Enable: {option}")
            self.options[option] = cb
            tasks_layout.addWidget(cb)
        tasks_group.setLayout(tasks_layout)
        layout.addWidget(tasks_group)
        
        files_group = QGroupBox("File Input")
        files_layout = QHBoxLayout()
        self.dataFileBox = FileDropBox("Data Files", "Drag and Drop Data Files Here", (250, 200))
        self.cacheFileBox = FileDropBox("Cache Files", "Drag and Drop Cache Files Here", (200, 150))
        self.logFileBox = FileDropBox("Log Files", "Drag and Drop Log Files Here", (200, 150))
        files_layout.addWidget(self.dataFileBox)
        files_layout.addWidget(self.cacheFileBox)
        files_layout.addWidget(self.logFileBox)
        files_group.setLayout(files_layout)
        layout.addWidget(files_group)
        
        self.runButton = QPushButton("Run")
        self.runButton.clicked.connect(self.run_process)
        layout.addWidget(self.runButton, alignment=Qt.AlignmentFlag.AlignCenter)

    def show_about(self):

        ''' 
        Show an About message box.
        
        Arguments:
        - None
        
        Returns:
        - None
        '''

        QMessageBox.information(self, "About", "Glider Guidance System Data Toolbox\nVersion 2.0")

    def load_configuration(self):

        ''' 
        Load configuration from JSON file and update GUI fields.
        
        Arguments:
        - None
        
        Returns:
        - None
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
            self.statusBar.showMessage("Configuration loaded", 3000)

    def save_configuration(self):

        ''' 
        Save current configuration to JSON file.
        
        Arguments:
        - None
        
        Returns:
        - None
        '''

        sensor_text = self.sensor_list.toPlainText().strip()
        if sensor_text == "":
            sensor_list_value = []
        else:
            try:
                sensor_list_value = json.loads(sensor_text)
            except Exception as e:
                print(f"Error decoding sensor_list: {e}")
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
        self.statusBar.showMessage("Configuration saved", 3000)

    def run_process(self):

        ''' 
        Save configuration, copy files to working directories, run main process, and add product tabs.
        
        Arguments:
        - None
        
        Returns:
        - None
        '''

        self.save_configuration()
        self.dataFileBox.copy_files(os.path.join(self.script_directory, "DBD_Files"))
        self.cacheFileBox.copy_files(os.path.join(self.script_directory, "cache"))
        self.logFileBox.copy_files(os.path.join(self.script_directory, "DBD_Files", "Logfiles"))
        self.statusBar.showMessage("Files copied. Running main process...", 3000)
        dataframe, root_directory = GGS_DataToolbox_Main(config_name="config")
        if self.options["Plot"].isChecked() and dataframe is not None:
            self.add_plot_tab(dataframe)
        if self.options["Energy Evaluation"].isChecked() and dataframe is not None:
            self.add_energy_tab(dataframe, root_directory)

    def add_plot_tab(self, dataframe):

        ''' 
        Add a Plot tab if not already present.
        
        Arguments:
        - dataframe (pd.DataFrame): Merged sensor data.
        
        Returns:
        - None
        '''

        if "Plot" not in self.product_tabs:
            plot_widget = PlotGUI(dataframe)
            self.product_tabs["Plot"] = plot_widget
            self.tabs.addTab(plot_widget, "Plot")
            self.statusBar.showMessage("Plot tab added", 3000)
        else:
            self.statusBar.showMessage("Plot tab already exists", 3000)

    def add_energy_tab(self, dataframe, root_directory):

        ''' 
        Add an Energy Evaluation tab if not already present.
        
        Arguments:
        - dataframe (pd.DataFrame): Merged sensor data.
        - root_directory (str): The output directory from the configuration.
        
        Returns:
        - None
        '''

        if "Energy Evaluation" not in self.product_tabs:
            glider_info = (self.glider_unit.text(), self.glider_version.text(), self.glider_type.text())
            energy_widget = run_energy_evaluation(root_directory, glider_info, dataframe)
            self.product_tabs["Energy Evaluation"] = energy_widget
            self.tabs.addTab(energy_widget, "Energy Eval")
            self.statusBar.showMessage("Energy Evaluation tab added", 3000)
        else:
            self.statusBar.showMessage("Energy Evaluation tab already exists", 3000)

if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = MainGUI()
    window.show()
    sys.exit(app.exec())
