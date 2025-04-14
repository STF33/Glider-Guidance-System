"""
Author: Salvatore Fricano
Repository: https://github.com/STF33/Glider-Guidance-System
"""

# =========================
# IMPORTS
# =========================

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from PyQt6.QtWidgets import (QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QComboBox, QListWidget, QListWidgetItem, QMessageBox, QSplitter, QScrollArea, QLineEdit, QDialog, QGridLayout, QDial)
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QColor
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
import os
from openpyxl import load_workbook
import shutil
import re
from datetime import datetime
import shutil
import pandas as pd
import json

# =========================

class plot(QWidget):

    ''' 
    Create an interactive plot GUI using PyQt and Matplotlib.
    
    Arguments:
    - dataframe (pd.DataFrame): The merged sensor data containing a "time" column and sensor values.
      
    Returns:
    - None
    '''

    def __init__(self, dataframe):
        super().__init__()
        self.dataframe = dataframe
        self.selected_variables = []
        self.variable_colors = {}
        self.available_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        self.plot_init()
    
    def plot_init(self):

        ''' 
        Initialize the interactive plot GUI layout, including the plot area and control panel.
        
        Arguments:
        - None
        
        Returns:
        - None
        '''

        if "time" not in self.dataframe.columns:
            QMessageBox.warning(self, "Missing Time Column", "The merged dataframe must have a 'time' column.")
            return

        self.x_variable = "time"
        
        main_layout = QVBoxLayout()
        splitter = QSplitter(Qt.Orientation.Vertical)
        
        plot_widget = QWidget()
        plot_layout = QVBoxLayout(plot_widget)
        self.figure, self.ax = plt.subplots()
        self.canvas = FigureCanvas(self.figure)
        plot_layout.addWidget(self.canvas)
        self.toolbar = NavigationToolbar(self.canvas, self)
        plot_layout.addWidget(self.toolbar)
        splitter.addWidget(plot_widget)
        splitter.setStretchFactor(0, 4)
        
        control_widget = QWidget()
        control_layout = QVBoxLayout(control_widget)
        top_control_layout = QHBoxLayout()
        self.variable_selector = QComboBox()
        variables = sorted([col for col in self.dataframe.columns if col != "time"])
        self.variable_selector.addItems(variables)
        top_control_layout.addWidget(QLabel("Select Variable:"))
        top_control_layout.addWidget(self.variable_selector)
        self.add_remove_button = QPushButton("Add/Remove")
        self.add_remove_button.clicked.connect(self.plot_toggle_variable)
        top_control_layout.addWidget(self.add_remove_button)
        control_layout.addLayout(top_control_layout)
        self.variable_list = QListWidget()
        control_layout.addWidget(self.variable_list)
        self.plot_button = QPushButton("Update Plot")
        self.plot_button.clicked.connect(self.plot_update_plot)
        control_layout.addWidget(self.plot_button)
        splitter.addWidget(control_widget)
        splitter.setStretchFactor(1, 1)
        
        main_layout.addWidget(splitter)
        self.setLayout(main_layout)
        self.setWindowTitle("Interactive Plot GUI")
    
    def plot_toggle_variable(self):

        ''' 
        Toggle the inclusion of a variable for plotting when the user clicks the "Add/Remove" button.
        
        Arguments:
        - None
        
        Returns:
        - None
        '''

        var = self.variable_selector.currentText()
        if var in self.selected_variables:
            self.selected_variables.remove(var)
            del self.variable_colors[var]
            for i in range(self.variable_list.count()):
                if self.variable_list.item(i).text() == var:
                    self.variable_list.takeItem(i)
                    break
        else:
            self.selected_variables.append(var)
            used_colors = [self.variable_colors[v] for v in self.selected_variables if v in self.variable_colors]
            available_colors = [color for color in self.available_colors if color not in used_colors]
            color = available_colors[0] if available_colors else self.available_colors[len(self.selected_variables) % len(self.available_colors)]
            self.variable_colors[var] = color
            item = QListWidgetItem(var)
            item.setForeground(QColor(color))
            self.variable_list.addItem(item)
    
    def plot_update_plot(self):

        ''' 
        Update the plot based on the selected variables, grouping y-axes by unit and adjusting layout.
        
        Arguments:
        - None
        
        Returns:
        - None
        '''

        self.figure.clf()
        
        left_margin = 0.1
        right_margin = 0.95
        top_margin = 0.95
        bottom_margin = 0.1
        
        host = self.figure.add_subplot(111)
        x_data = self.dataframe["time"]
        host.set_xlabel("time")
        
        units_dict = self.dataframe.attrs.get("units", {})
        unit_axes = {}
        axes_list = [host]
        variable_handles = []
        
        for idx, var in enumerate(self.selected_variables):
            unit = units_dict.get(var, "unknown")
            if unit in unit_axes:
                ax = unit_axes[unit]
            else:
                if not unit_axes:
                    ax = host
                else:
                    ax = host.twinx()
                    pos = 1 + 0.1 * len(unit_axes)
                    ax.spines["right"].set_position(("axes", pos))
                unit_axes[unit] = ax
                axes_list.append(ax)
            
            y_data = pd.to_numeric(self.dataframe[var], errors='coerce').dropna()
            x_data_cleaned = x_data[y_data.index]
            
            if var == 'm_depth':
                ax.invert_yaxis()
            
            color = self.variable_colors[var]
            handle = ax.scatter(x_data_cleaned, y_data, color=color, label=var, marker='o')
            ax.plot(x_data_cleaned, y_data, color=color)
            variable_handles.append(handle)
            
            if ax.get_ylabel() == "":
                ax.set_ylabel(f"{unit}")
        
        host.xaxis_date()
        locator = mdates.AutoDateLocator()
        formatter = mdates.ConciseDateFormatter(locator)
        host.xaxis.set_major_locator(locator)
        host.xaxis.set_major_formatter(formatter)
        
        legend_handles = [handle for handle in variable_handles]
        legend_labels = [handle.get_label() for handle in variable_handles]
        host.legend(
            handles=legend_handles, labels=legend_labels, loc="upper center",
            bbox_to_anchor=(0.5, 1.15), ncol=len(legend_handles),
            borderaxespad=0, frameon=True
        )
        
        n_extra = len(axes_list) - 1
        new_right = max(right_margin - 0.05 * n_extra, 0.5)
        self.figure.subplots_adjust(left=left_margin, right=new_right, top=top_margin, bottom=bottom_margin)
        
        self.canvas.draw()

def plot_run(dataframe):

    ''' 
    Create and return an instance of plot with the provided dataframe.
    This widget can then be embedded within a tab of the main GUI.
    
    Arguments:
    - dataframe (pd.DataFrame): The merged sensor data containing a "time" column and sensor values.
      
    Returns:
    - plot_widget (plot): The plot widget ready for embedding.
    '''

    return plot(dataframe)

# =========================
# EXCEL PROCESSING
# =========================

def excel_run(root_directory, glider_info, data_frame):

    ''' 
    Save the DataFrame to an Excel file.
    
    Arguments:
    - root_directory (str): The root directory from the configuration file.
    - glider_info (tuple): A tuple containing glider unit, version, and type information.
    - data_frame (pd.DataFrame): The DataFrame containing sensor data.
      
    Returns:
    - None
    '''
    
    print(f"\n### RUNNING: EXCEL ###\n")
    
    glider_unit, glider_version, glider_type = glider_info
    output_file = os.path.join(root_directory, f"{glider_unit}-{glider_version}-{glider_type}_DataOutput.xlsx")
    writer = pd.ExcelWriter(output_file, engine='openpyxl')
    data_frame.to_excel(writer, index=False)
    writer.close()
    
    workbook = load_workbook(output_file)
    worksheet = workbook.active
    for col in worksheet.columns:
        worksheet.column_dimensions[col[0].column_letter].width = 20
    workbook.save(output_file)
    
    print(f"All data saved to {output_file}")

# =========================
# DATA SORTING
# =========================

def datasorter_organize_data_files(input_data_folder, output_folder):

    ''' 
    Organize data files into mission-specific folders.
    
    Arguments:
    - input_data_folder (str): The directory containing the input data files.
    - output_folder (str): The directory where the organized data files will be saved.
      
    Returns:
    - None
    '''

    file_groups = []
    file_groups = {}

    sys_log_path = os.path.join(input_data_folder, 'sys.log')
    sys_log_exists = os.path.exists(sys_log_path)

    for filename in os.listdir(input_data_folder):
        if filename.endswith('.dbd'):
            with open(os.path.join(input_data_folder, filename), 'r', encoding='utf-8', errors='ignore') as file:
                lines = file.readlines()
                mission_name_line = lines[7]
                mission_name = mission_name_line.split('mission_name:')[1].strip()
            base_name = os.path.splitext(filename)[0]
            if mission_name not in file_groups:
                file_groups[mission_name] = []
            file_groups[mission_name].append(base_name)
        elif filename.endswith('.mbd'):
            with open(os.path.join(input_data_folder, filename), 'r', encoding='utf-8', errors='ignore') as file:
                lines = file.readlines()
                mission_name_line = lines[7]
                mission_name = mission_name_line.split('mission_name:')[1].strip()
            base_name = os.path.splitext(filename)[0]
            if mission_name not in file_groups:
                file_groups[mission_name] = []
            file_groups[mission_name].append(base_name)
        elif filename.endswith('.sbd'):
            with open(os.path.join(input_data_folder, filename), 'r', encoding='utf-8', errors='ignore') as file:
                lines = file.readlines()
                mission_name_line = lines[7]
                mission_name = mission_name_line.split('mission_name:')[1].strip()
            base_name = os.path.splitext(filename)[0]
            if mission_name not in file_groups:
                file_groups[mission_name] = []
            file_groups[mission_name].append(base_name)
        elif filename.endswith('.ebd'):
            with open(os.path.join(input_data_folder, filename), 'r', encoding='utf-8', errors='ignore') as file:
                lines = file.readlines()
                mission_name_line = lines[7]
                mission_name = mission_name_line.split('mission_name:')[1].strip()
            base_name = os.path.splitext(filename)[0]
            if mission_name not in file_groups:
                file_groups[mission_name] = []
            file_groups[mission_name].append(base_name)
        elif filename.endswith('.nbd'):
            with open(os.path.join(input_data_folder, filename), 'r', encoding='utf-8', errors='ignore') as file:
                lines = file.readlines()
                mission_name_line = lines[7]
                mission_name = mission_name_line.split('mission_name:')[1].strip()
            base_name = os.path.splitext(filename)[0]
            if mission_name not in file_groups:
                file_groups[mission_name] = []
            file_groups[mission_name].append(base_name)
        elif filename.endswith('.tbd'):
            with open(os.path.join(input_data_folder, filename), 'r', encoding='utf-8', errors='ignore') as file:
                lines = file.readlines()
                mission_name_line = lines[7]
                mission_name = mission_name_line.split('mission_name:')[1].strip()
            base_name = os.path.splitext(filename)[0]
            if mission_name not in file_groups:
                file_groups[mission_name] = []
            file_groups[mission_name].append(base_name)

    for mission_name, base_names in file_groups.items():
        mission_folder = os.path.join(output_folder, mission_name)
        os.makedirs(mission_folder, exist_ok=True)
        for base_name in base_names:
            for filename in os.listdir(input_data_folder):
                if filename.startswith(base_name):
                    src_file = os.path.join(input_data_folder, filename)
                    shutil.move(src_file, mission_folder)
        if sys_log_exists:
            shutil.copy(sys_log_path, mission_folder)

def datasorter_organize_log_files(input_log_folder, output_folder):

    ''' 
    Organize log files by matching them with corresponding data files.
    
    Arguments:
    - input_log_folder (str): The directory containing the input log files.
    - output_folder (str): The directory where the organized log files will be saved.
      
    Returns:
    - None
    '''

    for filename in os.listdir(input_log_folder):
        file_path = os.path.join(input_log_folder, filename)
        if os.path.getsize(file_path) == 0 or os.path.getsize(file_path) < 3072:
            os.remove(file_path)
    data_file_pattern = re.compile(r'\b\w+-\d{4}-\d{3}-\d+-\d+\b')
    for filename in os.listdir(input_log_folder):
        if filename.endswith('.log'):
            with open(os.path.join(input_log_folder, filename), 'r', encoding='utf-8', errors='ignore') as file:
                log_content = file.read()
                match = data_file_pattern.search(log_content)
                if match:
                    data_filename = match.group()
                    for mission_name in os.listdir(output_folder):
                        mission_folder = os.path.join(output_folder, mission_name)
                        if os.path.isdir(mission_folder):
                            for data_file in os.listdir(mission_folder):
                                if data_file.startswith(data_filename):
                                    src_file = os.path.join(input_log_folder, filename)
                                    dest_file = os.path.join(mission_folder, filename)
                                    file.close()
                                    shutil.move(src_file, dest_file)
                                    break

def datasorter_run(root_directory, glider_info):

    ''' 
    Run the data sorting process to organize data and log files.
    
    Arguments:
    - root_directory (str): The root directory from the configuration file.
    - glider_info (tuple): A tuple containing glider unit, version, and type information.
      
    Returns:
    - None
    '''

    print(f"\n### RUNNING: DATA SORTER ###\n")
    datetime_string = datetime.now().strftime("%Y%m%dT%H%M%S")
    glider_unit, glider_version, glider_type = glider_info
    data_sort_directory = os.path.join(root_directory, f"DataSorter")
    runtime_directory = os.path.join(data_sort_directory, f"{glider_unit}-{glider_type}-{glider_version}-{datetime_string}")
    os.makedirs(runtime_directory, exist_ok=True)
    current_directory = os.path.abspath(os.path.dirname(__file__))
    logfile_directory = os.path.join(current_directory, 'DBD_Files', 'Logfiles')
    data_directory = os.path.join(current_directory, 'DBD_Files')
    datasorter_organize_data_files(data_directory, runtime_directory)
    datasorter_organize_log_files(logfile_directory, runtime_directory)
    print(f"\nGlider Data/Log Files saved to: {runtime_directory}\n")

# =========================
# LOGFILE SEARCH
# =========================

class logsearch(QDialog):

    ''' 
    Create a popup dialog to collect multiple logfile search strings from the user.
    
    Arguments:
    - parent (QWidget, optional): The parent widget. Defaults to None.
      
    Returns:
    - search_strings (list): List of search strings entered by the user.
    '''

    def __init__(self, parent=None):
        super().__init__(parent)
        self.search_strings = []
        self.filter_rows = []
        self.logsearch_init()
    
    def logsearch_init(self):

        ''' 
        Initialize the dialog layout.
        
        Arguments:
        - None
        
        Returns:
        - None
        '''

        self.setWindowTitle("Logfile Search")
        self.setMinimumSize(500, 300)
        
        self.layout = QVBoxLayout(self)
        
        self.scrollArea = QScrollArea()
        self.scrollArea.setWidgetResizable(True)
        self.scrollWidget = QWidget()
        self.scrollLayout = QVBoxLayout(self.scrollWidget)
        self.scrollWidget.setLayout(self.scrollLayout)
        self.scrollArea.setWidget(self.scrollWidget)
        self.layout.addWidget(self.scrollArea)
        
        self.addButton = QPushButton("+ Add Search String")
        self.addButton.clicked.connect(self.logsearch_add_search)
        self.layout.addWidget(self.addButton)
        
        self.runButton = QPushButton("Run Logfile Search")
        self.runButton.clicked.connect(self.logsearch_apply_search)
        self.layout.addWidget(self.runButton)
        
        self.logsearch_add_search()
    
    def logsearch_add_search(self):

        ''' 
        Add a new row containing a QLineEdit for entering a search string.
        
        Arguments:
        - None
        
        Returns:
        - None
        '''

        rowWidget = QWidget()
        rowLayout = QHBoxLayout(rowWidget)
        searchEdit = QLineEdit()
        searchEdit.setPlaceholderText("Enter search string")
        rowLayout.addWidget(searchEdit)
        self.scrollLayout.addWidget(rowWidget)
        self.filter_rows.append(searchEdit)
    
    def logsearch_apply_search(self):

        ''' 
        Gather non-empty search strings from all rows and close the dialog.
        
        Arguments:
        - None
        
        Returns:
        - search_strings (list): List of non-empty search strings.
        '''

        self.search_strings = []
        for searchEdit in self.filter_rows:
            text = searchEdit.text().strip()
            if text:
                self.search_strings.append(text)
        if not self.search_strings:
            QMessageBox.warning(self, "No Search Strings", "Please enter at least one search string.")
            return
        self.accept()

def logsearch_run(root_directory):

    ''' 
    Run the logfile search process.
    
    Arguments:
    - root_directory (str): The main output directory where grouped log files are saved.
      
    Returns:
    - None
    '''

    print(f"\n### RUNNING: LOGFILE SEARCH ###\n")

    dialog = logsearch()
    if dialog.exec() == QDialog.DialogCode.Accepted:
        search_strings = dialog.search_strings
    else:
        print("Logfile search cancelled by user.")
        return

    current_directory = os.path.abspath(os.path.dirname(__file__))
    logfiles_dir = os.path.join(current_directory, "DBD_Files", "Logfiles")
    if not os.path.exists(logfiles_dir):
        print(f"Logfiles directory not found: {logfiles_dir}")
        return

    for search_str in search_strings:
        output_folder = os.path.join(root_directory, f"LogfileSearch-{search_str}")
        os.makedirs(output_folder, exist_ok=True)
        print(f"Searching for '{search_str}' in log files...")
        for file_name in os.listdir(logfiles_dir):
            file_path = os.path.join(logfiles_dir, file_name)
            if os.path.isfile(file_path):
                try:
                    with open(file_path, 'r', encoding='utf-8', errors='ignore') as f:
                        content = f.read()
                    if search_str in content:
                        shutil.copy(file_path, output_folder)
                        print(f"Copied {file_name} to {output_folder}")
                except Exception as e:
                    print(f"Error reading {file_path}: {e}")
        print(f"Search for '{search_str}' complete. Files copied to {output_folder}.")

# =========================
# ENERGY EVALUATION
# =========================

class energy(QWidget):

    ''' 
    A widget for energy evaluation that integrates into the main GUI as a tab.
    The user selects the glider type and clicks "Run Energy Evaluation" to compute energy metrics.
    The results are both saved to a JSON file and displayed as gas gauges.
    
    Arguments:
    - root_directory (str): The directory where the output file will be saved.
    - glider_info (tuple): A tuple containing glider unit, version, and type information.
    - dataframe (pd.DataFrame): The input dataframe containing power variables.
    - parent (QWidget, optional): Parent widget.
      
    Returns:
    - None
    '''

    def __init__(self, root_directory, glider_info, dataframe, parent=None):
        super().__init__(parent)
        self.root_directory = root_directory
        self.glider_info = glider_info
        self.dataframe = dataframe.copy()
        self.results = {}
        self.energy_init()
    
    def energy_init(self):

        ''' 
        Initialize the energy evaluation tab layout with controls and gauge display area.
        
        Arguments:
        - None
        
        Returns:
        - None
        '''

        self.setMinimumSize(400, 300)
        main_layout = QVBoxLayout(self)
        
        control_layout = QHBoxLayout()
        self.type_label = QLabel("Select Glider Type:")
        control_layout.addWidget(self.type_label)
        
        self.comboBox = QComboBox()
        self.comboBox.addItems(["Slocum G3s (.dbd)", "Slocum G3s (.sbd)", "Slocum Sentinel (.dbd)", "Slocum Sentinel (.sbd)"])
        control_layout.addWidget(self.comboBox)
        
        self.runButton = QPushButton("Run Energy Evaluation")
        self.runButton.clicked.connect(self.energy_evaluation)
        control_layout.addWidget(self.runButton)
        main_layout.addLayout(control_layout)
        
        self.resultsArea = QWidget()
        self.resultsLayout = QVBoxLayout(self.resultsArea)
        self.resultsArea.setLayout(self.resultsLayout)
        main_layout.addWidget(self.resultsArea)
        
        self.setLayout(main_layout)
    
    def energy_evaluation(self):

        ''' 
        Run the energy evaluation based on the selected glider type and update the gauge display.
        
        Arguments:
        - None
        
        Returns:
        - None
        '''

        glider_type = self.comboBox.currentText()
        energy_dictionary = {}
        
        if glider_type == "Slocum G3s (.dbd)":
            required_columns = ['m_coulomb_amphr_total', 'm_coulomb_current', 'm_battery', 'time']
            if not all(col in self.dataframe.columns for col in required_columns):
                QMessageBox.warning(self, "Missing Data", "Required columns are missing in the dataframe for Slocum G3s.")
                return
            df = self.dataframe.copy()
            df['time'] = pd.to_datetime(df['time'])
            df.set_index('time', inplace=True)
            daily_amp_hr_usage = df['m_coulomb_amphr_total'].diff().resample('D').sum().mean()
            daily_watt_hr_usage = (df['m_coulomb_amphr_total'].diff() * df['m_battery']).resample('D').sum().mean()
            avg_coulomb_current = df['m_coulomb_current'].mean()
            energy_dictionary = {
                'Average Amp-Hour Usage per Day': daily_amp_hr_usage,
                'Watt-Hour Usage per Day': daily_watt_hr_usage,
                'Average Coulomb Current': avg_coulomb_current
            }
        
        elif glider_type == "Slocum Sentinel (.dbd)":
            required_columns = [
                'm_bms1_batt_pack_0_inst_current', 'm_bms1_batt_pack_1_inst_current',
                'm_bms2_batt_pack_0_inst_current', 'm_bms2_batt_pack_1_inst_current',
                'm_bms3_batt_pack_0_inst_current', 'm_bms3_batt_pack_1_inst_current',
                'm_battery', 'm_battery_amp_hours_remaining', 'm_battery_watt_hours_remaining'
            ]
            if not all(col in self.dataframe.columns for col in required_columns):
                QMessageBox.warning(self, "Missing Data", "Required columns are missing in the dataframe for Slocum Sentinel.")
                return
            df = self.dataframe.copy()
            df['time'] = pd.to_datetime(df['time'])
            df.set_index('time', inplace=True)
            daily_amp_hr_usage = df['m_battery_amp_hours_remaining'].diff().resample('D').sum().mean()
            daily_watt_hr_usage = df['m_battery_watt_hours_remaining'].diff().resample('D').sum().mean()
            avg_bms1_pack0 = df['m_bms1_batt_pack_0_inst_current'].mean()
            avg_bms1_pack1 = df['m_bms1_batt_pack_1_inst_current'].mean()
            avg_bms2_pack0 = df['m_bms2_batt_pack_0_inst_current'].mean()
            avg_bms2_pack1 = df['m_bms2_batt_pack_1_inst_current'].mean()
            avg_bms3_pack0 = df['m_bms3_batt_pack_0_inst_current'].mean()
            avg_bms3_pack1 = df['m_bms3_batt_pack_1_inst_current'].mean()
            energy_dictionary = {
                'Average Amp-Hour Usage per Day': daily_amp_hr_usage,
                'Watt-Hour Usage per Day': daily_watt_hr_usage,
                'Average BMS1 Pack0 Current': avg_bms1_pack0,
                'Average BMS1 Pack1 Current': avg_bms1_pack1,
                'Average BMS2 Pack0 Current': avg_bms2_pack0,
                'Average BMS2 Pack1 Current': avg_bms2_pack1,
                'Average BMS3 Pack0 Current': avg_bms3_pack0,
                'Average BMS3 Pack1 Current': avg_bms3_pack1,
            }
        else:
            QMessageBox.warning(self, "Unknown Glider Type", "The selected glider type is not recognized.")
            return
        
        energy_json(self.root_directory, self.glider_info, energy_dictionary)
        
        for i in reversed(range(self.resultsLayout.count())):
            widget_to_remove = self.resultsLayout.itemAt(i).widget()
            self.resultsLayout.removeWidget(widget_to_remove)
            widget_to_remove.setParent(None)
        
        for metric, value in energy_dictionary.items():
            label = QLabel(f"{metric}: {value:.2f}")
            self.resultsLayout.addWidget(label)

def energy_json(root_directory, glider_info, data_dict):

    ''' 
    Save the energy evaluation results to a JSON file.
    
    Arguments:
    - root_directory (str): The directory where the output file will be saved.
    - glider_info (tuple): A tuple containing glider unit, version, and type information.
    - data_dict (dict): The dictionary containing energy evaluation results.
      
    Returns:
    - None
    '''

    print("\n### RUNNING: JSON ###\n")
    glider_unit, glider_version, glider_type = glider_info
    output_file = os.path.join(root_directory, f"{glider_unit}-{glider_version}-{glider_type}_EnergyOutput.json")
    json_data = {glider_unit: data_dict}
    with open(output_file, 'w') as json_file:
        json.dump(json_data, json_file, indent=4)
    print(f"All data saved to {output_file}")


def energy_run(root_directory, glider_info, dataframe):

    ''' 
    Create and return an instance of energy with the provided parameters.
    This widget is designed to be embedded as a tab in the main GUI.
    
    Arguments:
    - root_directory (str): The directory where the output file will be saved.
    - glider_info (tuple): A tuple containing glider unit, version, and type information.
    - dataframe (pd.DataFrame): The input dataframe containing power variables.
      
    Returns:
    - energy (energy): The energy evaluation widget.
    '''

    return energy(root_directory, glider_info, dataframe)
