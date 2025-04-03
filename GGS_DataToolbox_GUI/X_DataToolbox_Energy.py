"""
Author: Salvatore Fricano
Repository: https://github.com/STF33/Glider-Guidance-System
"""

# =========================
# IMPORTS
# =========================

import pandas as pd
import os
import json
from PyQt6.QtWidgets import QDialog, QVBoxLayout, QLabel, QComboBox, QPushButton, QApplication

# =========================

### CLASS:
class EnergyEvaluation(QDialog):

    '''
    Create a dialog for running energy evaluation based on glider type selection.

    Args:
    - dataframe (pd.DataFrame): The input dataframe containing power variables.
    - root_directory (str): The directory where the output file will be saved.
    - parent (QWidget, optional): The parent widget. Defaults to None.

    Returns:
    - None
    '''
    
    def __init__(self, root_directory, glider_info, dataframe, parent=None):
        super().__init__(parent)
        self.root_directory = root_directory
        self.glider_info = glider_info
        self.dataframe = dataframe.copy()
        self.initUI()
    
    def initUI(self):

        '''
        Initialize the energy evaluation dialog layout with a combo box and button.

        Args:
        - None

        Returns:
        - None
        '''

        self.setWindowTitle("Energy Evaluation")
        self.setMinimumSize(300, 150)
        self.layout = QVBoxLayout(self)
        
        self.label = QLabel("Select Glider Type:")
        self.layout.addWidget(self.label)
        
        self.comboBox = QComboBox()
        self.comboBox.addItems(["Slocum G3s", "Slocum Sentinel"])
        self.layout.addWidget(self.comboBox)
        
        self.runButton = QPushButton("Run Energy Evaluation")
        self.runButton.clicked.connect(self.run_evaluation)
        self.layout.addWidget(self.runButton)
    
    def run_evaluation(self):

        '''
        Run the energy evaluation based on the selected glider type.

        Args:
        - None

        Returns:
        - None
        '''

        glider_type = self.comboBox.currentText()
        
        if glider_type == "Slocum G3s":
            required_columns = ['m_coulomb_amphr_total', 'm_coulomb_current', 'm_battery', 'time']
    
            if not all(col in self.dataframe.columns for col in required_columns):
                print("Warning: Required columns are missing in the dataframe.")
                return self.dataframe
    
            self.dataframe['time'] = pd.to_datetime(self.dataframe['time'])
            self.dataframe.set_index('time', inplace=True)

            daily_amp_hr_usage = self.dataframe['m_coulomb_amphr_total'].diff().resample('D').sum().mean()
            
            daily_watt_hr_usage = (self.dataframe['m_coulomb_amphr_total'].diff() * self.dataframe['m_battery']).resample('D').sum().mean()
            
            avg_coulomb_current = self.dataframe['m_coulomb_current'].mean()
            
            results_dict = {
                'Average Amp-Hour Usage per Day': daily_amp_hr_usage,
                'Watt-Hour Usage per Day': daily_watt_hr_usage,
                'Average Coulomb Current': avg_coulomb_current
            }
            
            run_json(self.root_directory, self.glider_info, results_dict)
        
        self.accept()

def run_json(root_directory, glider_info, data_dict):
    
    '''
    Save the energy evaluation results to a JSON file.

    Args:
    - root_directory (str): The directory where the output file will be saved.
    - glider_info (tuple): A tuple containing glider unit, version, and type information.
    - data_dict (dict): The dictionary containing energy evaluation results.

    Returns:
    - None
    '''

    print(f"\n### RUNNING: JSON ###\n")
    
    glider_unit, glider_version, glider_type = glider_info
    output_file = os.path.join(root_directory, f"{glider_unit}-{glider_version}-{glider_type}_EnergyOutput.json")
    
    print(glider_info)

    json_data = {
        glider_unit: data_dict
    }
    
    with open(output_file, 'w') as json_file:
        json.dump(json_data, json_file, indent=4)

    print(f"All data saved to {output_file}")

def run_energy_evaluation(root_directory, glider_info, dataframe):

    '''
    Run the energy evaluation dialog on the input dataframe.

    Args:
    - dataframe (pd.DataFrame): The input dataframe containing power variables.
    - root_directory (str): The directory where the output file will be saved.

    Returns:
    - dataframe (pd.DataFrame): The dataframe after running the energy evaluation.
    '''

    app = QApplication.instance() or QApplication([])
    dialog = EnergyEvaluation(root_directory, glider_info, dataframe)
    if dialog.exec() == QDialog.DialogCode.Accepted:
        return dialog.dataframe
    else:
        return dataframe
