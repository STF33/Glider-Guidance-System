"""
Author: Salvatore Fricano
Repository: https://github.com/STF33/Glider-Guidance-System

This module provides functions for parsing ASCII data files and compiling sensor data into a DataFrame.
"""

# =========================
# IMPORTS
# =========================

import pandas as pd
from PyQt6.QtWidgets import QDialog, QVBoxLayout, QHBoxLayout, QLineEdit, QComboBox, QPushButton, QScrollArea, QWidget

# =========================

OPERATORS = {
    '<': lambda series, val: series < val,
    '<=': lambda series, val: series <= val,
    '=': lambda series, val: series == val,
    '>=': lambda series, val: series >= val,
    '>': lambda series, val: series > val
}

### CLASS:
class DataFilter(QDialog):

    '''
    Create a dialog for adding multiple data filter criteria.
    
    Args:
    - dataframe (pd.DataFrame): The input dataframe to be filtered.
      
    Returns:
    - None
    '''

    def __init__(self, dataframe, parent=None):
        super().__init__(parent)
        self.dataframe = dataframe.copy()
        self.filter_rows = []
        self.filtered_df = None
        self.initUI()
    
    def initUI(self):

        '''
        Initialize the filter dialog layout with a scroll area and buttons.
        
        Args:
        - None
          
        Returns:
        - None
        '''

        self.setWindowTitle("Data Filter")
        self.setMinimumSize(500, 300)
        self.layout = QVBoxLayout(self)
        
        self.filterWidget = QWidget()
        self.filterLayout = QVBoxLayout(self.filterWidget)
        self.filterWidget.setLayout(self.filterLayout)
        
        self.scroll = QScrollArea()
        self.scroll.setWidgetResizable(True)
        self.scroll.setWidget(self.filterWidget)
        self.layout.addWidget(self.scroll)
        
        self.addFilterButton = QPushButton("+ Add Filter")
        self.addFilterButton.clicked.connect(self.add_filter_row)
        self.layout.addWidget(self.addFilterButton)
        
        self.runFilterButton = QPushButton("Run Data Filter")
        self.runFilterButton.clicked.connect(self.apply_filters)
        self.layout.addWidget(self.runFilterButton)
    
    def add_filter_row(self):

        '''
        Add a new filter row to the dialog. Each row contains a QLineEdit for the variable, a QComboBox for the operator, and a QLineEdit for the value.
        
        Args:
        - None
          
        Returns:
        - None
        '''

        row_widget = QWidget()
        row_layout = QHBoxLayout(row_widget)
        
        var_edit = QLineEdit()
        var_edit.setPlaceholderText("Variable")
        row_layout.addWidget(var_edit)
        
        op_combo = QComboBox()
        op_combo.addItems(list(OPERATORS.keys()))
        row_layout.addWidget(op_combo)
        
        val_edit = QLineEdit()
        val_edit.setPlaceholderText("Value")
        row_layout.addWidget(val_edit)
        
        self.filterLayout.addWidget(row_widget)
        self.filter_rows.append((var_edit, op_combo, val_edit))
    
    def apply_filters(self):

        '''
        Apply all filter rows to the dataframe and close the dialog.
        
        Args:
        - None
          
        Returns:
        - None
        '''

        df_filtered = self.dataframe.copy()
        for var_edit, op_combo, val_edit in self.filter_rows:
            var = var_edit.text().strip()
            op = op_combo.currentText().strip()
            val_str = val_edit.text().strip()
            if var == "" or val_str == "":
                continue
            try:
                try:
                    val = float(val_str)
                except ValueError:
                    val = val_str
                if var not in df_filtered.columns:
                    print(f"Warning: Variable '{var}' not found in dataframe. Skipping filter.")
                    continue
                condition = OPERATORS[op](pd.to_numeric(df_filtered[var], errors='coerce'), val)
                df_filtered = df_filtered[condition]
            except Exception as e:
                print(f"Error applying filter on {var} {op} {val_str}: {e}")
        self.filtered_df = df_filtered
        self.accept()

### FUNCTION:
def run_data_filter(dataframe):

    '''
    Run the data filter dialog on the input dataframe.
    
    Args:
    - dataframe (pd.DataFrame): The input dataframe to filter.
      
    Returns:
    - filtered_df (pd.DataFrame): The filtered dataframe after applying user-selected criteria.
    '''

    dialog = DataFilter(dataframe)
    if dialog.exec() == QDialog.DialogCode.Accepted:
        return dialog.filtered_df
    else:
        return dataframe
