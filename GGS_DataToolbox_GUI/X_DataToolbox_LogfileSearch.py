"""
Author: Salvatore Fricano
Repository: https://github.com/STF33/Glider-Guidance-System
"""

# =========================
# IMPORTS
# =========================

import os
import shutil
from PyQt6.QtWidgets import QDialog, QVBoxLayout, QHBoxLayout, QLineEdit, QPushButton, QLabel, QScrollArea, QWidget, QMessageBox
from PyQt6.QtCore import Qt

# =========================

class LogfileSearchDialog(QDialog):

    '''
    Create a popup dialog to collect multiple logfile search strings from the user.

    Args:
    - parent (QWidget, optional): The parent widget. Defaults to None.
      
    Returns:
    - search_strings (list): List of search strings entered by the user.
    '''


    def __init__(self, parent=None):
        super().__init__(parent)
        self.search_strings = []
        self.filter_rows = []
        self.initUI()
    
    def initUI(self):
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
        self.addButton.clicked.connect(self.add_search_row)
        self.layout.addWidget(self.addButton)
        
        self.runButton = QPushButton("Run Logfile Search")
        self.runButton.clicked.connect(self.apply_search)
        self.layout.addWidget(self.runButton)
        
        self.add_search_row()
    
    def add_search_row(self):

        '''
        Add a new row containing a QLineEdit for entering a search string.

        Args:
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
    
    def apply_search(self):

        '''
        Gather non-empty search strings from all rows and close the dialog.

        Args:
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

def run_logfile_search(root_directory):
    
    '''
    Run the logfile search process.
 
    Args:
    - root_directory (str): The main output directory where grouped log files are saved.
      
    Returns:
    - None
    '''

    dialog = LogfileSearchDialog()
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
