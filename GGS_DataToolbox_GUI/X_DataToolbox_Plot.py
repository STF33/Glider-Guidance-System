"""
Author: Salvatore Fricano
Repository: https://github.com/STF33/Glider-Guidance-System
"""

# =========================
# IMPORTS
# =========================

import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from PyQt6.QtWidgets import (QApplication, QWidget, QVBoxLayout, QHBoxLayout, QLabel, QPushButton, QComboBox, QListWidget, QListWidgetItem, QMessageBox, QSplitter)
from PyQt6.QtCore import Qt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

# =========================

### CLASS:
class InteractivePlotGUI(QWidget):

    ### FUNCTION:    
    def __init__(self, dataframe):
        super().__init__()
        self.dataframe = dataframe
        self.selected_variables = []
        self.x_variable = None
        self.initUI()

    ### FUNCTION:    
    def initUI(self):
        possible_x_vars = ['m_present_time', 'sci_m_present_time']
        for var in possible_x_vars:
            if var in self.dataframe.columns:
                self.x_variable = var
                break
        if self.x_variable is None:
            QMessageBox.warning(self, "Missing Time Variable",
                                "Required time variable (m_present_time or sci_m_present_time) not found in the dataframe.")
            self.close()
            return

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
        variables = [col for col in self.dataframe.columns if col != self.x_variable]
        self.variable_selector.addItems(variables)
        top_control_layout.addWidget(QLabel("Select Variable:"))
        top_control_layout.addWidget(self.variable_selector)
        
        self.add_remove_button = QPushButton("Add/Remove")
        self.add_remove_button.clicked.connect(self.toggle_variable)
        top_control_layout.addWidget(self.add_remove_button)
        control_layout.addLayout(top_control_layout)
        
        self.variable_list = QListWidget()
        control_layout.addWidget(self.variable_list)
        
        self.plot_button = QPushButton("Update Plot")
        self.plot_button.clicked.connect(self.update_plot)
        control_layout.addWidget(self.plot_button)
        
        splitter.addWidget(control_widget)
        splitter.setStretchFactor(1, 1)
        
        main_layout.addWidget(splitter)
        self.setLayout(main_layout)
        self.setWindowTitle("Interactive Plot GUI")

    ### FUNCTION:    
    def toggle_variable(self):
        var = self.variable_selector.currentText()
        if var in self.selected_variables:
            self.selected_variables.remove(var)
            for i in range(self.variable_list.count()):
                if self.variable_list.item(i).text() == var:
                    self.variable_list.takeItem(i)
                    break
        else:
            self.selected_variables.append(var)
            self.variable_list.addItem(QListWidgetItem(var))

    ### FUNCTION:    
    def update_plot(self):
        self.figure.clf()
        if not self.selected_variables:
            self.canvas.draw()
            return
        
        host = self.figure.add_subplot(111)
        
        x_data = pd.to_datetime(self.dataframe[self.x_variable].astype(float), unit='s')
        host.set_xlabel(self.x_variable)
        
        axes = [host]
        colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        
        variable_handles = []
        
        for idx, var in enumerate(self.selected_variables):
            y_data = pd.to_numeric(self.dataframe[var], errors='coerce').dropna()
            x_data_cleaned = x_data[y_data.index]
            
            if var == 'm_depth':
                ax = host
                ax.invert_yaxis()
                color = colors[idx % len(colors)]
                handle = ax.scatter(x_data_cleaned, y_data, color=color, label=var, marker='o')
                ax.plot(x_data_cleaned, y_data, color=color)
                variable_handles.append(handle)

                y_min = y_data.min()
                y_max = y_data.max()
                buffer = 0.05 * (y_max - y_min) if y_max != y_min else 1
                ax.set_ylim(y_max + buffer, y_min - buffer)
                ax.set_ylabel(var, color=color)
                ax.tick_params(axis='y', colors=color)

            else:
                ax = host.twinx()
                ax.spines["right"].set_position(("axes", 1 + 0.1*(len(axes)-1)))
                axes.append(ax)

                color = colors[idx % len(colors)]
                handle = ax.scatter(x_data_cleaned, y_data, color=color, label=var, marker='o')
                ax.plot(x_data_cleaned, y_data, color=color)
                variable_handles.append(handle)

                y_min = y_data.min()
                y_max = y_data.max()
                buffer = 0.05 * (y_max - y_min) if y_max != y_min else 1
                ax.set_ylim(y_min - buffer, y_max + buffer)
                ax.set_ylabel(var, color=color)
                ax.tick_params(axis='y', colors=color)

        
        host.xaxis_date()
        locator = mdates.AutoDateLocator()
        formatter = mdates.ConciseDateFormatter(locator)
        host.xaxis.set_major_locator(locator)
        host.xaxis.set_major_formatter(formatter)

        host.legend(handles=variable_handles, labels=[handle.get_label() for handle in variable_handles],
                    loc="upper center", bbox_to_anchor=(0.5, 1.15),
                    ncol=len(variable_handles), borderaxespad=0)

        self.canvas.draw()

plot_window = None

### FUNCTION:
def run_plot(dataframe):
    global plot_window
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    plot_window = InteractivePlotGUI(dataframe)
    plot_window.show()
