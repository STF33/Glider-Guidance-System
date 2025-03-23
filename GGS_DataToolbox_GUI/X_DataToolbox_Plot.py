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
from PyQt6.QtGui import QColor
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

# =========================

### CLASS:
class InteractivePlotGUI(QWidget):

    """
    Create an interactive plot GUI using PyQt and Matplotlib.
    
    Args:
      dataframe (pd.DataFrame): The merged sensor data containing a "time" column and sensor values.
      
    Returns:
      None
    """

    def __init__(self, dataframe):
        super().__init__()
        self.dataframe = dataframe
        self.selected_variables = []
        self.variable_colors = {}
        self.available_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
        self.initUI()
    
    ### FUNCTION:
    def initUI(self):

        """
        Initialize the interactive plot GUI layout, including the plot area and control panel.
        
        Args:
          None
          
        Returns:
          None
        """

        if "time" not in self.dataframe.columns:
            QMessageBox.warning(self, "Missing Time Column", "The merged dataframe must have a 'time' column.")
            self.close()
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

      """
      Toggle the inclusion of a variable for plotting when the user clicks the "Add/Remove" button.
      
      Args:
        None
        
      Returns:
        None
      """

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

    ### FUNCTION:
    def update_plot(self):

        """
        Update the plot based on the selected variables, grouping y-axes by unit and adjusting layout.
        
        Args:
          None
          
        Returns:
          None
        """

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
        host.legend(handles=legend_handles, labels=legend_labels, loc="upper center", bbox_to_anchor=(0.5, 1.15),
                    ncol=len(legend_handles), borderaxespad=0, frameon=True)
        
        n_extra = len(axes_list) - 1
        new_right = max(right_margin - 0.05 * n_extra, 0.5)
        self.figure.subplots_adjust(left=left_margin, right=new_right, top=top_margin, bottom=bottom_margin)
        
        self.canvas.draw()

plot_window = None

### FUNCTION:
def run_plot(dataframe):

    """
    Launch the interactive plot GUI with the provided dataframe.

    Args:
      dataframe (pd.DataFrame): The merged sensor data containing a "time" column and sensor values.
      
    Returns:
      None
    """

    global plot_window
    app = QApplication.instance()
    if app is None:
        app = QApplication(sys.argv)
    plot_window = InteractivePlotGUI(dataframe)
    plot_window.show()
    
