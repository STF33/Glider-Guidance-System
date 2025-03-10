
# Visual Studio Code Installation and Environment Setup Guide

## Overview
This guide provides instructions for installing Visual Studio Code, setting up Python/Conda, and importing an existing Conda environment from a `.yml` file.

## Installation Section
This section outlines the steps to install Visual Studio Code and the necessary tools for Python.
- **Visual Studio Code**: Download and install Visual Studio Code from the [official website](https://code.visualstudio.com/). Follow the installation prompts suitable for your operating system (Windows, macOS, Linux).
- **Python**: Install Python by downloading it from the [official Python website](https://www.python.org/downloads/). Ensure to check the box "Add to PATH" during installation.
- **Conda**: Install Conda by by downloading it from [Anaconda](https://www.anaconda.com/download/success). Ensure to check the box "Add to PATH" during installation.

## Extension Section
This section helps you set up Visual Studio Code for Python development.
- **Python Extension for Visual Studio Code**: Open Visual Studio Code, go to Extensions, and search for `Python`. Install the extension published by Microsoft.

## Creating an Environment
- **Creating a Conda Environment**: Open your terminal (this can be done in Visual Studio Code at the top of the window by pressing: Terminal --> New Terminal) and create a new environment by running `conda create --name myenv`, replacing `myenv` with your desired environment name. For the Glider Guidance System: Data Toolbox, it is recommended to simply name the enfironment `ggs_dt`. This would make the full command to type: `conda create --name ggs_dt`.
- **Activating the Environment**: Activate the newly created environment by running `conda activate ggs_dt` (or whatever the environment name you wish to activate is saved as).
- Install packages for the environmeny using `pip install`. The packages you will need to install are:
    1. PyQt6 (`pip install PyQt6`)
    2. openpyxl (`pip install openpyxl`)
    3. pandas (`pip install pandas`)
    4. matplotlib (`pip install matplotlib`)
    5. numpy (`pip install numpy`)
- Type all of the following commands into your terminal while the environment is active.
