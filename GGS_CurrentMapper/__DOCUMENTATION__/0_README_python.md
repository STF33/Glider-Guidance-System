
# Visual Studio Code Installation and Environment Setup Guide

## Overview
This guide provides instructions for installing Visual Studio Code on any computer, setting up Python and Conda, and importing an existing Conda environment from a `.yml` file.

## Installation Section
This section outlines the steps to install Visual Studio Code and the necessary tools for Python.
- **Visual Studio Code**: Download and install Visual Studio Code from the [official website](https://code.visualstudio.com/). Follow the installation prompts suitable for your operating system (Windows, macOS, Linux).
- **Python**: Install Python by downloading it from the [official Python website](https://www.python.org/downloads/). Ensure to check the box "Add to PATH" during installation.
- **Conda**: Install Conda by choosing either Anaconda or Miniconda. Anaconda includes a suite of pre-installed packages suitable for scientific computing and data science. Miniconda is a minimal installer for Conda. Download from [Anaconda](https://www.anaconda.com/products/individual) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) websites. Ensure to check the box "Add to PATH" during installation.

## Extension Section
This section helps you set up Visual Studio Code for Python development.
- **Python Extension for Visual Studio Code**: Open Visual Studio Code, go to Extensions, and search for `Python`. Install the extension published by Microsoft.

## Environment Setup Section
After installing the necessary tools, you can setup a new environment in Python using Conda.
- **Creating a Conda Environment**: Open your terminal (this can be done in Visual Studio Code at the top of the window by pressing: Terminal --> New Terminal) and create a new environment by running `conda create --name myenv`, replacing `myenv` with your desired environment name and `3.12` with your preferred Python version.
- **Activating the Environment**: Activate the newly created environment by running `conda activate myenv`.

## Installing Packages to Environment
If you have a `.yml` file specifying an environment, you can easily import it.
- **Importing the Environment**: Ensure the `.yml` file is accessible on your computer. Open your terminal (this can be done in Visual Studio Code at the top of the window by pressing: Terminal --> New Terminal), navigate to the directory containing the `.yml` file, and run `conda env create -f environment.yml`, replacing `environment.yml` with the name of your file.
If you are working with a new environment, use the following syntax to install packages.
- **Installing Packages**: Ensure that your environment is active by typing `conda activate myenv`, replacing `myenv` with the name of your environment. Next, type in `conda install -c conda-forge` followed by the package names you want to install such as `pandas` or `xarray`. An example of this command in full may look like: `conda install -c conda-forge pandas xarray`.
