"""
Author: Salvatore Fricano
Repository: https://github.com/STF33/Glider-Guidance-System
"""

# =========================
# IMPORTS
# =========================

import os
import pandas as pd
from openpyxl import load_workbook

# =========================

### FUNCTION:
def run_excel(root_directory, glider_info, data_frame):
    
    '''
    Save the DataFrame to an Excel file.
    
    Args:
    - data_frame (pd.DataFrame): The DataFrame containing sensor data.
    - root_directory (str): The root directory from the configuration file.
    - glider_info (tuple): A tuple containing glider unit, version, and type information.

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