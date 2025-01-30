"""
Author: Salvatore Fricano
Repository: https://github.com/STF33/Glider-Guidance-System
"""

# =========================
# IMPORTS
# =========================

import pandas as pd
import os
from openpyxl import load_workbook

# =========================

### FUNCTION:
def calculate_statistics(df, pump_value):

    '''
    Calculates statistics for the low power data evaluation
    
    Args:
    - df (pd.DataFrame): The dataframe interval of glider data to evaluate.
    - pump_value (str): The column name for the pump data (either 'm_ballast_pumped' or 'm_de_oil_vol').

    Returns:
    - glider_statistics (dict): A dictionary containing the calculated statistics.
    '''
    
    initial_time = df['m_present_secs_into_mission'].iloc[0]
    final_time = df['m_present_secs_into_mission'].iloc[-1]
    initial_amp_hr = df['m_coulomb_amphr_total'].iloc[0]
    final_amp_hr = df['m_coulomb_amphr_total'].iloc[-1]
    
    time_change = final_time - initial_time
    amp_hr_change = final_amp_hr - initial_amp_hr
    drain_rate = (amp_hr_change / time_change) * 86400 if time_change != 0 else float('nan')
    avg_current = df['m_coulomb_current'].mean()
    min_depth = df['m_depth'].min()
    max_depth = df['m_depth'].max()
    
    delta_ccs_pumped = df[pump_value].diff().abs().sum()
    ratio_amp_hr_to_ccs_pumped = amp_hr_change / delta_ccs_pumped if delta_ccs_pumped != 0 else float('nan')
    
    return {
        'initial_time': initial_time,
        'final_time': final_time,
        'time_change': time_change,
        'initial_amp_hr': initial_amp_hr,
        'final_amp_hr': final_amp_hr,
        'amp_hr_change': amp_hr_change,
        'drain_rate': drain_rate,
        'avg_current': avg_current,
        'min_depth': min_depth,
        'max_depth': max_depth,
        'delta_ccs_pumped': delta_ccs_pumped,
        'ratio_amp_hr_to_ccs_pumped': ratio_amp_hr_to_ccs_pumped
    }

### FUNCTION:
def run_yo_eval(root_directory, glider_info):

    '''
    Executes the yo data evaluation.
    
    Args:
    - root_directory (str): The root directory from the configuration file.
    - glider_info (tuple): A tuple containing glider unit, version, and type information.

    Returns:
    - None
    '''

    print(f"\n### RUNNING: YO EVALUATION ###\n")

    glider_unit, glider_version, glider_type = glider_info

    df_filename = f"{glider_unit}-{glider_version}-{glider_type}_DataOutput.xlsx"
    df_file = os.path.join(root_directory, df_filename)

    df = pd.read_excel(df_file, engine='openpyxl')

    if 'm_ballast_pumped' in df.columns:
        pump_value = 'm_ballast_pumped'
    elif 'm_de_oil_vol' in df.columns:
        pump_value = 'm_de_oil_vol'
    else:
        raise ValueError("Unable to determine pump type from DataOutput file.")
    
    yo_intervals = []
    current_interval = []
    in_dive = False
    in_climb = False
    buffer_value = 20

    for i in range(1, len(df)):
        if not in_dive and df[pump_value].iloc[i] < df[pump_value].iloc[i-1] - buffer_value:
            in_dive = True
            current_interval.append(df.iloc[i-1])
            current_interval.append(df.iloc[i])
        elif in_dive and not in_climb and df[pump_value].iloc[i] > df[pump_value].iloc[i-1] + buffer_value:
            in_climb = True
            current_interval.append(df.iloc[i])
        elif in_climb and (df[pump_value].iloc[i] < df[pump_value].iloc[i-1] - buffer_value or df['c_air_pump'].iloc[i] == 1):
            current_interval.append(df.iloc[i])
            yo_intervals.append(pd.DataFrame(current_interval))
            current_interval = []
            in_dive = False
            in_climb = False
        elif in_dive or in_climb:
            current_interval.append(df.iloc[i])

    yo_stats = []

    for interval in yo_intervals:
        stats = calculate_statistics(interval, pump_value)
        yo_stats.append(stats)

    output_file = os.path.join(root_directory, f'{glider_unit}-{glider_version}_{glider_type}_YoStatistics.xlsx')
    stats_df = pd.DataFrame(yo_stats)
    stats_df.to_excel(output_file, index=False)

    workbook = load_workbook(output_file)
    worksheet = workbook.active
    for col in worksheet.columns:
        worksheet.column_dimensions[col[0].column_letter].width = 20
    workbook.save(output_file)

    print(f"Statistics calculated and saved to {output_file}")