"""
Author: Salvatore Fricano
Repository: https://github.com/STF33/Glider-Guidance-System
"""

# =========================
# IMPORTS
# =========================

import pandas as pd
import os

# =========================

### FUNCTION:
def calculate_statistics(sub_df):
    
    '''
    Calculates statistics for the low power data evaluation
    
    Args:
    - sub_df (pd.DataFrame): The dataframe interval of glider data to evaluate.

    Returns:
    - glider_statistics (dict): A dictionary containing the calculated statistics.
    '''

    initial_time = sub_df['m_present_secs_into_mission'].iloc[0]
    final_time = sub_df['m_present_secs_into_mission'].iloc[-1]
    initial_amp_hr = sub_df['m_coulomb_amphr_total'].iloc[0]
    final_amp_hr = sub_df['m_coulomb_amphr_total'].iloc[-1]
    
    time_change = final_time - initial_time
    amp_hr_change = final_amp_hr - initial_amp_hr
    drain_rate = (amp_hr_change / time_change) * 86400 if time_change != 0 else float('nan')
    avg_current = sub_df['m_coulomb_current'].mean()
    min_depth = sub_df['m_depth'].min()
    max_depth = sub_df['m_depth'].max()
    
    return {
        'initial_time': initial_time,
        'final_time': final_time,
        'time_change': time_change,
        'initial_amp_hr': initial_amp_hr,
        'final_amp_hr': final_amp_hr,
        'amp_hr_change': amp_hr_change,
        'drain_rate': drain_rate,
        'avg_current': avg_current,
        'sub_category': sub_df['sub_category'].iloc[0],
        'min_depth': min_depth,
        'max_depth': max_depth
    }

### FUNCTION:
def run_low_power_eval(root_directory, glider_info):

    '''
    Executes the low power data evaluation.
    
    Args:
    - root_directory (str): The root directory from the configuration file.
    - glider_info (tuple): A tuple containing glider unit, version, and type information.

    Returns:
    - None
    '''

    print(f"\n### RUNNING: LOW POWER EVALUATION ###\n")

    glider_unit, glider_version, glider_type = glider_info

    df_filename = f"{glider_unit}-{glider_version}-{glider_type}_DataOutput.xlsx"
    df_file = os.path.join(root_directory, df_filename)

    df = pd.read_excel(df_file)

    low_power_df = df[df['x_low_power_status'] == 0].copy()

    if glider_type == "shallow":
        low_power_df['sub_category'] = low_power_df['m_ballast_pumped'].apply(lambda x: 'diving' if x < 0 else 'climbing')
    if glider_type == "deep":
        low_power_df['sub_category'] = low_power_df['m_de_oil_vol'].apply(lambda x: 'diving' if x < 0 else 'climbing')

    low_power_intervals = []
    current_interval = []

    for i in range(len(low_power_df)):
        if i == 0 or low_power_df.index[i] == low_power_df.index[i-1] + 1:
            current_interval.append(low_power_df.iloc[i])
        else:
            low_power_intervals.append(pd.DataFrame(current_interval))
            current_interval = [low_power_df.iloc[i]]

    if current_interval:
        low_power_intervals.append(pd.DataFrame(current_interval))

    diving_stats = []
    climbing_stats = []

    for interval in low_power_intervals:
        stats = calculate_statistics(interval)
        if interval['sub_category'].iloc[0] == 'diving':
            diving_stats.append(stats)
        else:
            climbing_stats.append(stats)

    output_file = os.path.join(root_directory, f'{glider_unit}-{glider_version}_{glider_type}_LowPowerStatistics.xlsx')
    stats_df = pd.DataFrame(diving_stats + climbing_stats)
    stats_df.to_excel(output_file, index=False)

    print(f"Statistics calculated and saved to {output_file}")
