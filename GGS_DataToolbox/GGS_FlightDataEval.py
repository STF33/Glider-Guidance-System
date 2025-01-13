# =========================
# IMPORTS
# =========================

import os
import pandas as pd
import matplotlib.pyplot as plt

# =========================

def approve(message):
    """Function to get user approval"""
    response = input(f"{message} (y/n): ")
    return response.lower() == 'y'

def run_fde(df):

    print(f"\n### RUNNING: FLIGHT DATA EVAL ###\n")

    df['timestamp'] = pd.to_datetime(df['m_present_time'], unit='s')

    # Test 1: Average time difference between m_depth datapoints
    avg_time_diff = df['timestamp'].diff().mean()
    if not approve(f"Test 1: Average m_depth cycle time = {avg_time_diff}. Approve?"):
        return

    # Test 2: m_de_oil_vol vs c_de_oil_vol
    if 'm_de_oil_vol' in df.columns and 'c_de_oil_vol' in df.columns:
        plt.figure()
        plt.scatter(df['timestamp'], df['m_de_oil_vol'], label='m_de_oil_vol')
        plt.scatter(df['timestamp'], df['c_de_oil_vol'], label='c_de_oil_vol')
        plt.xlabel('Time (secs)')
        plt.ylabel('Oil Volume')
        plt.legend()
        plt.title('Oil Volume over Time')
        plt.show()
        if not approve("Test 2: m_de_oil_vol and c_de_oil_vol versus time. Approve?"):
            return

    # Test 3: Plot m_pitch and c_pitch over time
    if 'm_pitch' in df.columns and 'c_pitch' in df.columns:
        plt.figure()
        plt.scatter(df['timestamp'], df['m_pitch'], label='m_pitch')
        plt.scatter(df['timestamp'], df['c_pitch'], label='c_pitch')
        plt.xlabel('Time (secs)')
        plt.ylabel('Pitch')
        plt.legend()
        plt.title('Pitch over Time')
        plt.show()
        if not approve("Test 3: m_pitch and c_pitch versus time. Approve?"):
            return

    # Test 4: Plot m_battpos and c_battpos over time
    if 'm_battpos' in df.columns and 'c_battpos' in df.columns:
        plt.figure()
        plt.scatter(df['timestamp'], df['m_battpos'], label='m_battpos')
        plt.scatter(df['timestamp'], df['c_battpos'], label='c_battpos')
        plt.xlabel('Time (secs)')
        plt.ylabel('Battery Position')
        plt.legend()
        plt.title('Battery Position over Time')
        plt.show()
        if not approve("Test 4: m_battpos and c_battpos versus time. Approve?"):
            return

    # Test 5: Get the average value for m_roll
    avg_m_roll = df['m_roll'].mean()
    if not approve(f"Test 5: Average m_roll value = {avg_m_roll}. Approve?"):
        return
    
    # Test 6: Plot m_fin and c_fin over time
    if 'm_fin' in df.columns and 'c_fin' in df.columns:
        plt.figure()
        plt.scatter(df['timestamp'], df['m_fin'], label='m_fin')
        plt.scatter(df['timestamp'], df['c_fin'], label='c_fin')
        plt.xlabel('Time (secs)')
        plt.ylabel('Fin Position')
        plt.legend()
        plt.title('Fin Position over Time')
        plt.show()
        if not approve("Test 6: m_fin and c_fin versus time. Approve?"):
            return

    # Test 7: Plot m_heading and c_heading over time
    if 'm_heading' in df.columns and 'c_heading' in df.columns:
        plt.figure()
        plt.scatter(df['timestamp'], df['m_heading'], label='m_heading')
        plt.scatter(df['timestamp'], df['c_heading'], label='c_heading')
        plt.xlabel('Time (secs)')
        plt.ylabel('Heading')
        plt.legend()
        plt.title('Heading over Time')
        plt.show()
        if not approve("Test 7: m_heading and c_heading versus time. Approve?"):
            return

    # Test 8: Plot m_dist_to_wpt and x_hit_a_waypoint over time
    if 'm_dist_to_wpt' in df.columns and 'x_hit_a_waypoint' in df.columns:
        plt.figure()
        plt.scatter(df['timestamp'], df['m_dist_to_wpt'], label='m_dist_to_wpt')
        plt.scatter(df['timestamp'], df['x_hit_a_waypoint'], label='x_hit_a_waypoint')
        plt.xlabel('Time (secs)')
        plt.ylabel('Distance to Waypoint')
        plt.legend()
        plt.title('Distance to Waypoint over Time')
        plt.show()
        if not approve("Test 8: m_dist_to_wpt and x_hit_a_waypoint versus time. Approve?"):
            return

    # Test 9: Plot m_vacuum and m_depth over time
    if 'm_vacuum' in df.columns and 'm_depth' in df.columns:
        plt.figure()
        plt.scatter(df['timestamp'], df['m_vacuum'], label='m_vacuum')
        plt.scatter(df['timestamp'], df['m_depth'], label='m_depth')
        plt.xlabel('Time (secs)')
        plt.ylabel('Vacuum and Depth')
        plt.legend()
        plt.title('Vacuum and Depth over Time')
        plt.show()
        if not approve("Test 9: m_vacuum and m_depth versus time. Approve?"):
            return

    # Test 10: Plot m_battery and m_depth over time
    if 'm_battery' in df.columns and 'm_depth' in df.columns:
        plt.figure()
        plt.scatter(df['timestamp'], df['m_battery'], label='m_battery')
        plt.scatter(df['timestamp'], df['m_depth'], label='m_depth')
        plt.xlabel('Time (secs)')
        plt.ylabel('Battery and Depth')
        plt.legend()
        plt.title('Battery and Depth over Time')
        plt.show()
        if not approve("Test 10: m_battery and m_depth versus time. Approve?"):
            return

    # Test 11: Get the first and last value for m_coulomb_amphr_total and calculate the delta change
    if 'm_coulomb_amphr_total' in df.columns:
        first_value = df['m_coulomb_amphr_total'].dropna().iloc[0]
        last_value = df['m_coulomb_amphr_total'].dropna().iloc[-1]
        delta_change = last_value - first_value
        if not approve(f"Test 11: m_coulomb_amphr_total - First value: {first_value}, Last value: {last_value}, Delta change: {delta_change}. Approve?"):
            return

    # Test 12: Plot m_coulomb_amphr_total against time
    if 'm_coulomb_amphr_total' in df.columns:
        plt.figure()
        plt.scatter(df['timestamp'], df['m_coulomb_amphr_total'])
        plt.xlabel('Time (secs)')
        plt.ylabel('Coulomb Amp-Hour Total')
        plt.title('Coulomb Amp-Hour Total over Time')
        plt.show()
        if not approve("Test 12: m_coulomb_amphr_total versus time. Approve?"):
            return

    # Test 13: Plot any variables that have names following this format and wildcard: m_bms_*_current over time
    bms_columns = [col for col in df.columns if col.startswith('m_bms_') and col.endswith('_current')]
    if bms_columns:
        plt.figure()
        for col in bms_columns:
            plt.scatter(df['timestamp'], df[col], label=col)
        plt.xlabel('Time (secs)')
        plt.ylabel('Current')
        plt.legend()
        plt.title('BMS Currents over Time')
        plt.show()
        if not approve("Test 13: BMS currents versus time. Approve?"):
            return

    # Test 14: Get the first and last value for m_tot_num_inflections and calculate the delta change
    if 'm_tot_num_inflections' in df.columns:
        first_value = df['m_tot_num_inflections'].dropna().iloc[0]
        last_value = df['m_tot_num_inflections'].dropna().iloc[-1]
        delta_change = last_value - first_value
        if not approve(f"Test 14: m_tot_num_inflections - First value: {first_value}, Last value: {last_value}, Delta change: {delta_change}. Approve?"):
            return

    # Test 15: Plot m_tot_num_inflections against time
    if 'm_tot_num_inflections' in df.columns:
        plt.figure()
        plt.scatter(df['timestamp'], df['m_tot_num_inflections'])
        plt.xlabel('Time (secs)')
        plt.ylabel('Total Number of Inflections')
        plt.title('Total Number of Inflections over Time')
        plt.show()
        if not approve("Test 15: m_tot_num_inflections versus time. Approve?"):
            return

    # Test 16: Plot m_altitude, m_raw_altitude, and m_water_depth against time
    if 'm_altitude' in df.columns and 'm_raw_altitude' in df.columns and 'm_water_depth' in df.columns:
        plt.figure()
        plt.scatter(df['timestamp'], df['m_altitude'], label='m_altitude')
        plt.scatter(df['timestamp'], df['m_raw_altitude'], label='m_raw_altitude')
        plt.scatter(df['timestamp'], df['m_water_depth'], label='m_water_depth')
        plt.xlabel('Time (secs)')
        plt.ylabel('Altitude and Water Depth')
        plt.legend()
        plt.title('Altitude and Water Depth over Time')
        plt.show()
        if not approve("Test 16: m_altitude, m_raw_altitude, and m_water_depth versus time. Approve?"):
            return

    # Test 17: Plot m_depth, c_argos_on, and c_air_pump against time
    if 'm_depth' in df.columns and 'c_argos_on' in df.columns and 'c_air_pump' in df.columns:
        plt.figure()
        plt.scatter(df['timestamp'], df['m_depth'], label='m_depth')
        plt.scatter(df['timestamp'], df['c_argos_on'], label='c_argos_on')
        plt.scatter(df['timestamp'], df['c_air_pump'], label='c_air_pump')
        plt.xlabel('Time (secs)')
        plt.ylabel('Depth, Argos On, and Air Pump')
        plt.legend()
        plt.title('Depth, Argos On, and Air Pump over Time')
        plt.show()
        if not approve("Test 17: m_depth, c_argos_on, and c_air_pump versus time. Approve?"):
            return

    # Test 18: Display the average value for m_mission_avg_speed_climbing and m_mission_avg_speed_diving
    if 'm_mission_avg_speed_climbing' in df.columns and 'm_mission_avg_speed_diving' in df.columns:
        avg_speed_climbing = df['m_mission_avg_speed_climbing'].mean()
        avg_speed_diving = df['m_mission_avg_speed_diving'].mean()
        if not approve(f"Test 18: Average speed climbing: {avg_speed_climbing}, Average speed diving: {avg_speed_diving}. Approve?"):
            return

