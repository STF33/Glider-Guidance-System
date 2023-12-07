# =========================
# X - IMPORTS
# =========================

### /// ROUTE ANALYSIS ///
import os
from SUB_functions import calculate_distance, calculate_heading

# =========================
# ROUTE ANALYSIS
# =========================

### FUNCTION:
def route_analysis(config, waypoints):
    
    '''
    Analyze each leg of the route to compute metrics such as distance, time, and battery drain.
    
    Args:
    - config (dict): Glider Guidance System configuration
    - waypoints (list): list of waypoints

    Returns:
    - results (list): list of dictionaries containing the analysis results for each leg
    '''
    
    results = []
    
    for i in range(len(waypoints) - 1):
        leg_distance = calculate_distance(waypoints[i], waypoints[i + 1])
        leg_time_seconds = leg_distance / config['avg_velocity']
        leg_time_hours = leg_time_seconds / 3600
        battery_drain = config['battery_drain'] * (leg_time_seconds / 86400)
        
        leg_data = {
            "Leg Start": waypoints[i],
            "Leg End": waypoints[i + 1],
            "Leg Distance (meters)": leg_distance,
            "Time Taken (seconds)": leg_time_seconds,
            "Time Taken (hours)": leg_time_hours,
            "Expected Battery Drain (amp-hours)": battery_drain
        }
        
        results.append(leg_data)
    
    return results

### FUNCTION:
def route_analysis_output(config, directory, analysis_results):
    
    '''
    Print a summary of the route analysis and save the details in a text file.
    
    Args:
    - config (dict): Glider Guidance System configuration
    - directory (str): path to the directory containing the config files
    - analysis_results (list): list of dictionaries containing the analysis results for each leg

    Returns:
    - total_distance (float): total distance of the route
    - total_time_seconds (float): total time of the route in seconds
    - total_time_hours (float): total time of the route in hours
    - total_battery_drain (float): total battery drain of the route
    '''
    
    output_str = ""

    point_names = [f"Waypoint {i+1}" for i in range(len(config['waypoints']))]
    
    total_distance = 0
    total_time_seconds = 0
    total_time_hours = 0
    total_battery_drain = 0

    for idx, result in enumerate(analysis_results):
        leg_description = f"Leg {idx+1}: {point_names[idx]} {result['Leg Start']} to {point_names[idx+1]} {result['Leg End']}"
        leg_heading = calculate_heading(result['Leg Start'], result['Leg End'])
        
        output_str += (
            f"{leg_description}\n"
            f"Heading (degrees): {leg_heading}\n"
            f"Leg Distance (meters): {result['Leg Distance (meters)']}\n"
            f"Time Taken (seconds): {result['Time Taken (seconds)']}\n"
            f"Time Taken (hours): {result['Time Taken (hours)']}\n"
            f"Expected Battery Drain (amp-hours): {result['Expected Battery Drain (amp-hours)']}\n\n"
        )
        
        total_distance += result['Leg Distance (meters)']
        total_time_seconds += result['Time Taken (seconds)']
        total_time_hours += result['Time Taken (hours)']
        total_battery_drain += result['Expected Battery Drain (amp-hours)']

    output_str += (
        "Round Trip Details:\n"
        f"Total Distance (meters): {total_distance}\n"
        f"Total Time (seconds): {total_time_seconds}\n"
        f"Total Time (hours): {total_time_hours}\n"
        f"Total Expected Battery Drain (amp-hours): {total_battery_drain}\n"
    )
    
    print(output_str)

    analysis_filename = f"{config['glider_name']}_leg_analysis.txt"
    analysis_path = os.path.join(directory, analysis_filename)
    with open(analysis_path, 'w') as file:
        file.write(output_str)

    return total_distance, total_time_seconds, total_time_hours, total_battery_drain

# =========================
#
# analysis_results = route_analysis(config, waypoints)
#
# route_analysis_output(config, directory, analysis_results)
#  
# =========================