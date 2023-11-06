# =========================
# X - Imports
# =========================

import os
from datetime import datetime

# =========================
# X - Slocum Glider Files
# =========================

### FUNCTION:
def goto_l10(config, waypoints, directory, calculate_heading, DD_to_DDMM):
    
    """
    Generate the 'goto' file using the given configuration and save it to the specified location.
    
    Args:
    - config (dict): The GGS configuration dictionary.
    - directory (str): The path to the directory where the file should be saved.
    - calculate_heading (func): A function to calculate the heading.
    - DD_to_DDMM (func): A function to convert from decimal degrees to DDMM format.
    """
    
    utc_datetime = datetime.utcnow()

    goto = []

    goto.append("behavior_name=goto_list\n")
    goto.append("#==========================================================================================")
    goto.append(f"# Glider: {config['glider_name']} --- goto_l10.ma")
    goto.append(f"# Satisfying Radius: {config['satisfying_radius']}")
    goto.append(f"# Generated the Glider Guidance System (GGS) on {datetime.utcnow().strftime('%m-%d-%Y')} at {datetime.utcnow().strftime('%H:%M')} UTC")
    goto.append("#==========================================================================================\n")

    goto.append("<start:b_arg>")
    goto.append("\n<end:b_arg>\n")

    goto.append("<start:waypoints>")
    goto.append("#   LON        LAT        Heading | Waypoint")

    for i, wpt in enumerate(waypoints):
        lon_ddmm = DD_to_DDMM(wpt[1])
        lat_ddmm = DD_to_DDMM(wpt[0])

        if i == 0:
            goto.append(f"    {lon_ddmm}   {lat_ddmm}   # N/A | Waypoint {i} /// Start Point")
        elif i == len(waypoints) - 1:
            leg_heading = calculate_heading(waypoints[i-1], wpt)
            goto.append(f"    {lon_ddmm}   {lat_ddmm}   # {leg_heading:.0f} | Waypoint {i} /// End Point")
        else:
            leg_heading = calculate_heading(waypoints[i-1], wpt)
            goto.append(f"    {lon_ddmm}   {lat_ddmm}   # {leg_heading:.0f} | Waypoint {i}")

    goto.append("<end:waypoints>")
    goto_content = '\n'.join(goto)

    filename = f"{config['glider_name']}_goto_l10.ma"
    
    goto_file_path = os.path.join(directory, filename)
    
    try:
        with open(goto_file_path, 'w') as GotoFile:
            GotoFile.write(goto_content)
        print(f"{filename} saved at {directory}")
    except IOError as e:
        print(f"Error writing to file: {e}")