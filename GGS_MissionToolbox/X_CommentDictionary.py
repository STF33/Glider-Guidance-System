"""
Author: Salvatore Fricano
Repository: https://github.com/STF33/Glider-Guidance-System
"""

# =========================

### FUNCTION:
def comment_dictionary(config_dictionary):

    '''
    A repository of Slocum glider behavior argument comments for various behaviors and their dependencies.

    Args:
    config_dictionary (dict): A dictionary of configuration behaviors and behavior arguments.

    Returns:
    b_arg_comments (dict): A dictionary of the comments for each behaviors arguments.
    '''
    
    b_arg_comments = {}
    behaviors = config_dictionary.get('behaviors', {})

    for behavior, details in behaviors.items():
        b_arg_comments[behavior] = {}
        for sub_behavior, sub_details in details.items():
            b_arg_comments[behavior][sub_behavior] = {}
            b_args = sub_details.get('b_args', {})
            # print(f"b_args = {b_args}")
            for b_arg, b_arg_value in b_args.items():
                try:
                    b_arg_value = float(b_arg_value)
                except ValueError:
                    b_arg_value = b_arg_value

                #########
                # ABEND #
                #########
                if behavior == "abend":
                    if b_arg == "overdepth(m)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# MS_ABORT_OVERDEPTH: Overdepth threshold = {b_arg_value} meters OR F_MAX_WORKING_DEPTH"
                    elif b_arg == "overdepth_sample_time(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# MS_ABORT_OVERDEPTH: Overdepth sample time = {b_arg_value} seconds"
                    elif b_arg == "overtime(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# MS_ABORT_OVERTIME: Overtime threshold = {b_arg_value} seconds"
                    elif b_arg == "undervolts(volts)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# MS_ABORT_UNDERVOLTS = Unvervolts threshold = {b_arg_value} volts"
                    elif b_arg == "undervolts_sample_time(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# MS_ABORT_UNDERVOLTS = Undervolts sample time = {b_arg_value} seconds"
                    elif b_arg == "samedepth_for(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# MS_ABORT_SAMEDEPTH: Samedepth threshold = {b_arg_value} seconds"
                    elif b_arg == "samedepth_for_sample_time(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# MS_ABORT_SAMEDEPTH: Samedepth sample time = {b_arg_value} seconds"
                    elif b_arg == "stalled_for(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# MS_ABORT_STALLED: Stalled threshold = {b_arg_value} seconds"
                    elif b_arg == "stalled_for_sample_time(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# MS_ABORT_STALLED: Stalled sample time = {b_arg_value} seconds"
                    elif b_arg == "no_cop_tickle_for(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# MS_ABORT_NO_TICKLE: No COP tickle threshold = {b_arg_value} seconds"
                    elif b_arg == "no_cop_tickle_percent(%)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# MS_ABORT_NO_TICKLE: No COP tickle threshold = {b_arg_value} %"
                    elif b_arg == "no_comms_tickle_for(hours)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# MS_ABORT_NO_COMMS_TICKLE: No comms tickle threshold = {b_arg_value} hours"
                    elif b_arg == "max_wpt_distance(m)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# MS_ABORT_WPT_TOOFAR: Max waypoint distance = {b_arg_value} meters"
                    elif b_arg == "chk_sensor_reasonableness(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" MS_ABORT_UNREASONABLE_SETTINGS: Check sensor reasonableness = False"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" MS_ABORT_UNREASONABLE_SETTINGS: Check sensor reasonableness = True"
                    elif b_arg == "reqd_free_heap(bytes)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# MS_ABORT_NO_HEAP: Required free heap = {b_arg_value} bytes"
                    elif b_arg == "leakdetect_sample_time(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# MS_ABORT_LEAK: Leakdetect sample time = {b_arg_value} seconds"
                    elif b_arg == "vacuum_min(inHg)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# MS_ABORT_VACUUM: Vacuum min = {b_arg_value} inHg"
                    elif b_arg == "vacuum_max(inHg)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# MS_ABORT_VACUUM: Vacuum max = {b_arg_value} inHg"
                    elif b_arg == "vacuum_sample_time(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# MS_ABORT_VACUUM: Vacuum sample time = {b_arg_value} seconds"
                    elif b_arg == "oil_volume_sample_time(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# MS_ABORT_VACUUM: Oil volume sample time = {b_arg_value} seconds"
                    elif b_arg == "remaining_charge_min(%)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# MS_ABORT_CHARGE_MIN: Remaining charge min = {b_arg_value} %"
                    elif b_arg == "remaining_charge_sample_time(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# MS_ABORT_CHARGE_MIN: Remaining charge sample time = {b_arg_value} seconds"
                    elif b_arg == "invalid_gps(nodim)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# MS_ABORT_INVALID_GPS: Invalid GPS threshold = {b_arg_value}"
                    elif b_arg == "check_emergency_battery_active(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" MS_ABORT_EMERGENCY_BATTERY_ACTIVE: Check emergency battery active = False"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" MS_ABORT_EMERGENCY_BATTERY_ACTIVE: Check emergency battery active = True"
                    elif b_arg == "samedepth_for_surfacing(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# MS_ABORT_SURFACE_BLOCKED: Samedepth for surfacing threshold = {b_arg_value} seconds"

                ###########
                # SURFACE #
                ###########
                elif behavior == "surface":
                    if b_arg == "args_from_file(enum)":
                        b_arg_value = int(float(b_arg_value))
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# mafiles/surfac{b_arg_value}.ma"
                    elif b_arg == "start_when(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 0 = BAW_IMMEDIATELY"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 1 = BAW_STK_IDLE"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 2 = BAW_PITCH_IDLE"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 3 = BAW_HEADING_IDLE"
                        elif b_arg_value == 6:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 6 = BAW_WHEN_SECS"
                        elif b_arg_value == 7:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 7 = BAW_WHEN_WPT_DIST"
                        elif b_arg_value == 8:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 8 = BAW_WHEN_HIT_WAYPOINT"
                        elif b_arg_value == 9:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 9 = BAW_EVERY_SECS"
                        elif b_arg_value == 10:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 10 = BAW_EVERY_SECS_UPDWN_IDLE"
                        elif b_arg_value == 11:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 11 = BAW_SCI_SURFACE"
                        elif b_arg_value == 12:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 12 = BAW_NOCOMM_SECS"
                        elif b_arg_value == 13:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 13 = BAW_WHEN_UTC_TIME"
                        elif b_arg_value == 16:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 16 = BAW_NUM_INFLECTIONS"
                    elif b_arg == "when_secs(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        start_when_value = b_args.get("start_when(enum)", None)
                        if start_when_value is not None:
                            start_when_value = int(start_when_value)
                        if start_when_value in [6, 9, 12]:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Surface every {b_arg_value} seconds"
                        elif start_when_value == 13:
                            if b_arg_value == -1:
                                b_arg_comments[behavior][sub_behavior][b_arg] += f" Surface triggers anytime after the when_utc_timestamp"
                            elif b_arg_value == 0:
                                b_arg_comments[behavior][sub_behavior][b_arg] += f" Surface once at the when_utc_timestamp"
                            elif b_arg_value > 0:
                                b_arg_comments[behavior][sub_behavior][b_arg] += f" Surface on a repetition period after the when_utc_timestamp"
                    elif b_arg == "when_wpt_dist(m)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Surface when waypoint distance is within {b_arg_value} meters"
                    elif b_arg == "when_num_inflections(nodim)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        start_when_value = b_args.get("start_when(enum)", None)
                        if start_when_value is not None:
                            start_when_value = int(start_when_value)
                        if start_when_value == 16:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f"# Surface after {b_arg_value} inflections"
                    elif b_arg == "end_action(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" end_action 0 = quit"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" end_action 1 = wait for a Ctrl-C to quit"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" end_action 2 = resume"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" end_action 3 = drift until END_WPT_DIST"
                        elif b_arg_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" end_action 4 = wait for a Ctrl-C once"
                        elif b_arg_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" end_action 5 = wait for a Ctrl-C quit on timeout"
                    elif b_arg == "report_all(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Report just GPS"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Report all sensors once"
                    elif b_arg == "gps_wait_time(s)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Wait {b_arg_value} seconds for a GPS fix"
                    elif b_arg == "keystroke_wait_time(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Wait {b_arg_value} seconds for a keystroke"
                    elif b_arg == "end_wpt_dist(m)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        end_action_value = b_args.get("end_action(enum)", None)
                        if end_action_value is not None:
                            end_action_value = int(end_action_value)
                        if end_action_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f"# Stop when M_DIST_TO_WPT is < {b_arg_value} meters"
                    elif b_arg == "c_use_bpump(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Autoballast"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Reserved, Use C_BPUMP_VALUE(X) around 0"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Absolute, use C_BPUMP_VALUE(X) as the total difference between dives and climbs"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Absolute, use C_BPUMP_VALUE(X) during the surface interval"
                        elif b_arg_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Speed relative to water current"
                    elif b_arg == "c_bpump_value(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        c_use_bpump_value = b_args.get("c_use_bpump(enum)", None)
                        if c_use_bpump_value is not None:
                            c_use_bpump_value = int(c_use_bpump_value)
                        if c_use_bpump_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" C_BPUMP_VALUE(X) is ignored"
                        elif c_use_bpump_value in [1, 2, 3]:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" climb with {b_arg_value} cc (clipped to max legal)"
                    elif b_arg == "c_use_pitch(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 1 = battpos"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 2 = pitch (set once from curve)"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 3 = servo on pitch"
                        elif b_arg_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 4 = fluid pumped absolute"
                    elif b_arg == "c_pitch_value(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        c_use_pitch_value = b_args.get("c_use_pitch(enum)", None)
                        if c_use_pitch_value is not None:
                            c_use_pitch_value = int(c_use_pitch_value)
                        if c_use_pitch_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" battpos = {b_arg_value} inches"
                        elif c_use_pitch_value in [2, 3]:
                            degrees = b_arg_value * (180 / 3.14159)
                            degrees = round(degrees, 2)
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" pitch = {b_arg_value} radians = {degrees} degrees"
                        elif c_use_pitch_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" fluid pumped absolute = {b_arg_value} cc"
                    elif b_arg == "c_stop_when_air_pump(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Terminate climb once air pump has been inflated, for use with thruster only"
                    elif b_arg == "c_use_thruster(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Thruster not in use"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command input voltage (as % of glider voltage)"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command input voltage (as % of max thruster voltage)"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command depth change rate"
                        elif b_arg_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command input wattage (between 1 and 9)"
                    elif b_arg == "c_thruster_value(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = "#"
                        c_use_thruster_value = b_args.get("c_use_thruster(enum)", None)
                        if c_use_thruster_value is not None:
                            c_use_thruster_value = int(c_use_thruster_value)
                        if c_use_thruster_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Thruster is not in use"
                        elif c_use_thruster_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded input voltage = {b_arg_value} % of glider voltage"
                        elif c_use_thruster_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded input voltage = {b_arg_value} % of max thruster voltage"
                        elif c_use_thruster_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded depth change rate = {b_arg_value} m/s"
                        elif c_use_thruster_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded input wattage = {b_arg_value} watts"
                    elif b_arg == "printout_cycle_time(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Printout dialog every {b_arg_value} seconds"
                    elif b_arg == "gps_postfix_wait_time(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Wait {b_arg_value} seconds afterinitial GPS fix to turn on iridium"
                    elif b_arg == "force_iridium_use(nodim)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Force iridium even if FreeWave if present"
                    elif b_arg == "min_time_between_gps_fixes(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Iridium will be hung up every {b_arg_value} seconds to get a GPS fix"
                    elif b_arg == "sensor_input_wait_time(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Wait {b_arg_value} seconds for input sensors at the surface"
                    elif b_arg == "when_utc_timestamp(dtime)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        yymmddhhmm = str(b_arg_value)
                        year = yymmddhhmm[0:2]
                        month = yymmddhhmm[2:4]
                        if month == "01":
                            month = "JAN"
                        elif month == "02":
                            month = "FEB"
                        elif month == "03":
                            month = "MAR"
                        elif month == "04":
                            month = "APR"
                        elif month == "05":
                            month = "MAY"
                        elif month == "06":
                            month = "JUN"
                        elif month == "07":
                            month = "JUL"
                        elif month == "08":
                            month = "AUG"
                        elif month == "09":
                            month = "SEP"
                        elif month == "10":
                            month = "OCT"
                        elif month == "11":
                            month = "NOV"
                        elif month == "12":
                            month = "DEC"
                        day = yymmddhhmm[4:6]
                        hour = yymmddhhmm[6:8]
                        minute = yymmddhhmm[8:10]
                        b_arg_comments[behavior][sub_behavior][b_arg] += f" Surface on {day} {month} {year} at {hour}:{minute} UTC, with reference to U_MISSION_YEAR_BASE"
                    elif b_arg == "when_utc_on_surface(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Adjust WHEN_UTC_TIMESTAMP to get the glider to the surface by the WHEN_UTC_TIMESTAMP"
                    elif b_arg == "strobe_on(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Strobe off"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Strobe on"
                    elif b_arg == "thruster_burst(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Thruster burst off"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Thruster burst on"

                ############
                # GOTO_WPT #
                ############
                elif behavior == "goto_wpt":
                    if b_arg == "start_when(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 0 = BAW_IMMEDIATELY"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 1 = BAW_STK_IDLE"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 2 = BAW_PITCH_IDLE"
                        elif b_arg_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 4 = BAW_UPDWN_IDLE"
                    elif b_arg == "stop_when(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 0 = complete"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 1 = BAW_STK_IDLE"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 2 = BAW_PITCH_IDLE"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 3 = BAW_HEADING_IDLE"
                        elif b_arg_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 4 = BAW_UPDWN_IDLE"
                        elif b_arg_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 5 = BAW_NEVER"
                        elif b_arg_value == 6:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 6 = BAW_WHEN_SECS"
                        elif b_arg_value == 7:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 7 = BAW_WHEN_WPT_DIST"
                        elif b_arg_value == 8:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 8 = BAW_WHEN_HIT_WAYPOINT"
                        elif b_arg_value == 9:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 9 = BAW_EVERY_SECS"
                        elif b_arg_value == 10:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 10 = BAW_EVERY_SECS_UPDWN_IDLE"
                        elif b_arg_value == 11:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 11 = BAW_SCI_SURFACE"
                        elif b_arg_value == 12:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 12 = BAW_NOCOMM_SECS"
                        elif b_arg_value == 13:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 13 = BAW_WHEN_UTC_TIME"
                        elif b_arg_value == 14:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 14 = BAW_HOVER_ACTIVE"
                        elif b_arg_value == 15:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 15 = BAW_WHEN_WPT_EXCEEDS_DIST"
                        elif b_arg_value == 16:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 16 = BAW_NUM_INFLECTIONS"
                        elif b_arg_value == 17:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 17 = BAW_WHEN_PRIMARY_WPT_DIST"
                        elif b_arg_value == 18:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 18 = BAW_WHEN_PRIMARY_WPT_EXCEEDS_DIST"
                    elif b_arg == "when_wpt_dist(m)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        stop_when_value = b_args.get("stop_when_value(enum)", None)
                        if stop_when_value is not None:
                            stop_when_value = int(stop_when_value)
                        if stop_when_value == 7:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Stop when waypoint distance is within {b_arg_value} meters"
                    elif b_arg == "wpt_units(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" wpt_units 0 = LMC"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" wpt_units 1 = UTM"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" wpt_units 2 = LAT/LON"
                    elif b_arg == "wpt_x(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# The waypoint x-coordinate is {b_arg_value}"
                    elif b_arg == "wpt_y(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# The waypoint y-coordinate is {b_arg_value}"
                    elif b_arg == "utm_zd(byte)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# The UTM zone digit is {b_arg_value}"
                    elif b_arg == "utm_zc(byte)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# The UTM zone character is {b_arg_value}"
                    elif b_arg == "end_action(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" end_action 0 = quit"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" end_action 1 = wait for a Ctrl-C to quit"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" end_action 2 = resume"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" end_action 3 = drift until END_WPT_DIST"
                        elif b_arg_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" end_action 4 = wait for a Ctrl-C once"
                        elif b_arg_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" end_action 5 = wait for a Ctrl-C quit on timeout"

                #############
                # GOTO_LIST #
                #############
                elif behavior == "goto_list":
                    if b_arg == "args_from_file(enum)":
                        b_arg_value = int(float(b_arg_value))
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# mafiles/goto_l{b_arg_value}.ma"
                    elif b_arg == "start_when(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 0 = BAW_IMMEDIATELY"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 1 = BAW_STK_IDLE"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 2 = BAW_PITCH_IDLE"
                        elif b_arg_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 4 = BAW_UPDWN_IDLE"
                    elif b_arg == "num_waypoints(nodim)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Number of waypoints = {b_arg_value}"
                    elif b_arg == "num_legs_to_run(nodim)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == -2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Traverse waypoint list once"
                        elif b_arg_value == -1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Loop forever"
                        elif b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" ILLEGAL VALUE"
                        elif b_arg_value >= 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Sequence through {b_arg_value} waypoints"
                    elif b_arg == "initial_wpt(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == -2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Start at the closest waypoint"
                        elif b_arg_value == -1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Start at the next waypoint after the last one achieved"
                        elif b_arg_value >= 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Start at waypoint {b_arg_value}"
                    elif b_arg == "list_stop_when(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 1 = BAW_STK_IDLE"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 2 = BAW_PITCH_IDLE"
                        elif b_arg_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 5 = BAW_NEVER"
                        elif b_arg_value == 7:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 7 = BAW_WHEN_WPT_DIST"
                        elif b_arg_value == 15:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 15 = BAW_WHEN_WPT_EXCEEDS_DIST"
                    elif b_arg == "list_when_wpt_dist(m)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        list_stop_when_value = b_args.get("list_stop_when(enum)", None)
                        if list_stop_when_value is not None:
                            list_stop_when_value = int(list_stop_when_value)
                        if list_stop_when_value == 7:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Stop when waypoint distance is within {b_arg_value} meters"
                    elif b_arg == "end_action(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" end_action 0 = quit"
                        elif b_arg_value == 6:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" end_action 6 = wait for a Ctrl-f to re-read mafiles"
                    elif b_arg == "primary_wpt(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# The final waypoint in the waypoint list is the primary waypoint"
                    elif b_arg == "primary_stop_when(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = "#"
                        if b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 1 = BAW_STK_IDLE"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 2 = BAW_PITCH_IDLE"
                        elif b_arg_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 5 = BAW_NEVER"
                        elif b_arg_value == 7:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 7 = BAW_WHEN_WPT_DIST"
                        elif b_arg_value == 15:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 15 = BAW_WHEN_WPT_EXCEEDS_DIST"
                    elif b_arg == "primary_when_wpt_dist(m)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = "#"
                        list_stop_when_value = b_args.get("list_stop_when(enum)", None)
                        if list_stop_when_value is not None:
                            list_stop_when_value = int(list_stop_when_value)
                        if list_stop_when_value in [7, 15]:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Stop when primary waypoint distance is within {b_arg_value} meters"
                        
                ######
                # YO #
                ######
                elif behavior == "yo":
                    if b_arg == "args_from_file(enum)":
                        b_arg_value = int(float(b_arg_value))
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# mafiles/yo{b_arg_value}.ma"
                    elif b_arg == "start_when(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 0 = BAW_IMMEDIATELY"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 1 = BAW_STK_IDLE"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 2 = BAW_PITCH_IDLE"
                        elif b_arg_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 4 = BAW_UPDWN_IDLE"
                        elif b_arg_value == 7:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 7 = BAW_WHEN_WPT_DIST"
                        elif b_arg_value == 8:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 8 = BAW_WHEN_HIT_WAYPOINT"
                        elif b_arg_value == 15:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 15 = BAW_WHEN_WPT_EXCEEDS_DIST"
                        elif b_arg_value == 17:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 17 = BAW_WHEN_PRIMARY_WPT_DIST"
                        elif b_arg_value == 18:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 18 = BAW_WHEN_PRIMARY_WPT_EXCEEDS_DIST"
                    elif b_arg == "start_diving(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Climb first"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Dive first"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Continue present"
                    elif b_arg == "num_half_cycles_to_do(nodim)":
                        num_yos_to_do = b_arg_value/2
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Number of half cycles to do = {b_arg_value} = {num_yos_to_do} yos"
                    elif b_arg == "d_target_depth(m)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if 3 < b_arg_value < 1000:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Dive to {b_arg_value} meters" 
                        else:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" WARNING: OUT OF LIMITS"
                    elif b_arg == "d_target_altitude(m)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == -1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Altimeter disabled on dive"
                        elif 0 < b_arg_value < 100:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Target altitude = {b_arg_value} meters"
                        else:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" WARNING: OUT OF LIMITS"
                    elif b_arg == "d_use_bpump(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Autoballast"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Reserved, Use D_BPUMP_VALUE(X) around 0"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Absolute, use D_BPUMP_VALUE(X) as the total difference between dives and climbs"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Absolute, use D_BPUMP_VALUE(X) during the surface interval"
                        elif b_arg_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Speed relative to water current"
                    elif b_arg == "d_bpump_value(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        d_use_bpump_value = b_args.get("d_use_bpump(enum)", None)
                        if d_use_bpump_value is not None:
                            d_use_bpump_value = int(d_use_bpump_value)
                        if d_use_bpump_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" D_BPUMP_VALUE(X) is ignored"
                        elif d_use_bpump_value in [1, 2, 3]:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Dive with {b_arg_value} cc (clipped to max legal)"
                    elif b_arg == "d_use_pitch(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 1 = battpos"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 2 = pitch (set once from curve)"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 3 = servo on pitch"
                        elif b_arg_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 4 = fluid pumped absolute"
                    elif b_arg == "d_pitch_value(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        d_use_pitch_value = b_args.get("d_use_pitch(enum)", None)
                        if d_use_pitch_value is not None:
                            d_use_pitch_value = int(d_use_pitch_value)
                        if d_use_pitch_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" battpos = {b_arg_value} inches"
                        elif d_use_pitch_value in [2, 3]:
                            degrees = b_arg_value * (180 / 3.14159)
                            degrees = round(degrees, 2)
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" pitch = {b_arg_value} radians = {degrees} degrees"
                        elif d_use_pitch_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" fluid pumped absolute = {b_arg_value} cc"
                    elif b_arg == "d_stop_when_hover_for(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Stop when hovering for {b_arg_value} seconds"
                    elif b_arg == "d_stop_when_stalled_for(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Stop when stalled for {b_arg_value} seconds"
                    elif b_arg == "d_speed_min(m/s)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Minimum speed = {b_arg_value} m/s"
                    elif b_arg == "d_speed_max(m/s)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Maximum speed = {b_arg_value} m/s"
                    elif b_arg == "d_use_thruster(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Thruster not in use"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command input voltage (as % of glider voltage)"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command input voltage (as % of max thruster voltage)"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command depth change rate"
                        elif b_arg_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command input wattage (between 1 and 9)"
                    elif b_arg == "d_thruster_value(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        d_use_thruster_value = b_args.get("d_use_thruster(enum)", None)
                        if d_use_thruster_value is not None:
                            d_use_thruster_value = int(d_use_thruster_value)
                        if d_use_thruster_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Thruster is not in use"
                        elif d_use_thruster_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded input voltage = {b_arg_value} % of glider voltage"
                        elif d_use_thruster_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded input voltage = {b_arg_value} % of max thruster voltage"
                        elif d_use_thruster_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded depth change rate = {b_arg_value} m/s"
                        elif d_use_thruster_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded input wattage = {b_arg_value} watts"
                    elif b_arg == "d_depth_rate_method(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Depth rate method 0 = raw m_depth_rate"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Depth rate method 1 = m_depth_rate_subsample"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Depth rate method 2 = ILLEGAL"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Depth rate method 3 = running average (m_depth_rate_avg_final)"
                    elif b_arg == "d_wait_for_pitch(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Wait for pitch"
                    elif b_arg == "d_wait_for_ballast(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Wait for {b_arg_value} seconds for ballast"
                    elif b_arg == "d_delta_bpump_speed(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Delta ballast pump speed = {b_arg_value}"
                    elif b_arg == "d_delta_bpump_ballast(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Delta ballast pump ballast = {b_arg_value}"
                    elif b_arg == "d_time_ratio(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Time ratio = {b_arg_value}"
                    elif b_arg == "d_use_sc_model(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Use SC model off"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Use SC model on"
                    elif b_arg == "d_max_thermal_charge_time(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Maximum thermal charge time = {b_arg_value} seconds"
                    elif b_arg == "d_max_pumping_charge_time(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Maximum pumping charge time = {b_arg_value} seconds"
                    elif b_arg == "d_thr_reqd_pres_mul(nodim)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Thruster required pressure multiplier = {b_arg_value}"
                    elif b_arg == "c_target_depth(m)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if 3 < b_arg_value < 1000:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Climb to {b_arg_value} meters"
                        else:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" WARNING: OUT OF LIMITS"
                    elif b_arg == "c_target_altitude(m)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == -1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Altimeter disabled on climb"
                        else:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Target altitude = {b_arg_value} meters"
                    elif b_arg == "c_use_bpump(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Autoballast"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Reserved, Use C_BPUMP_VALUE(X) around 0"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Absolute, use C_BPUMP_VALUE(X) as the total difference between dives and climbs"
                        elif b_arg_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Speed relative to water current"
                    elif b_arg == "c_bpump_value(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        c_use_bpump_value = b_args.get("c_use_bpump(enum)", None)
                        if c_use_bpump_value is not None:
                            c_use_bpump_value = int(c_use_bpump_value)
                        if c_use_bpump_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" C_BPUMP_VALUE(X) is ignored"
                        elif c_use_bpump_value in [1, 2]:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Climb with {b_arg_value} cc (clipped to max legal)"
                        elif c_use_bpump_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Climb with speed relative to water current"
                    elif b_arg == "c_use_pitch(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 1 = battpos"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 2 = pitch (set once from curve)"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 3 = servo on pitch"
                        elif b_arg_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 4 = fluid pumped absolute"
                    elif b_arg == "c_pitch_value(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        c_use_pitch_value = b_args.get("c_use_pitch(enum)", None)
                        if c_use_pitch_value is not None:
                            c_use_pitch_value = int(c_use_pitch_value)
                        if c_use_pitch_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" battpos = {b_arg_value} inches"
                        elif c_use_pitch_value in [2, 3]:
                            degrees = b_arg_value * (180 / 3.14159)
                            degrees = round(degrees, 2)
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" pitch = {b_arg_value} radians = {degrees} degrees"
                        elif c_use_pitch_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" fluid pumped absolute = {b_arg_value} cc"
                    elif b_arg == "c_stop_when_hover_for(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Stop when hovering for {b_arg_value} seconds"
                    elif b_arg == "c_stop_when_stalled_for(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Stop when stalled for {b_arg_value} seconds"
                    elif b_arg == "c_speed_min(m/s)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Minimum speed = {b_arg_value} m/s"
                    elif b_arg == "c_speed_max(m/s)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Maximum speed = {b_arg_value} m/s"
                    elif b_arg == "c_use_thruster(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Thruster not in use"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command input voltage (as % of glider voltage)"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command input voltage (as % of max thruster voltage)"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command depth change rate"
                        elif b_arg_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command input wattage (between 1 and 9)"
                    elif b_arg == "c_thruster_value(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        c_use_thruster_value = b_args.get("c_use_thruster(enum)", None)
                        if c_use_thruster_value is not None:
                            c_use_thruster_value = int(c_use_thruster_value)
                        if c_use_thruster_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Thruster is not in use"
                        elif c_use_thruster_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded input voltage = {b_arg_value} % of glider voltage"
                        elif c_use_thruster_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded input voltage = {b_arg_value} % of max thruster voltage"
                        elif c_use_thruster_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded depth change rate = {b_arg_value} m/s"
                        elif c_use_thruster_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded input wattage = {b_arg_value} watts"
                    elif b_arg == "end_action(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" end_action 0 = quit"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" end_action 2 = resume"
                    elif b_arg == "stop_when(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 5 = BAW_NEVER"
                        elif b_arg_value == 6:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 6 = BAW_WHEN_SECS"
                        elif b_arg_value == 7:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 7 = BAW_WHEN_WPT_DIST"
                        elif b_arg_value == 8:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 8 = BAW_WHEN_HIT_WAYPOINT"
                        elif b_arg_value == 15:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 15 = BAW_WHEN_WPT_EXCEEDS_DIST"
                        elif b_arg_value == 17:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 17 = BAW_WHEN_PRIMARY_WPT_DIST"
                        elif b_arg_value == 18:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 18 = BAW_WHEN_PRIMARY_WPT_EXCEEDS_DIST"
                    elif b_arg == "when_secs(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        start_when_value = b_args.get("start_when(enum)", None)
                        if start_when_value is not None:
                            start_when_value = int(start_when_value)
                        if start_when_value in [6, 9, 10]:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Stop after {b_arg_value} seconds"
                    elif b_arg == "when_wpt_dist(m)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Stop when waypoint distance is within {b_arg_value} meters"

                ###################
                # PREPARE_TO_DIVE #
                ###################
                elif behavior == "prepare_to_dive":
                    if b_arg == "args_from_file(enum)":
                        b_arg_value = int(float(b_arg_value))
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# mafiles/prepar{b_arg_value}.ma"
                    elif b_arg == "start_when(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 0 = BAW_IMMEDIATELY"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 1 = BAW_STK_IDLE"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 2 = BAW_PITCH_IDLE"
                    elif b_arg == "wait_time(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Wait {b_arg_value} seconds for a GPS fix"
                    elif b_arg == "max_thermal_charge_time(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Maximum thermal charge time = {b_arg_value} seconds"
                    elif b_arg == "max_pumping_charge_time(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Maximum pumping charge time = {b_arg_value} seconds"

                ##############
                # SENSORS_IN #
                ##############
                elif behavior == "sensors_in":
                    if b_arg == "c_att_time(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = "#"
                        if b_arg_value < 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Off"
                        elif b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Sample as fast as possible"
                        elif b_arg_value > 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Sample every {b_arg_value} seconds"
                    elif b_arg == "c_pressure_time(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value < 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Off"
                        elif b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Sample as fast as possible"
                        elif b_arg_value > 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Sample every {b_arg_value} seconds"
                    elif b_arg == "c_alt_time(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value < 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Off"
                        elif b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Sample as fast as possible"
                        elif b_arg_value > 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Sample every {b_arg_value} seconds"
                    elif b_arg == "u_battery_time(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value < 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Off"
                        elif b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Sample as fast as possible"
                        elif b_arg_value > 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Sample every {b_arg_value} seconds"
                    elif b_arg == "u_vacuum_time(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value < 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Off"
                        elif b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Sample as fast as possible"
                        elif b_arg_value > 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Sample every {b_arg_value} seconds"
                    elif b_arg == "c_leakdetect_time(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value < 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Off"
                        elif b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Sample as fast as possible"
                        elif b_arg_value > 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Sample every {b_arg_value} seconds"
                    elif b_arg == "c_gps_on(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value < 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Off"
                        elif b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Sample as fast as possible"
                        elif b_arg_value > 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Sample every {b_arg_value} seconds"
                    elif b_arg == "c_science_all_on(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value < 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Off"
                        elif b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Sample as fast as possible"
                        elif b_arg_value > 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Sample every {b_arg_value} seconds"
                    elif b_arg == "c_profile_on(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value < 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Off"
                        elif b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Sample as fast as possible"
                        elif b_arg_value > 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Sample every {b_arg_value} seconds"

                ###############
                # SET_HEADING #
                ###############
                elif behavior == "set_heading":
                    if b_arg == "args_from_file(enum)":
                        b_arg_value = int(float(b_arg_value))
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# mafiles/set_he{b_arg_value}.ma"
                    elif b_arg == "use_heading(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = "#"
                        if b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_heading 2 = HM_HEADING"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_heading 2 = HM_ROLL"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_heading 3 = HM_BATTROLL"
                        elif b_arg_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_heading 4 = HM_FIN"
                        elif b_arg_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_heading 5 = HM_CURRENT"
                    elif b_arg == "heading_value(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        use_heading_value = b_args.get("use_heading(bool)", None)
                        if use_heading_value is not None:
                            use_heading_value = int(use_heading_value)
                        if use_heading_value == 1:
                            degrees = b_arg_value * (180 / 3.14159)
                            degrees = round(degrees, 2)
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Heading = {b_arg_value} radians = {degrees} degrees"
                        elif use_heading_value == 2:
                            degrees = b_arg_value * (180 / 3.14159)
                            degrees = round(degrees, 2)
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Roll = {b_arg_value} radians = {degrees} degrees"
                        elif use_heading_value == 3:
                            degrees = b_arg_value * (180 / 3.14159)
                            degrees = round(degrees, 2)
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Battroll = {b_arg_value} radians = {degrees} degrees"
                        elif use_heading_value == 4:
                            degrees = b_arg_value * (180 / 3.14159)
                            degrees = round(degrees, 2)
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Fin = {b_arg_value} radians = {degrees} degrees"
                        elif use_heading_value == 5:
                            degrees = b_arg_value * (180 / 3.14159)
                            degrees = round(degrees, 2)
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Heading = {b_arg_value} radians offset from the current = {degrees} degrees offset from the current"
                    elif b_arg == "start_when(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 0 = BAW_IMMEDIATELY"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 1 = BAW_STK_IDLE"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 2 = BAW_PITCH_IDLE"
                        elif b_arg_value == 6:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 6 = BAW_WHEN_SECS"
                        elif b_arg_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 5 = BAW_NEVER"
                        elif b_arg_value == 7:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 7 = BAW_WHEN_WPT_DIST"
                        elif b_arg_value == 8:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 8 = BAW_WHEN_HIT_WAYPOINT"
                        elif b_arg_value == 15:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 15 = BAW_WHEN_WPT_EXCEEDS_DIST"
                    elif b_arg == "stop_when(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 1 = BAW_STK_IDLE"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 2 = BAW_PITCH_IDLE"
                        elif b_arg_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 5 = BAW_NEVER"
                        elif b_arg_value == 6:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 6 = BAW_WHEN_SECS"
                        elif b_arg_value == 7:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 7 = BAW_WHEN_WPT_DIST"
                        elif b_arg_value == 8:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 8 = BAW_WHEN_HIT_WAYPOINT"
                        elif b_arg_value == 15:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" stop_when 15 = BAW_WHEN_WPT_EXCEEDS_DIST"
                    elif b_arg == "when_secs(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        start_when_value = b_args.get("start_when(enum)", None)
                        if start_when_value is not None:
                            start_when_value = int(start_when_value)
                        if start_when_value in [6, 9, 12]:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f"  Start after {b_arg_value} seconds"
                    elif b_arg == "end_action(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" end_action 0 = quit"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" end_action 2 = resume"

                ###########
                # DIVE_TO #
                ###########
                elif behavior == "dive_to":
                    if b_arg == "target_depth(m)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Dive to {b_arg_value} meters"
                    elif b_arg == "target_altitude(m)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Target altitude = {b_arg_value} meters"
                    elif b_arg == "use_bpump(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Autoballast"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Reserved, Use BPUMP_VALUE(X) around 0"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Absolute, use BPUMP_VALUE(X) as the total difference between dives and climbs"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Absolute, use BPUMP_VALUE(X) during the surface interval"
                        elif b_arg_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Speed relative to water current"
                    elif b_arg == "bpump_value(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        use_bpump_value = b_args.get("use_bpump(enum)", None)
                        if use_bpump_value is not None:
                            use_bpump_value = int(use_bpump_value)
                        if use_bpump_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" BPUMP_VALUE(X) is ignored"
                        elif use_bpump_value in [1, 2, 3]:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Dive with {b_arg_value} cc (clipped to max legal)"
                        elif use_bpump_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Dive with speed relative to water current"
                    elif b_arg == "use_pitch(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 1 = battpos"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 2 = pitch (set once from curve)"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 3 = servo on pitch"
                        elif b_arg_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 4 = fluid pumped absolute"
                    elif b_arg == "pitch_value(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = "#"
                        use_pitch_value = b_args.get("use_pitch(enum)", None)
                        if use_pitch_value is not None:
                            use_pitch_value = int(use_pitch_value)
                        if use_pitch_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" battpos = {b_arg_value} inches"
                        elif use_pitch_value in [2, 3]:
                            degrees = b_arg_value * (180 / 3.14159)
                            degrees = round(degrees, 2)
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" pitch = {b_arg_value} radians = {degrees} degrees"
                        elif use_pitch_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" fluid pumped absolute = {b_arg_value} cc"
                    elif b_arg == "start_when(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 0 = BAW_IMMEDIATELY"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 1 = BAW_STK_IDLE"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 2 = BAW_PITCH_IDLE"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 3 = BAW_HEADING_IDLE"
                        elif b_arg_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 4 = BAW_UPDWN_IDLE"
                        elif b_arg_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 5 = BAW_NEVER"
                        elif b_arg_value == 6:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 6 = BAW_WHEN_SECS"
                        elif b_arg_value == 7:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 7 = BAW_WHEN_WPT_DIST"
                        elif b_arg_value == 8:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 8 = BAW_WHEN_HIT_WAYPOINT"
                        elif b_arg_value == 9:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 9 = BAW_EVERY_SECS"
                        elif b_arg_value == 10:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 10 = BAW_EVERY_SECS_UPDWN_IDLE"
                        elif b_arg_value == 11:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 11 = BAW_SCI_SURFACE"
                        elif b_arg_value == 12:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 12 = BAW_NOCOMM_SECS"
                        elif b_arg_value == 13:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 13 = BAW_WHEN_UTC_TIME"
                        elif b_arg_value == 14:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 14 = BAW_HOVER_ACTIVE"
                        elif b_arg_value == 15:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 15 = BAW_WHEN_WPT_EXCEEDS_DIST"
                        elif b_arg_value == 16:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 16 = BAW_NUM_INFLECTIONS"
                        elif b_arg_value == 17:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 17 = BAW_WHEN_PRIMARY_WPT_DIST"
                        elif b_arg_value == 18:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 18 = BAW_WHEN_PRIMARY_WPT_EXCEEDS_DIST"
                    elif b_arg == "stop_when_hover_for(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Stop when hovering for {b_arg_value} seconds"
                    elif b_arg == "stop_when_stalled_for(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Stop when stalled for {b_arg_value} seconds"
                    elif b_arg == "stop_when_air_pump(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Stop when air pump is {b_arg_value}"
                    elif b_arg == "initial_inflection(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        use_bpump_value = b_args.get("use_bpump(enum)", None)
                        if use_bpump_value is not None:
                            use_bpump_value = int(use_bpump_value)
                        if use_bpump_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Start with an initil inflection"
                        else:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Ignored for USE_BPUMP != 0"
                    elif b_arg == "speed_min(m/s)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Minimum speed = {b_arg_value} m/s"
                    elif b_arg == "speed_max(m/s)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Maximum speed = {b_arg_value} m/s"
                    elif b_arg == "use_thruster(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Thruster not in use"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command input voltage (as % of glider voltage)"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command input voltage (as % of max thruster voltage)"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command depth change rate"
                        elif b_arg_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command input wattage (between 1 and 9)"
                    elif b_arg == "thruster_value(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        use_thruster_value = b_args.get("use_thruster(enum)", None)
                        if use_thruster_value is not None:
                            use_thruster_value = int(use_thruster_value)
                        if use_thruster_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Thruster is not in use"
                        elif use_thruster_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded input voltage = {b_arg_value} % of glider voltage"
                        elif use_thruster_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded input voltage = {b_arg_value} % of max thruster voltage"
                        elif use_thruster_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded depth change rate = {b_arg_value} m/s"
                        elif use_thruster_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded input wattage = {b_arg_value} watts"
                    elif b_arg == "depth_rate_method(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" depth_rate_method 0 = raw m_depth_rate"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" depth_rate_method 1 = m_depth_rate_subsample"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" ILLEGAL"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" depth_rate_method 3 = running average (m_depth_rate_avg_final)"
                    elif b_arg == "wait_for_pitch(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Do not wait for pitch/battpos to settle before enabling speed control"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Wait for pitch/battpos to settle before enabling speed control"
                    elif b_arg == "wait_for_ballast(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Do not wait for ballast to settle before enabling speed control"
                        elif b_arg_value > 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Wait for ballast to settle before enabling speed control"
                    elif b_arg == "delta_bpump_speed(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Amount of ballast to add to bpump for achieving desired speed"
                    elif b_arg == "delta_bpump_ballast(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Amount of ballast to add to bpump for converging on ballast"
                    elif b_arg == "time_ratio(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Ratio nof climb/dive times that must be maintained for speed control"
                    elif b_arg == "use_sc_model(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Do not use the speed control model"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Use the speed control model"
                    elif b_arg == "max_thermal_charge_time(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Maximum thermal charge time = {b_arg_value} seconds"
                    elif b_arg == "max_pumping_charge_time(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Maximum pumping charge time = {b_arg_value} seconds"
                    elif b_arg == "thr_reqd_pres_mul(nodim)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Engine pressure must be {b_arg_value} times the ocean pressure at depth before the dive is started"

                ############
                # CLIMB_TO #
                ############
                elif behavior == "climb_to":
                    if b_arg == "target_depth(m)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Climb to {b_arg_value} meters"
                    elif b_arg == "target_altitude(m)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Target altitude = {b_arg_value} meters"
                    elif b_arg == "use_bpump(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Autoballast"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Reserved, Use BPUMP_VALUE(X) around 0"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Absolute, use BPUMP_VALUE(X) as the total difference between dives and climbs"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Absolute, use BPUMP_VALUE(X) during the surface interval"
                        elif b_arg_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Speed relative to water current"
                    elif b_arg == "bpump_value(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        use_bpump_value = b_args.get("use_bpump(enum)", None)
                        if use_bpump_value is not None:
                            use_bpump_value = int(use_bpump_value)
                        if use_bpump_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" BPUMP_VALUE(X) is ignored"
                        elif use_bpump_value in [1, 2, 3]:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Climb with {b_arg_value} cc (clipped to max legal)"
                        elif use_bpump_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Climb with speed relative to water current"
                    elif b_arg == "use_pitch(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 1 = battpos"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 2 = pitch (set once from curve)"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 3 = servo on pitch"
                        elif b_arg_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 5 = fluid pumped absolute"
                    elif b_arg == "pitch_value(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        use_pitch_value = b_args.get("use_pitch(enum)", None)
                        if use_pitch_value is not None:
                            use_pitch_value = int(use_pitch_value)
                        if use_pitch_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" battpos = {b_arg_value} inches"
                        elif use_pitch_value in [2, 3]:
                            degrees = b_arg_value * (180 / 3.14159)
                            degrees = round(degrees, 2)
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" pitch = {b_arg_value} radians = {degrees} degrees"
                        elif use_pitch_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" fluid pumped absolute = {b_arg_value} cc"
                    elif b_arg == "start_when(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 0 = BAW_IMMEDIATELY"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 1 = BAW_STK_IDLE"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 2 = BAW_PITCH_IDLE"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 3 = BAW_HEADING_IDLE"
                        elif b_arg_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 4 = BAW_UPDWN_IDLE"
                        elif b_arg_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 5 = BAW_NEVER"
                        elif b_arg_value == 6:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 6 = BAW_WHEN_SECS"
                        elif b_arg_value == 7:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 7 = BAW_WHEN_WPT_DIST"
                        elif b_arg_value == 8:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 8 = BAW_WHEN_HIT_WAYPOINT"
                        elif b_arg_value == 9:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 9 = BAW_EVERY_SECS"
                        elif b_arg_value == 10:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 10 = BAW_EVERY_SECS_UPDWN_IDLE"
                        elif b_arg_value == 11:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 11 = BAW_SCI_SURFACE"
                        elif b_arg_value == 12:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 12 = BAW_NOCOMM_SECS"
                        elif b_arg_value == 13:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 13 = BAW_WHEN_UTC_TIME"
                        elif b_arg_value == 14:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 14 = BAW_HOVER_ACTIVE"
                        elif b_arg_value == 15:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 15 = BAW_WHEN_WPT_EXCEEDS_DIST"
                        elif b_arg_value == 16:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 16 = BAW_NUM_INFLECTIONS"
                        elif b_arg_value == 17:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 17 = BAW_WHEN_PRIMARY_WPT_DIST"
                        elif b_arg_value == 18:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 18 = BAW_WHEN_PRIMARY_WPT_EXCEEDS_DIST"
                    elif b_arg == "stop_when_hover_for(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Stop when hovering for {b_arg_value} seconds"
                    elif b_arg == "stop_when_stalled_for(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Stop when stalled for {b_arg_value} seconds"
                    elif b_arg == "stop_when_air_pump(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Terminate the climb once the air pump has been inflated"
                    elif b_arg == "initial_inflection(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        use_bpump_value = b_args.get("use_bpump(enum)", None)
                        if use_bpump_value is not None:
                            use_bpump_value = int(use_bpump_value)
                        if use_bpump_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Start with an initil inflection"
                        else:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Ignored for USE_BPUMP != 0"
                    elif b_arg == "speed_min(m/s)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Minimum speed = {b_arg_value} m/s"
                    elif b_arg == "speed_max(m/s)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Maximum speed = {b_arg_value} m/s"
                    elif b_arg == "use_thruster(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Thruster not in use"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command input voltage (as % of glider voltage)"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command input voltage (as % of max thruster voltage)"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command depth change rate"
                        elif b_arg_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command input wattage (between 1 and 9)"
                    elif b_arg == "thruster_value(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        use_thruster_value = b_args.get("use_thruster(enum)", None)
                        if use_thruster_value is not None:
                            use_thruster_value = int(use_thruster_value)
                        if use_thruster_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Thruster is not in use"
                        elif use_thruster_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded input voltage = {b_arg_value} % of glider voltage"
                        elif use_thruster_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded input voltage = {b_arg_value} % of max thruster voltage"
                        elif use_thruster_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded depth change rate = {b_arg_value} m/s"
                        elif use_thruster_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded input wattage = {b_arg_value} watts"

                ##################
                # DRIFT_AT_DEPTH #
                ##################
                elif behavior == "drift_at_depth":
                    if b_arg == "args_from_file(enum)":
                        b_arg_value = int(float(b_arg_value))
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# mafiles/drift{b_arg_value}.ma"
                    elif b_arg == "start_when(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 0 = BAW_IMMEDIATELY"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 1 = BAW_STK_IDLE"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 2 = BAW_PITCH_IDLE"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 3 = BAW_HEADING_IDLE"
                        elif b_arg_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 4 = BAW_UPDWN_IDLE"
                        elif b_arg_value == 6:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 6 = BAW_WHEN_SECS"
                        elif b_arg_value == 7:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 7 = BAW_WHEN_WPT_DIST"
                        elif b_arg_value == 8:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 8 = BAW_WHEN_HIT_WAYPOINT"
                        elif b_arg_value == 9:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 9 = BAW_EVERY_SECS"
                        elif b_arg_value == 10:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 10 = BAW_EVERY_SECS_UPDWN_IDLE"
                        elif b_arg_value == 13:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" start_when 13 = BAW_WHEN_UTC_TIME"
                    elif b_arg == "when_secs(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        start_when_value = b_args.get("start_when(enum)", None)
                        if start_when_value is not None:
                            start_when_value = int(start_when_value)
                        if start_when_value in [6, 9, 10]:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f"  Start after {b_arg_value} seconds"
                    # CONFIRM THIS IS CORRECT
                    elif b_arg == "when_wpt_dist(m)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        start_when_value = b_args.get("start_when(enum)", None)
                        if start_when_value is not None:
                            start_when_value = int(start_when_value)
                        if start_when_value == 7:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f"  Start when waypoint is {b_arg_value} meters away"
                    # CONFIRM THIS IS CORRECT
                    elif b_arg == "when_utc_timestamp(dtime)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        start_when_value = b_args.get("start_when(enum)", None)
                        if start_when_value is not None:
                            start_when_value = int(start_when_value)
                        if start_when_value == 13:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f"  Start at {b_arg_value}"
                    elif b_arg == "end_action(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" end_action 0 = quit"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" end_action 2 = resume"
                    elif b_arg == "stop_when_hover_for(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Stop when hovering for {b_arg_value} seconds"
                    elif b_arg == "est_time_to_settle(s)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Estimated time to settle = {b_arg_value} seconds"
                    elif b_arg == "target_depth(m)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Target depth = {b_arg_value} meters"
                    elif b_arg == "target_altitude(m)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Target altitude = {b_arg_value} meters"
                    elif b_arg == "alt_time(s)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Altimeter measurement frequency = {b_arg_value} seconds"
                    elif b_arg == "target_deadband(m)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Deadband = {b_arg_value} meters around target depth"
                    elif b_arg == "start_dist_from_target(m)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Start drifting {b_arg_value} meters from target depth/altitude"
                    elif b_arg == "depth_ctrl(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" depth_ctrl 0 = Defauly mode for buoyancy drive"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" depth_ctrl 1 = Use pitch to control depth"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" depth_ctrl 2 = Pitch-based depth control"
                    elif b_arg == "bpump_delta_value(cc)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Increments by {b_arg_value} cc to adjust X_HOVER_BALLAST to obtain neutral buoyancy"
                    elif b_arg == "bpump_delay(s)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Minimum time between making buoyancy adjustments"
                    elif b_arg == "bpump_deadz_width(cc)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# For temporarily adjusting the buoyancy pump"
                    elif b_arg == "bpump_db_frac_dz(nodim)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Deadband during the drift behavior"
                    elif b_arg == "use_pitch(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 1 = battpos"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 2 = pitch (set once from curve)"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 3 = servo on pitch"
                        elif b_arg_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 4 = fluid pumped absolute"
                    elif b_arg == "pitch_value(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        use_pitch_value = b_args.get("use_pitch(enum)", None)
                        if use_pitch_value is not None:
                            use_pitch_value = int(use_pitch_value)
                        if use_pitch_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" battpos = {b_arg_value} inches"
                        elif use_pitch_value in [2, 3]:
                            degrees = b_arg_value * (180 / 3.14159)
                            degrees = round(degrees, 2)
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" pitch = {b_arg_value} radians = {degrees} degrees"
                        elif use_pitch_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" fluid pumped absolute = {b_arg_value} cc"
                    elif b_arg == "wait_for_pitch(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Do not wait for pitch/battpos to settle before enabling speed control"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Wait for pitch/battpos to settle before enabling speed control"
                    elif b_arg == "use_thruster(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Thruster not in use"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command input voltage (as % of glider voltage)"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command input voltage (as % of max thruster voltage)"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command depth change rate"
                        elif b_arg_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command input wattage (between 1 and 9)"
                    elif b_arg == "thruster_value(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        use_thruster_value = b_args.get("use_thruster(enum)", None)
                        if use_thruster_value is not None:
                            use_thruster_value = int(use_thruster_value)
                        if use_thruster_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Thruster is not in use"
                        elif use_thruster_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded input voltage = {b_arg_value} % of glider voltage"
                        elif use_thruster_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded input voltage = {b_arg_value} % of max thruster voltage"
                        elif use_thruster_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded depth change rate = {b_arg_value} m/s"
                        elif use_thruster_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded input wattage = {b_arg_value} watts"
                    elif b_arg == "enable_steering(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Disable steering while hovering"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Enable steering while hovering"
                    elif b_arg == "d_use_bpump(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Autoballast"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Reserved, Use BPUMP_VALUE(X) around 0"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Absolute, use BPUMP_VALUE(X) as the total difference between dives and climbs"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Absolute, use BPUMP_VALUE(X) during the surface interval"
                        elif b_arg_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Speed relative to water current"
                    elif b_arg == "d_bpump_value(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        d_use_bpump_value = b_args.get("d_use_bpump(enum)", None)
                        if use_thruster_value is not None:
                            use_thruster_value = int(use_thruster_value)
                        if d_use_bpump_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" BPUMP_VALUE(X) is ignored"
                        elif d_use_bpump_value in [1, 2, 3]:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Dive with {b_arg_value} cc (clipped to max legal)"
                        elif d_use_bpump_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Dive with speed relative to water current"
                    elif b_arg == "d_use_pitch(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 1 = battpos"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 2 = pitch (set once from curve)"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 3 = servo on pitch"
                        elif b_arg_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 5 = fluid pumped absolute"
                    elif b_arg == "d_pitch_value(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        d_use_pitch_value = b_args.get("d_use_pitch(enum)", None)
                        if d_use_pitch_value is not None:
                            d_use_pitch_value = int(d_use_pitch_value)
                        if d_use_pitch_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" battpos = {b_arg_value} inches"
                        elif d_use_pitch_value in [2, 3]:
                            degrees = b_arg_value * (180 / 3.14159)
                            degrees = round(degrees, 2)
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" pitch = {b_arg_value} radians = {degrees} degrees"
                        elif d_use_pitch_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" fluid pumped absolute = {b_arg_value} cc"
                    elif b_arg == "d_use_thruster(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Thruster not in use"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command input voltage (as % of glider voltage)"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command input voltage (as % of max thruster voltage)"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command depth change rate"
                        elif b_arg_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command input wattage (between 1 and 9)"
                    elif b_arg == "d_thruster_value(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        d_use_thruster_value = b_args.get("d_use_thruster(enum)", None)
                        if d_use_thruster_value is not None:
                            d_use_thruster_value = int(d_use_thruster_value)
                        if d_use_thruster_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Thruster is not in use"
                        elif d_use_thruster_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded input voltage = {b_arg_value} % of glider voltage"
                        elif d_use_thruster_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded input voltage = {b_arg_value} % of max thruster voltage"
                        elif d_use_thruster_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded depth change rate = {b_arg_value} m/s"
                        elif d_use_thruster_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded input wattage = {b_arg_value} watts"
                    elif b_arg == "c_use_bpump(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Autoballast"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Reserved, Use BPUMP_VALUE(X) around 0"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Absolute, use BPUMP_VALUE(X) as the total difference between dives and climbs"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Absolute, use BPUMP_VALUE(X) during the surface interval"
                        elif b_arg_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Speed relative to water current"
                    elif b_arg == "c_bpump_value(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        c_use_bpump_value = b_args.get("c_use_bpump(enum)", None)
                        if c_use_bpump_value is not None:
                            c_use_bpump_value = int(c_use_bpump_value)
                        if c_use_bpump_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" BPUMP_VALUE(X) is ignored"
                        elif c_use_bpump_value in [1, 2, 3]:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Climb with {b_arg_value} cc (clipped to max legal)"
                        elif c_use_bpump_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Climb with speed relative to water current"
                    elif b_arg == "c_use_pitch(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 1 = battpos"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 2 = pitch (set once from curve)"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 3 = servo on pitch"
                        elif b_arg_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" use_pitch 5 = fluid pumped absolute"
                    elif b_arg == "c_pitch_value(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        c_use_pitch_value = b_args.get("c_use_pitch(enum)", None)
                        if c_use_pitch_value is not None:
                            c_use_pitch_value = int(c_use_pitch_value)
                        if c_use_pitch_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" battpos = {b_arg_value} inches"
                        elif c_use_pitch_value in [2, 3]:
                            degrees = b_arg_value * (180 / 3.14159)
                            degrees = round(degrees, 2)
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" pitch = {b_arg_value} radians = {degrees} degrees"
                        elif c_use_pitch_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" fluid pumped absolute = {b_arg_value} cc"
                    elif b_arg == "c_use_thruster(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Thruster not in use"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command input voltage (as % of glider voltage)"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command input voltage (as % of max thruster voltage)"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command depth change rate"
                        elif b_arg_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command input wattage (between 1 and 9)"
                    elif b_arg == "c_thruster_value(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        c_use_thruster_value = b_args.get("c_use_thruster(enum)", None)
                        if c_use_thruster_value is not None:
                            c_use_thruster_value = int(c_use_thruster_value)
                        if c_use_thruster_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Thruster is not in use"
                        elif c_use_thruster_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded input voltage = {b_arg_value} % of glider voltage"
                        elif c_use_thruster_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded input voltage = {b_arg_value} % of max thruster voltage"
                        elif c_use_thruster_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded depth change rate = {b_arg_value} m/s"
                        elif c_use_thruster_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Commanded input wattage = {b_arg_value} watts"
                    elif b_arg == "depth_pitch_limit(rad)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        degrees = b_arg_value * (180 / 3.14159)
                        degrees = round(degrees, 2)
                        b_arg_comments[behavior][sub_behavior][b_arg] += f" Limit pitch response to {b_arg_value} radians = {degrees} degrees"
                    elif b_arg == "depth_gain_scale(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Determines whether to use X_HOVER_DEPTH_P_GAIN = m * speed + b"
                    elif b_arg == "depth_p_gain(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Proportional gain = {b_arg_value}"
                    elif b_arg == "depth_i_gain(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Integral gain = {b_arg_value}"
                    elif b_arg == "depth_d_gain(X)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Derivative gain = {b_arg_value}"
                    elif b_arg == "depth_pitch_deadband(m/s)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Restricts buoyancy adjustment until depth rate is < {b_arg_value} m/s"
                    elif b_arg == "depth_pitch_max_time(s)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Maximum time allowed at maximum pitch = {b_arg_value} seconds"
                    elif b_arg == "pressure_median(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Enables or disables pressure median"
                    elif b_arg == "battpos_db(nodim)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# How much to scale from deadband (f_battpos_db_frac_dz)"

                ##########
                # SAMPLE #
                ##########
                elif behavior == "sample":
                    if b_arg == "args_from_file(enum)":
                        b_arg_value = int(float(b_arg_value))
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# mafiles/sample{b_arg_value}.ma"
                    elif b_arg == "sensor_type(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 0 = C_SCIENCE_ALL_ON"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 1 = C_PROFILE_ON"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 2 = C_HS2_ON !!REMOVED!!"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 3 = C_BB2F_ON"
                        elif b_arg_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 4 = C_BB2C_ON"
                        elif b_arg_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 5 = C_BB2LSS_ON"
                        elif b_arg_value == 6:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 6 = C_SAM_ON"
                        elif b_arg_value == 7:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 7 = C_WHPAR_ON !!REMOVED!!"
                        elif b_arg_value == 8:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 8 = C_WHGPBM_ON !!REMOVED!!"
                        elif b_arg_value == 9:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 9 = C_MOTEOPD_ON"
                        elif b_arg_value == 10:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 10 = C_BBFL2S_ON"
                        elif b_arg_value == 11:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 11 = C_FL3SLO_ON"
                        elif b_arg_value == 12:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 12 = C_BB3SLO_ON"
                        elif b_arg_value == 13:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 13 = C_OXY3835_ON"
                        elif b_arg_value == 14:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 14 = C_WHFCTD_ON"
                        elif b_arg_value == 15:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 15 = C_BAM_ON"
                        elif b_arg_value == 16:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 16 = C_OCR504R_ON"
                        elif b_arg_value == 17:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 17 = C_OCR504I_ON"
                        elif b_arg_value == 18:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 18 = C_BADD_ON"
                        elif b_arg_value == 19:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 19 = C_FLNTU_ON"
                        elif b_arg_value == 20:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 20 = C_FL3SLOV2_ON"
                        elif b_arg_value == 21:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 21 = C_BB3SLOV2_ON"
                        elif b_arg_value == 22:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 22 = C_OCR507R_ON"
                        elif b_arg_value == 23:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 23 = C_OCR507I_ON"
                        elif b_arg_value == 24:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 24 = C_BB3SLOV3_ON"
                        elif b_arg_value == 25:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 25 = C_BB2FLS_ON"
                        elif b_arg_value == 26:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 26 = C_BB2FLSV2_ON"
                        elif b_arg_value == 27:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 27 = C_OXY3835_WPHASE_ON"
                        elif b_arg_value == 28:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 28 = C_AUVB_ON"
                        elif b_arg_value == 29:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 29 = C_BB2FV2_ON"
                        elif b_arg_value == 30:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 30 = C_TARR_ON"
                        elif b_arg_value == 31:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 31 = C_BBFL2SV2_ON"
                        elif b_arg_value == 32:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 32 = C_GLBPS_ON"
                        elif b_arg_value == 33:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 33 = C_SSCSD_ON"
                        elif b_arg_value == 34:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 34 = C_BB2FLSV3_ON"
                        elif b_arg_value == 35:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 35 = C_FIRE_ON"
                        elif b_arg_value == 36:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 36 = C_OHF_ON !!REMOVED!!"
                        elif b_arg_value == 37:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 37 = C_BB2FLSV4_ON"
                        elif b_arg_value == 38:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 38 = C_BB2FLSV5_ON"
                        elif b_arg_value == 39:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 39 = C_LOGGER_ON"
                        elif b_arg_value == 40:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 40 = C_BBAM_ON"
                        elif b_arg_value == 41:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 41 = C_UMODEM_ON"
                        elif b_arg_value == 42:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 42 = C_RINKOII_ON"
                        elif b_arg_value == 43:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 43 = C_DVL_ON"
                        elif b_arg_value == 44:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 44 = C_BB2FLSV6_ON"
                        elif b_arg_value == 45:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 45 = C_FLBBRH_ON"
                        elif b_arg_value == 46:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 46 = C_FLUR_ON"
                        elif b_arg_value == 47:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 47 = C_BB2FLSV7_ON"
                        elif b_arg_value == 48:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 48 = C_FLBBCD_ON"
                        elif b_arg_value == 49:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 49 = C_DMON_ON"
                        elif b_arg_value == 50:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 50 = C_C3SFL_ON"
                        elif b_arg_value == 51:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 51 = C_SUNA_ON"
                        elif b_arg_value == 52:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 52 = C_SATPAR_ON"
                        elif b_arg_value == 53:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 53 = C_VSF_ON"
                        elif b_arg_value == 54:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 54 = C_OXY4_ON"
                        elif b_arg_value == 55:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 55 = C_GAMMA_RAD5_ON !!REMOVED!!"
                        elif b_arg_value == 56:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 56 = C_BSIPAR_ON"
                        elif b_arg_value == 57:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 57 = C_FLBB_ON"
                        elif b_arg_value == 58:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 58 = C_VR2C_ON"
                        elif b_arg_value == 59:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 59 = C_CTD41CP2_ON"
                        elif b_arg_value == 60:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 60 = C_ECHOSNDR853_ON"
                        elif b_arg_value == 61:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 61 = C_FLRH_ON"
                        elif b_arg_value == 62:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 62 = C_BB2FLSV8_ON"
                        elif b_arg_value == 63:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 63 = C_UVILUXPAH_ON"
                        elif b_arg_value == 64:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 64 = C_AD2CP_ON"
                        elif b_arg_value == 65:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 65 = C_MINIPROCO2_ON"
                        elif b_arg_value == 66:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 66 = C_PCO2_ON"
                        elif b_arg_value == 67:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 67 = C_SEAOWL_ON"
                        elif b_arg_value == 68:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 68 = C_AZFP_ON"
                        elif b_arg_value == 69:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 69 = C_UBAT_ON"
                        elif b_arg_value == 70:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 70 = C_LISST_ON"
                        elif b_arg_value == 71:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 71 = C_LMS_ON"
                        elif b_arg_value == 72:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 72 = C_SVS603_ON"
                        elif b_arg_value == 73:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 73 = C_MICRORIDER_ON"
                        elif b_arg_value == 74:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 74 = C_BB2FLSV9_ON"
                        elif b_arg_value == 75:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 75 = C_SBE41N_PH_ON"
                        elif b_arg_value == 76:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 76 = C_FL2URRH_ON"
                        elif b_arg_value == 77:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 77 = C_FLBBBBV1_ON"
                        elif b_arg_value == 78:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 78 = C_FLBBBBV2_ON"
                        elif b_arg_value == 79:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 79 = C_OBSVR_ON"
                        elif b_arg_value == 80:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 80 = C_FL2PECDOM_ON"
                        elif b_arg_value == 81:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 81 = C_WETLABSA_ON"
                        elif b_arg_value == 82:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 82 = C_WETLABSB_ON"
                        elif b_arg_value == 83:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 83 = C_WETLABSC_ON"
                        elif b_arg_value == 84:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 84 = C_ECHODROID_ON"
                        elif b_arg_value == 85:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 85 = C_TAU_ON"
                        elif b_arg_value == 86:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 86 = C_RBRODO_ON"
                        elif b_arg_value == 87:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 87 = C_SOLOCAM_ON"
                        elif b_arg_value == 88:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 88 = C_AMAR_ON"
                        elif b_arg_value == 89:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 89 = C_VRO_ON"
                        elif b_arg_value == 90:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" sensor_type 90 = C_EK80_ON"
                    elif b_arg == "state_to_sample(enum)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" state_to_sample 0 = none"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" state_to_sample 1 = diving"
                        elif b_arg_value == 2:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" state_to_sample 2 = hovering"
                        elif b_arg_value == 3:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" state_to_sample 3 = diving / hovering"
                        elif b_arg_value == 4:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" state_to_sample 4 = climbing"
                        elif b_arg_value == 5:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" state_to_sample 5 = diving / climbing"
                        elif b_arg_value == 6:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" state_to_sample 6 = hovering / climbing"
                        elif b_arg_value == 7:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" state_to_sample 7 = diving / hovering / climbing"
                        elif b_arg_value == 8:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" state_to_sample 8 = on surface"
                        elif b_arg_value == 9:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" state_to_sample 9 = diving / on surface"
                        elif b_arg_value == 10:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" state_to_sample 10 = hovering / on surface"
                        elif b_arg_value == 11:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" state_to_sample 11 = diving / hovering / on surface"
                        elif b_arg_value == 12:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" state_to_sample 12 = climbing / on surface"
                        elif b_arg_value == 13:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" state_to_sample 13 = diving / climbing / on surface"
                        elif b_arg_value == 14:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" state_to_sample 14 = hovering / climbing / on surface"
                        elif b_arg_value == 15:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" state_to_sample 15 = diving / hovering / climbing / on surface"
                    elif b_arg == "sample_time_after_state_change(s)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Time to wait after state change before sampling = {b_arg_value} seconds"
                    elif b_arg == "intersample_time(s)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Sample as fast as possible"
                        elif b_arg_value > 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Sample every {b_arg_value} second"
                    elif b_arg == "nth_yo_to_sample(nodim)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# After the first sample, sample every {b_arg_value} yo(s)"
                    elif b_arg == "intersample_depth(m)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value <= 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Disabled"
                        elif b_arg_value > 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Dyamically adjust sample rate to sample every {b_arg_value} meters"
                    elif b_arg == "min_depth(m)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Minimum depth to sample = {b_arg_value} meters"
                    elif b_arg == "max_depth(m)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"# Maximum depth to sample = {b_arg_value} meters"

                ############
                # NOP_CMDS #
                ############
                elif behavior == "nop_cmds":
                    if b_arg == "nop_pitch(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Disabled"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command pitch to _IGNORE to keep stack busy"
                    elif b_arg == "nop_bpump(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Disabled"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command bpump to _IGNORE to keep stack busy"
                    elif b_arg == "nop_heading(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Disabled"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command heading to _IGNORE to keep stack busy"
                    elif b_arg == "nop_threng(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Disabled"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command threng to _IGNORE to keep stack busy"
                    elif b_arg == "secs_to_run(sec)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Run forever"
                        elif b_arg_value > 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Run for {b_arg_value} seconds"
                    elif b_arg == "nop_air_pump(bool)":
                        b_arg_comments[behavior][sub_behavior][b_arg] = f"#"
                        if b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Disabled"
                        elif b_arg_value == -1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Do nothing"
                        elif b_arg_value == 0:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command air pump OFF"
                        elif b_arg_value == 1:
                            b_arg_comments[behavior][sub_behavior][b_arg] += f" Command air pump ON"
    
    return b_arg_comments