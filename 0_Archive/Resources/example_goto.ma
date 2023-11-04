behavior_name=goto_list
#============================================================
# --- goto_l10.ma
# Generated for RU23 by the Glider Guidance System (GGS) at: Oct 16 2023 08:49:37
#============================================================

<start:b_arg>
    b_arg: num_legs_to_run(nodim)   -1    # Loop, Run all waypoints
    b_arg: start_when(enum)          0    # BAW_IMMEDIATELY
    b_arg: list_stop_when(enum)      7    # BAW_WHEN_WPT_DIST
    b_arg: num_waypoints(nodim)      5    # Number of waypoints in list
    b_arg: initial_wpt(enum)	     0    # 0 to n-1, -1 first after last, -2 closest
<end:b_arg>

<start:waypoints>
#     LON     LAT      °T         name
    -7411.587    3927.437   #  |  Start Point
    -7356.607    3917.583   #  |  Waypoint 1
    -7339.926    3907.744   #  |  Waypoint 2
    -7320.627    3855.899   #  |  Waypoint 3
    -7304.705    3846.046   #  |  Waypoint 4, Turn Point
    -7320.628    3855.900   #  |  Waypoint 5 (3)
    -7339.927    3907.745   #  |  Waypoint 6 (2)
    -7356.608    3917.584   #  |  Waypoint 7 (1)

<end:waypoints>
