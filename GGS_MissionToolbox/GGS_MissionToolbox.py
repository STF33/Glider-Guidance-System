import os
import json
import re
from X_CommentDictionary import comment_dictionary

def create_mission(config_name="unknown"):

    '''
    Create a set of mission files for Slocum gliders.
    
    Args:
    - config_name (str): Name of the configuration file to import.
    
    Returns:
    - None
    '''

    current_directory = os.path.dirname(os.path.abspath(__file__))
    config_path = os.path.join(current_directory, 'config', f"{config_name}.json")
    mission_name = config_name

    output_directory = os.path.join(os.path.expanduser('~'), 'Downloads', mission_name)
    os.makedirs(output_directory, exist_ok=True)

    def load_config(file_path):

        '''
        Import a configuration file.
        
        Args:
        - file_path (str): Path to the specified configuration file.
        
        Returns:
        - config_dictionary (dict): The configuration dictionary with explicit keys.
        '''
        
        with open(file_path, 'r') as file:
            config = json.load(file)
        
        config_dictionary = {
            'behaviors': config.get('behaviors', {})
        }

        return config_dictionary

    def create_mi_file(config_dictionary, mission_name, output_directory, comments):

        '''
        Create the mission file.
        
        Args:
        - config_dictionary (dict): The configuration dictionary.
        - mission_name (str): The mission name.
        - output_directory (str): The path for the output mission file.
        - comments (dict): The comments dictionary.
        
        Returns:
        - None
        '''

        mission_content = f"##################################################\n"
        mission_content += f"# mission name = {mission_name}.mi\n"
        mission_content += f"##################################################\n\n"

        sensors = config_dictionary.get('sensors', {})
        for sensor, value in sensors.items():
            mission_content += f"sensor: {sensor} {value}\n"

        mission_content += "\n#########################\n\n"

        behaviors = config_dictionary.get('behaviors', {})

        if 'abend' in behaviors:
            mission_content += "behavior: abend\n"
            b_args = behaviors['abend'].get('b_args', {})
            for b_arg, b_arg_value in b_args.items():
                comment = comments.get('abend', {}).get(b_arg, "")
                mission_content += f"    b_arg: {b_arg} {b_arg_value} {comment}\n"
            mission_content += "\n#########################\n\n"

        for behavior, details in behaviors.items():
            if behavior == "abend":
                continue
            for sub_behavior, sub_details in details.items():
                b_args = sub_details.get('b_args', {})
                if "args_from_file(enum)" in b_args:
                    b_arg_value = b_args["args_from_file(enum)"]
                    comment = comments.get(behavior, {}).get(sub_behavior, {}).get("args_from_file(enum)", "")
                    mission_content += f"behavior: {behavior}\n"
                    mission_content += f"    b_arg: args_from_file(enum) {b_arg_value} # {comment}\n"
                    mission_content += "\n#########################\n\n"
                    create_ma_file(behavior, sub_behavior, sub_details, output_directory, comments)
                else:
                    mission_content += f"behavior: {behavior}\n"
                    for b_arg, b_arg_value in b_args.items():
                        comment = comments.get(behavior, {}).get(sub_behavior, {}).get(b_arg, "")
                        mission_content += f"    b_arg: {b_arg} {b_arg_value} {comment}\n"
                    mission_content += "\n#########################\n\n"

        with open(os.path.join(output_directory, f"{mission_name}.mi"), 'w') as file:
            file.write(mission_content)

    def create_ma_file(behavior, sub_behavior, sub_details, output_directory, comments):

        '''
        Create a single mission argument file.
        
        Args:
        - behavior (str): The behavior name.
        - sub_behavior (str): The sub-behavior name.
        - sub_details (dict): The sub-behavior details.
        - output_directory (str): The path for the output mission file.
        - comments (dict): The comments dictionary.
        
        Returns:
        - None
        '''

        file_name = os.path.join(output_directory, f"{sub_behavior}.ma")
        with open(file_name, 'w') as file:
            file.write(f"behavior_name={behavior}")
            file.write(f"\n\n<start:b_arg>\n\n")
            b_args = sub_details.get('b_args', {})
            for b_arg, b_arg_value in b_args.items():
                comment = comments.get(behavior, {}).get(sub_behavior, {}).get(b_arg, "")
                file.write(f"    b_arg: {b_arg} {b_arg_value} {comment}\n")
            file.write("\n<end:b_arg>\n")
            if behavior == "goto_list":
                waypoints = sub_details.get('waypoints', {})
                if waypoints:
                    file.write("\n<start:waypoints>\n\n")
                    for i, (waypoint, coords) in enumerate(waypoints.items()):
                        lon, lat = coords
                        file.write(f"{lon} {lat} # waypoint {i}\n")
                    file.write("\n<end:waypoints>\n")

    def create_data_lists(config_dictionary, output_directory):
        
        '''
        Create data list files based on the configuration.
        
        Args:
        - config_dictionary (dict): The configuration dictionary.
        - output_directory (str): The path for the output data list files.
        
        Returns:
        - None
        '''

        data_lists = ['sbdlist', 'mbdlist', 'tbdlist', 'nbdlist']
        
        for data_list in data_lists:
            if data_list in config_dictionary:
                with open(os.path.join(output_directory, f"{data_list}.dat"), 'w') as file:
                    sensors = config_dictionary[data_list]
                    for sensor, variables in sensors.items():
                        line = f"{sensor} {' '.join(variables)}\n"
                        file.write(line)

    config_dictionary = load_config(config_path)
    comments = comment_dictionary(config_dictionary)
    
    create_mi_file(config_dictionary, mission_name, output_directory, comments)
    create_data_lists(config_dictionary, output_directory)

    print(f"The mission files have been saved to: '{os.path.join(os.path.expanduser('~'), 'Downloads', config_name)}'")

if __name__ == "__main__":
    create_mission(config_name="spin_r")