import yaml

def load_settings(file_paths):
    """
    Loads and merges settings from one or more YAML configuration files.

    Input:
        file_paths (list of str): List of paths to YAML configuration files.

    Output:
        values (dict): Merged 'value' entries from all YAML files.
        descriptions (dict): Merged 'description' entries from all YAML files.
        levels (dict): Merged 'levels' entries from all YAML files.
        level_descriptions (dict): Merged 'level-description' entries from all YAML files.
    """
    # Initialize dictionaries
    values = {}
    descriptions = {}
    levels = {}
    level_descriptions = {}

    for file_path in file_paths:
        with open(file_path, 'r') as file:
            settings = yaml.safe_load(file)

        for key in settings:
            values[key] = settings[key]['value']

            if 'description' in settings[key]:
                descriptions[key] = settings[key]['description']

            if 'levels' in settings[key]:
                levels[key] = settings[key]['levels']

            if 'level-description' in settings[key]:
                level_descriptions[key] = settings[key]['level-description']

    return values, descriptions, levels, level_descriptions
