def get_screen_settings(main_dir, project_dir, subject, task_name):
    """
    Find screen settings for a subject using task-<task>_events.json files.
    Falls back to dataset-level JSON if no subject-specific one exists.
    """
    from bids import BIDSLayout
    from pathlib import Path
    import json

    layout = BIDSLayout(f"{main_dir}/{project_dir}", validate=False, index_metadata=True)

    # Get all JSONs for this task
    events_files = layout.get(
        suffix='events',
        extension='.json',
        task=task_name,
        return_type='file'
    )

    found = None
    # Look for subject-specific JSONs
    for ev in events_files:
        entities = layout.parse_file_entities(ev)
        sub = entities.get('subject', None)
        if sub == subject:
            found = ev
            break

    # Fallback: use root-level file if none found
    if found is None:
        root_json = Path(f"{main_dir}/{project_dir}") / f"task-{task_name}_events.json"
        if root_json.exists():
            found = str(root_json)
            print(f"[INFO] Using root-level {root_json}")
        else:
            raise FileNotFoundError(
                f"No subject-specific or root-level task-{task_name}_events.json found."
            )

    # Load and extract
    with open(found, 'r') as f:
        meta = json.load(f)

    stim_pres = meta.get('StimulusPresentation', None)
    if stim_pres is None:
        raise KeyError(f"No 'StimulusPresentation' entry found in {found}")

    # Convert from meters → centimeters
    screen_size_cm = [round(float(s) * 100, 2) for s in stim_pres['ScreenSize']]
    screen_distance_cm = round(float(stim_pres['ScreenDistance']) * 100, 2)

    print(f"for subject {subject}:")
    print(f"[INFO] Using screen info from {found}")
    print(f"[INFO] Screen size: {screen_size_cm} cm, distance: {screen_distance_cm} cm")

    return screen_size_cm, screen_distance_cm
