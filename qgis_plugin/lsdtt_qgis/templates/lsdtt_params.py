def _format_value(value):
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def write_parameter_file(file_path, parameter_map):
    ordered_keys = sorted(parameter_map.keys())
    with open(file_path, "w", encoding="utf-8") as dst:
        for key in ordered_keys:
            dst.write("{}: {}\n".format(key, _format_value(parameter_map[key])))
