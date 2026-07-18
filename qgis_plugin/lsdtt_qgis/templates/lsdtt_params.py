def _format_value(value):
    if isinstance(value, bool):
        return "true" if value else "false"
    return str(value)


def write_parameter_file(file_path, parameter_map, header_lines=None):
    """Write an lsdtt parameter file; header_lines become leading # comments."""
    ordered_keys = sorted(parameter_map.keys())
    with open(file_path, "w", encoding="utf-8") as dst:
        for line in header_lines or []:
            dst.write("# {}\n".format(line))
        if header_lines:
            dst.write("\n")
        for key in ordered_keys:
            dst.write("{}: {}\n".format(key, _format_value(parameter_map[key])))
