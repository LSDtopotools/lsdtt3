import os
import platform

from qgis.core import QgsSettings


SETTINGS_BACKEND_PATH = "lsdtt_qgis/backend_path"


class BackendManager:
    def __init__(self, settings=None):
        self.settings = settings or QgsSettings()

    def backend_path(self):
        configured = self.settings.value(SETTINGS_BACKEND_PATH, "", type=str)
        if configured:
            return configured
        return os.environ.get("LSDTT_BACKEND_PATH", "")

    def executable_path(self, cli_name):
        backend = self.backend_path()
        if not backend:
            raise RuntimeError("Backend path is not configured. Set it in plugin settings.")

        executable_name = cli_name
        if platform.system().lower().startswith("win") and not executable_name.endswith(".exe"):
            executable_name = "{}.exe".format(executable_name)

        candidates = [
            os.path.join(backend, executable_name),
            os.path.join(backend, "bin", executable_name),
        ]
        for candidate in candidates:
            if os.path.exists(candidate):
                return candidate

        raise RuntimeError("Executable not found for '{}': {}".format(cli_name, candidates))

    def runtime_environment(self):
        backend = self.backend_path()
        if not backend:
            raise RuntimeError("Backend path is not configured. Set it in plugin settings.")

        env = dict(os.environ)
        backend_bin = os.path.join(backend, "bin")
        if os.path.isdir(backend_bin):
            env["PATH"] = backend_bin + os.pathsep + env.get("PATH", "")

        gdal_data = os.path.join(backend, "share", "gdal")
        if os.path.isdir(gdal_data):
            env["GDAL_DATA"] = gdal_data

        proj_data = os.path.join(backend, "share", "proj")
        if os.path.isdir(proj_data):
            env["PROJ_DATA"] = proj_data
            env["PROJ_LIB"] = proj_data

        return env
