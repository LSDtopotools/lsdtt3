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
            return self._normalise_backend_path(configured)
        return self._normalise_backend_path(os.environ.get("LSDTT_BACKEND_PATH", ""))

    def _normalise_backend_path(self, backend_path):
        if not backend_path:
            return ""

        path = os.path.abspath(backend_path)
        if not os.path.isdir(path):
            return path

        # Common user mistake: selecting <bundle>/bin instead of the bundle root.
        # If so, transparently promote to the parent root when it looks valid.
        base = os.path.basename(path).lower()
        parent = os.path.dirname(path)
        if base == "bin":
            parent_share_proj = os.path.join(parent, "share", "proj")
            parent_share_gdal = os.path.join(parent, "share", "gdal")
            if os.path.isdir(parent_share_proj) and os.path.isdir(parent_share_gdal):
                return parent

        return path

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

        # Prevent GDAL from trying to load QGIS-installed plugin drivers inside
        # the backend subprocess. Those plugin DLLs are built against QGIS' GDAL
        # and can emit noisy "specified procedure could not be found" messages
        # when resolved against backend runtime DLLs.
        env["GDAL_DRIVER_PATH"] = "disable"
        env.pop("GDAL_PYTHON_DRIVER_PATH", None)
        env.pop("OGR_DRIVER_PATH", None)

        return env
