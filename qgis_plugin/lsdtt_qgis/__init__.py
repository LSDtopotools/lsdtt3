try:
    from .plugin import classFactory
except ImportError:
    # qgis is unavailable outside a QGIS install (e.g. when running the
    # unit tests); qgis-free submodules remain importable.
    pass
