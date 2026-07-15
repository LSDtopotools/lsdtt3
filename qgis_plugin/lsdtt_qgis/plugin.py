import os

from qgis.core import QgsApplication
from qgis.core import QgsSettings
from qgis.PyQt.QtGui import QIcon
from qgis.PyQt.QtWidgets import QAction

from .provider import LSDTTProvider
from .ui.settings_dialog import SettingsDialog


class LSDTTPlugin:
    def __init__(self, iface):
        self.iface = iface
        self.provider = None
        self.settings_action = None

    def initGui(self):
        self.provider = LSDTTProvider()
        QgsApplication.processingRegistry().addProvider(self.provider)

        icon = QIcon(os.path.join(os.path.dirname(__file__), "icon.png"))
        self.settings_action = QAction(icon, "LSDTT Backend Settings", self.iface.mainWindow())
        self.settings_action.triggered.connect(self.open_settings_dialog)
        self.iface.addPluginToMenu("lsdtt3", self.settings_action)

    def open_settings_dialog(self):
        dialog = SettingsDialog(QgsSettings(), self.iface.mainWindow())
        dialog.exec_()

    def unload(self):
        if self.settings_action is not None:
            self.iface.removePluginMenu("lsdtt3", self.settings_action)
            self.settings_action = None

        if self.provider is not None:
            QgsApplication.processingRegistry().removeProvider(self.provider)
            self.provider = None


def classFactory(iface):
    return LSDTTPlugin(iface)
