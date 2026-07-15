from qgis.PyQt.QtWidgets import (
    QDialog,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLineEdit,
    QPushButton,
    QVBoxLayout,
)

from ..backend.manager import SETTINGS_BACKEND_PATH


class SettingsDialog(QDialog):
    def __init__(self, settings, parent=None):
        super().__init__(parent)
        self.settings = settings
        self.setWindowTitle("LSDTT Backend Settings")
        self.resize(520, 120)

        self.backend_path = QLineEdit(self)
        self.backend_path.setText(self.settings.value(SETTINGS_BACKEND_PATH, "", type=str))

        browse = QPushButton("Browse", self)
        browse.clicked.connect(self._browse_backend_dir)

        path_row = QHBoxLayout()
        path_row.addWidget(self.backend_path)
        path_row.addWidget(browse)

        form = QFormLayout()
        form.addRow("Backend folder", path_row)

        save = QPushButton("Save", self)
        save.clicked.connect(self._save)

        root = QVBoxLayout()
        root.addLayout(form)
        root.addWidget(save)
        self.setLayout(root)

    def _browse_backend_dir(self):
        selected = QFileDialog.getExistingDirectory(self, "Select lsdtt backend folder")
        if selected:
            self.backend_path.setText(selected)

    def _save(self):
        self.settings.setValue(SETTINGS_BACKEND_PATH, self.backend_path.text().strip())
        self.accept()
