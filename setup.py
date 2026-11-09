"""Include the partner-only resource tree in wheels without duplicating its sources."""

import shutil
from pathlib import Path

from setuptools import setup
from setuptools.command.build_py import build_py


class BuildPy(build_py):
    def run(self):
        super().run()
        root = Path(__file__).parent
        target = Path(self.build_lib) / "appfl_bio_suite" / "resources"
        for source, name in [("docs/partner", "partner"), ("partners/templates", "templates")]:
            shutil.copytree(root / source, target / name, dirs_exist_ok=True)
        shutil.copy2(root / "constraints.txt", target / "constraints.txt")


setup(cmdclass={"build_py": BuildPy})
