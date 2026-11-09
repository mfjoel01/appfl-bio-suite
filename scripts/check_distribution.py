"""Check the built wheel's metadata and essential installed resources before upload."""

from email.parser import BytesParser
from pathlib import Path
from zipfile import ZipFile

from packaging.requirements import Requirement


def main():
    wheels = list(Path("dist").glob("*.whl"))
    assert len(wheels) == 1, f"Expected one release wheel, found {wheels}"
    with ZipFile(wheels[0]) as wheel:
        names = wheel.namelist()
        metadata = BytesParser().parsebytes(
            wheel.read(next(name for name in names if name.endswith(".dist-info/METADATA")))
        )
        for value in metadata.get_all("Requires-Dist", []):
            assert Requirement(value).url is None, f"PyPI rejects direct dependencies: {value}"
        for suffix in [
            "resources/constraints.txt",
            "resources/partner/endpoint-setup.md",
            "resources/templates/validate_mapping.py",
            "resources/templates/user_config_template.yaml.j2",
            "experiments/fine_mapping/fedfm/UPSTREAM.json",
        ]:
            assert f"appfl_bio_suite/{suffix}" in names, f"Missing wheel resource: {suffix}"
        for directory, extension in [
            ("core/ga4gh/data", ".json"),
            ("core/watch_assets", ".js"),
            ("core/plot_assets", ".css"),
            ("experiments/gwas/configs", ".yaml"),
            ("experiments/fine_mapping/configs", ".yaml"),
        ]:
            assert any(
                n.startswith(f"appfl_bio_suite/{directory}/") and n.endswith(extension)
                for n in names
            ), f"Missing packaged assets: {directory}"
        print(f"Validated {wheels[0]} ({metadata['Name']} {metadata['Version']})")


if __name__ == "__main__":
    main()
