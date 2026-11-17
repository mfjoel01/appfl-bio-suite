"""The site layout: the map at the root, one documentation page per section under docs/.

Run against the real website/ pages, so a section added, renamed or linked there is
checked the day it lands: every page exists, every local link and #anchor resolves, each
page marks its own menu entry, and old addresses still lead somewhere.
"""

from __future__ import annotations

import importlib.util
import re
from urllib.parse import urlsplit

import pytest

from appfl_bio_suite.core.experiments import repo_root

SITE = repo_root() / "website"


@pytest.fixture(scope="module")
def built(tmp_path_factory):
    spec = importlib.util.spec_from_file_location(
        "build_docs", repo_root() / "scripts" / "build_docs.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    out = tmp_path_factory.mktemp("site")
    (out / "map").mkdir()
    (out / "map" / "map.html").write_text(
        '<!doctype html><html lang="en"><head><title>viewer</title></head><body></body></html>'
    )
    (out / "assets").mkdir()
    for asset in (SITE / "assets").iterdir():
        (out / "assets" / asset.name).write_bytes(asset.read_bytes())
    (out / "assets" / "suite-logo.png").write_bytes(b"")
    module.build(SITE, out)
    return module, out


def test_every_section_has_its_page(built):
    module, out = built
    for slug in [*module.PAGES, module.NETWORK_PAGE]:
        assert (out / "docs" / slug / "index.html").is_file(), slug
    assert (out / "index.html").is_file() and (out / "network.html").is_file()


def test_every_local_link_and_anchor_resolves(built):
    _, out = built
    missing = []
    for page in sorted((out / "docs").rglob("index.html")):
        text = page.read_text(encoding="utf-8")
        for link in re.findall(r'(?:href|src)="([^"]+)"', text):
            if re.match(r"(https?:|mailto:|data:|javascript:)", link):
                continue
            parts = urlsplit(link)
            target = (page.parent / parts.path).resolve() if parts.path else page
            if parts.path.endswith("/") or target.is_dir():
                target = target / "index.html"
            if not target.exists():
                missing.append(f"{page.relative_to(out)}: {link}")
            elif parts.fragment and target.suffix == ".html":
                if f'id="{parts.fragment}"' not in target.read_text(encoding="utf-8"):
                    missing.append(f"{page.relative_to(out)}: {link} (no anchor)")
    assert not missing, "\n".join(missing)


def test_each_page_marks_its_own_menu_entry(built):
    module, out = built
    for slug in module.PAGES:
        text = (out / "docs" / slug / "index.html").read_text(encoding="utf-8")
        nav = re.search(r'<nav class="site-nav".*?</nav>', text, re.S).group(0)
        current = re.findall(r'<a href="([^"]+)" aria-current="page"', nav)
        assert len(current) == 1 and current[0].endswith(f"docs/{slug + '/' if slug else ''}")


def test_the_root_is_the_map_and_old_anchors_redirect(built):
    module, out = built
    root = (out / "index.html").read_text(encoding="utf-8")
    assert 'data-bio-default-metadata="map/network.map.json"' in root
    assert '"quickstart": "docs/installation/"' in root
    assert '"experiments": "docs/experiments/"' in root
    assert "url=./" in (out / "network.html").read_text(encoding="utf-8")


def test_a_section_with_no_page_stops_the_build(built, tmp_path):
    module, _ = built
    site = tmp_path / "website"
    site.mkdir()
    page = (SITE / "index.html").read_text(encoding="utf-8")
    extra = '<section id="new-section"><div class="wrap">New</div></section>\n</main>'
    (site / "index.html").write_text(page.replace("</main>", extra, 1), encoding="utf-8")
    (site / "network.html").write_text((SITE / "network.html").read_text(encoding="utf-8"))
    with pytest.raises(SystemExit, match="new-section"):
        module.build(site, tmp_path / "out")
