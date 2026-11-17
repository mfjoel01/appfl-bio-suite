#!/usr/bin/env python3
"""Lay out the project site: the map at its root, the documentation under docs/.

    python scripts/build_docs.py --site website --out _site

Run by scripts/build_site.sh after `watch export` has written the map to OUT/map. The
site's front door is the network map itself; the long project page in website/index.html
becomes one page per section, so each has its own address:

    /                       the map, full screen (a copy of map/map.html)
    /docs/                  About: the introduction, the map teaser and how to cite
    /docs/experiments/      ...and one page for each remaining section
    /docs/network-map/      website/network.html, how to read and export the map

The pages are still written in website/index.html and website/network.html, as one
document each; this script cuts them apart and rewrites their links, so the sources
stay easy to edit and nothing is maintained twice. A link to an anchor that now lives on
another page is pointed at that page. Old addresses keep working: network.html and the
old #section anchors on the root redirect to their new pages.
"""

from __future__ import annotations

import argparse
import html
import json
import re
from pathlib import Path

# Section id in website/index.html -> page under docs/. About collects the sections that
# introduce the project; every other section is a page of its own.
PAGES = {
    "": ("APPFL Bio Suite", ["hero", "about", "map-teaser", "citation"]),
    "experiments": ("Experiments", ["experiments"]),
    "design": ("Design", ["design"]),
    "standards": ("Standards", ["standards"]),
    "reproducibility": ("Reproducibility", ["reproducibility"]),
    "installation": ("Installation", ["quickstart"]),
}
NETWORK_PAGE = "network-map"
# The sections the long page's menu pointed at; a link to one now means its whole page.
PRIMARY = {"about", "experiments", "design", "standards", "reproducibility", "quickstart"}

_SECTION = re.compile(r"<section(?P<attrs>[^>]*)>.*?</section>", re.S)


def _section_key(attrs: str) -> str | None:
    found = re.search(r'id="([^"]+)"', attrs)
    if found:
        return found.group(1)
    return "hero" if 'class="hero"' in attrs else None


def split(page: str) -> tuple[str, str, str, dict[str, str]]:
    """(head, header, footer, {section key: html}) of a website/ page."""
    head = page[: page.index("<body>")]
    header = re.search(r"<header class=\"site-header\">.*?</header>", page, re.S).group(0)
    footer = re.search(r"<footer class=\"site-footer\">.*?</footer>", page, re.S).group(0)
    sections = {}
    for match in _SECTION.finditer(page):
        key = _section_key(match.group("attrs"))
        if key:
            sections[key] = match.group(0)
    return head, header, footer, sections


def anchor_pages(sections: dict[str, str]) -> dict[str, tuple[str, str]]:
    """Every id in the project page -> (docs page slug, section key that holds it)."""
    owner = {}
    for slug, (_, keys) in PAGES.items():
        for key in keys:
            for anchor in re.findall(r'\bid="([^"]+)"', sections.get(key, "")):
                owner[anchor] = (slug, key)
    return owner


def rewrite(text: str, root: str, owner: dict[str, tuple[str, str]], here: str) -> str:
    """Make a page's links right for where it now lives. ``root`` leads to the site root."""

    def page_href(slug: str) -> str:
        return f"{root}docs/{slug + '/' if slug else ''}"

    def anchor(match: re.Match) -> str:
        name = match.group(1)
        if name not in owner:
            return match.group(0)
        slug, _ = owner[name]
        if name in PRIMARY:
            return f'href="{page_href(slug)}"'  # a menu entry: always the page itself
        return f'href="{"" if slug == here else page_href(slug)}#{name}"'

    text = re.sub(r'href="#([\w-]+)"', anchor, text)
    text = text.replace('href="./#', f'href="{root}docs/#')
    text = text.replace('href="./network.html"', f'href="{root}"')
    text = text.replace('href="./"', f'href="{root}docs/"')
    text = text.replace("./assets/", f"{root}assets/")
    text = text.replace("./map/", f"{root}map/")
    # A nav entry that stood in for a section on the long page has no section to watch.
    return re.sub(r'\s+data-spy="[^"]*"', "", text)


def mark_current(header: str, href: str) -> str:
    """Highlight this page's own nav entry; the scroll highlighter has nothing to watch."""
    header = header.replace(' aria-current="page"', "")
    return header.replace(f'<a href="{href}"', f'<a href="{href}" aria-current="page"', 1)


def retitle(head: str, title: str) -> str:
    full = "APPFL Bio Suite" if title == "APPFL Bio Suite" else f"{title} — APPFL Bio Suite"
    return re.sub(r"<title>.*?</title>", f"<title>{html.escape(full)}</title>", head, count=1)


def assemble(head: str, header: str, body: str, footer: str) -> str:
    return (
        f"{head}<body>\n\n{header}\n\n<main>\n\n{body}\n\n</main>\n\n{footer}\n\n</body>\n</html>\n"
    )


def build(site: Path, out: Path) -> list[Path]:
    written = []
    head, header, footer, sections = split((site / "index.html").read_text(encoding="utf-8"))
    # A section added to website/index.html must be given a page here; dropping it from
    # the live site without a word is the failure this refuses.
    placed = {key for _, keys in PAGES.values() for key in keys}
    unplaced = sorted(set(sections) - placed)
    if unplaced:
        raise SystemExit(
            f"website/index.html has sections no docs page holds: {', '.join(unplaced)}. "
            "Add each to PAGES in scripts/build_docs.py (its own page, or About's list)."
        )
    owner = anchor_pages(sections)

    for slug, (title, keys) in PAGES.items():
        root = "../" if not slug else "../../"
        body = "\n\n".join(sections[key] for key in keys if key in sections)
        page = assemble(retitle(head, title), header, body, footer)
        page = rewrite(page, root, owner, slug)
        here = f"{root}docs/{slug + '/' if slug else ''}"
        page = mark_current(page, here)
        target = out / "docs" / slug / "index.html"
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_text(page, encoding="utf-8")
        written.append(target)

    # The map's own guide, with the same header as the rest of the documentation.
    network = (site / "network.html").read_text(encoding="utf-8")
    n_head = network[: network.index("<body>")]
    main = re.search(r"<main>(.*?)</main>", network, re.S).group(1)
    root = "../../"
    page = assemble(n_head, header, main.strip(), footer)
    page = rewrite(page, root, owner, NETWORK_PAGE)
    page = mark_current(page, root)
    target = out / "docs" / NETWORK_PAGE / "index.html"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(page, encoding="utf-8")
    written.append(target)

    written.append(_root_map(out, owner))
    written.append(_redirect(out / "network.html", "./"))
    return written


def _root_map(out: Path, owner: dict[str, tuple[str, str]]) -> Path:
    """The site's front door: the map viewer, reading its data from map/.

    The viewer takes its data from ?metadata_url=; with none, it falls back to the
    data-bio-default-metadata attribute set here, so the address stays bare. Old links
    to a section of the long page (#about, #quickstart, ...) are sent to its new page.
    """
    viewer = (out / "map" / "map.html").read_text(encoding="utf-8")
    viewer, n = re.subn(
        r"<html\b", '<html data-bio-default-metadata="map/network.map.json"', viewer, count=1
    )
    if n != 1:
        raise SystemExit("map/map.html has no <html> tag")
    old = {name: ("docs/" + (slug + "/" if slug else "") + ("" if name in PRIMARY else f"#{name}"))
           for name, (slug, key) in owner.items()}  # fmt: skip
    redirect = (
        "<script>(function () { var to = " + json.dumps(old) + "[location.hash.slice(1)];"
        " if (to) location.replace(to); })();</script>"
    )
    viewer = viewer.replace("<head>", "<head>\n" + redirect, 1)
    target = out / "index.html"
    target.write_text(viewer, encoding="utf-8")
    return target


def _redirect(path: Path, to: str) -> Path:
    path.write_text(
        '<!doctype html>\n<meta charset="utf-8">\n<title>APPFL Bio Suite</title>\n'
        f'<meta http-equiv="refresh" content="0; url={to}">\n'
        f'<link rel="canonical" href="{to}">\n'
        f'<p><a href="{to}">The network map has moved.</a></p>\n',
        encoding="utf-8",
    )
    return path


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--site", type=Path, required=True, help="website/ source directory")
    parser.add_argument("--out", type=Path, required=True, help="the assembled site")
    args = parser.parse_args()
    for path in build(args.site, args.out):
        print(path.relative_to(args.out))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
