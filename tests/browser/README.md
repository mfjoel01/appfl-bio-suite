# Watch viewer browser regression

Use Node.js and Playwright with Chromium installed in your development environment.
These are test tools, not dependencies of the exported viewer.

From the repository root, export any example map and serve it:

```bash
appfl-bio-suite watch export --federation federation.yaml.example --out local/watch-test
python -m http.server 7081 --bind 127.0.0.1 --directory local/watch-test
```

In a second terminal:

```bash
WATCH_VIEWER_URL=http://127.0.0.1:7081/map.html node tests/browser/watch_viewer.cjs
```

The browser test supplies synthetic metadata and a controlled SSE stream, then uses the
actual exported viewer. It checks default flat view, experiment unions and empty
selection, matching map layers/counts, null and zero coordinates, globe interaction and
reduced motion, mobile layout, JSONL playback, late joining an active run, stream
isolation, and escaped metadata. The flat map still loads Leaflet and tiles from the
internet, just as the product does.

Optional environment variables:

- `PLAYWRIGHT_MODULE`: path to an existing Playwright or playwright-core installation.
- `CHROMIUM_EXECUTABLE`: path to an existing Chromium executable.
- `WATCH_SCREENSHOT_DIR`: directory for the mobile screenshot.

The Python tests in `tests/test_watch.py` independently cover export packaging, private
data boundaries, server lifecycle, and real HTTP run discovery after live publication.
