/* appfl-bio-suite project site: theme toggle and section highlighting.
   Both pages load this; each feature is skipped when its markup is absent. The theme
   is already resolved onto <html data-theme> by the inline script in each page's head,
   so nothing here runs before first paint. */
(() => {
  'use strict';

  const root = document.documentElement;
  const STORE_KEY = 'appfl-bio-suite:theme';

  // --- theme ---------------------------------------------------------------

  const toggle = document.getElementById('theme-toggle');
  const frame = document.getElementById('map-frame');

  /* The embedded viewer is same-origin and its globe already watches <html data-theme>
     for the sake of its own theme button, so following it is one attribute write -- no
     reload, and the globe redraws itself. Wrapped anyway: a viewer that failed to load
     should not take the toggle down with it. */
  function applyToFrame(theme) {
    try {
      const doc = frame && frame.contentDocument;
      if (doc && doc.documentElement) doc.documentElement.dataset.theme = theme;
    } catch (error) {
      /* Cross-origin or not yet loaded. The frame's own src carries the theme. */
    }
  }

  function setTheme(theme, remember) {
    root.dataset.theme = theme;
    if (toggle) {
      const next = theme === 'light' ? 'dark' : 'light';
      toggle.setAttribute('aria-label', `Switch to ${next} theme`);
      toggle.title = `Switch to ${next} theme`;
      toggle.setAttribute('aria-pressed', String(theme === 'light'));
    }
    applyToFrame(theme);
    if (remember) {
      try { localStorage.setItem(STORE_KEY, theme); } catch (error) { /* private mode */ }
    }
  }

  setTheme(root.dataset.theme === 'light' ? 'light' : 'dark', false);

  if (toggle) {
    toggle.addEventListener('click', () => {
      setTheme(root.dataset.theme === 'light' ? 'dark' : 'light', true);
    });
  }

  // A frame that finishes loading after a toggle, or is restored from cache, still lands
  // on the page's theme rather than the one baked into its src.
  if (frame) frame.addEventListener('load', () => applyToFrame(root.dataset.theme));

  // Follow the OS only while the visitor has never chosen for themselves.
  const scheme = window.matchMedia('(prefers-color-scheme: light)');
  const listen = scheme.addEventListener ? scheme.addEventListener.bind(scheme, 'change')
    : scheme.addListener.bind(scheme);
  listen(event => {
    let stored = null;
    try { stored = localStorage.getItem(STORE_KEY); } catch (error) { /* private mode */ }
    if (stored !== 'light' && stored !== 'dark') setTheme(event.matches ? 'light' : 'dark', false);
  });

  // --- which section am I reading ------------------------------------------

  const header = document.querySelector('.site-header');
  /* A nav entry claims a section either by its own #hash or, where the entry points at
     another page, by data-spy naming the section that stands in for it here -- which is
     how "Network map" lights up over the block that links to it. Sorted by document
     position rather than by nav order, so reordering the nav cannot silently break the
     sequence the scan below depends on. */
  const spied = [...document.querySelectorAll('.site-nav a[href^="#"], .site-nav a[data-spy]')]
    .map(link => ({ link, section: document.getElementById(link.dataset.spy || link.hash.slice(1)) }))
    .filter(pair => pair.section)
    .sort((a, b) =>
      a.section.compareDocumentPosition(b.section) & Node.DOCUMENT_POSITION_FOLLOWING ? -1 : 1);

  if (spied.length) {
    /* Nearest section whose top has passed under the header, rather than whichever
       intersects: sections here are taller than the viewport and several are visible at
       once, so "the one I am in" is the last one I scrolled past. */
    function spy() {
      const cutoff = (header ? header.offsetHeight : 0) + 28;
      let current = null;
      for (const pair of spied) {
        if (pair.section.getBoundingClientRect().top <= cutoff) current = pair.link;
      }
      // The last section is usually too short to reach the cutoff before the page runs
      // out of scroll, so the bottom of the document counts as being in it. Requires an
      // actual scroll: on a document shorter than the viewport every position is "the
      // bottom", and the last entry would light up while the reader is still at the top.
      const bottom = window.scrollY > 0
        && window.innerHeight + window.scrollY >= document.body.scrollHeight - 2;
      if (bottom) current = spied[spied.length - 1].link;

      for (const pair of spied) {
        if (pair.link === current) pair.link.setAttribute('aria-current', 'location');
        else if (pair.link.getAttribute('aria-current') === 'location') {
          pair.link.removeAttribute('aria-current');
        }
      }
    }

    let queued = false;
    function schedule() {
      if (queued) return;
      queued = true;
      requestAnimationFrame(() => { queued = false; spy(); });
    }

    addEventListener('scroll', schedule, { passive: true });
    addEventListener('resize', schedule);
    spy();
  }
})();
