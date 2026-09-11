/* Suite interface for the pinned HiveWatch viewer. No federation data is sent elsewhere. */
(() => {
  'use strict';
  const $ = id => document.getElementById(id);
  const esc = escapeHtml;
  const ui = { selected: null, experiments: [], view: 'flat', tab: 'experiments',
    server: null, globe: null, runId: null, loading: false, loadToken: 0, focus: null,
    liveRuns: new Set(), queued: [], error: '', resultGroups: [], resultsRevision: 0, resultsSignature: '',
    source: !!(METADATA_URL || EVENTS_URL) };
  const labels = { 'fine-mapping': 'Fine-mapping', gwas: 'GWAS',
    'flamby-heart-disease': 'FLamby · Heart disease', caidf: 'CAIDF', cpg: 'CpG', tbd: 'Project TBD' };
  const label = name => labels[name] || name;
  const number = value => Number(value || 0).toLocaleString();
  const coordinate = v => v == null || (typeof v === 'string' && v.trim() === '') || typeof v === 'boolean' ? null
    : Number.isFinite(Number(v)) ? Number(v) : null;
  const placed = c => coordinate(c.lat) !== null && coordinate(c.lng) !== null
    && Math.abs(Number(c.lat)) <= 90 && Math.abs(Number(c.lng)) <= 180;
  function memberships(c) {
    const value = c.experiments ?? c.experiment ?? c.extra?.experiments ?? c.extra?.experiment ?? '';
    return [...new Set((Array.isArray(value) ? value : String(value).split(','))
      .map(x => String(x).trim()).filter(Boolean))];
  }
  const matches = c => ui.selected === null || memberships(c).some(e => ui.selected.has(e));
  const visible = () => Object.entries(state.clients).filter(([, c]) => matches(c));
  const color = c => ({ '0': '#8d9eaa', '1': '#669fe0', '2': '#d89e42', '4': '#15b9a6', X: '#db717e' }
    [String(c.partnership_stage || '').slice(0, 1)] || { idle: '#f5bd66', dropped: '#f08085', failed: '#f08085' }[c.status] || '#50d9c7');
  const statusLabel = c => c.partnership_stage || (ui.runId === 'network' ? `Member · ${c.status || 'active'}` : c.status || 'active');

  document.documentElement.dataset.theme = 'dark';
  document.body.classList.add('bio-viewer');
  document.title = 'Federation network · Hive Watch';
  const hiveLogo = document.querySelector('.logo-img');
  const branding = JSON.parse($('bio-branding-data')?.textContent || '{}');
  document.querySelector('.logo').innerHTML = `<img class="bio-suite-logo" alt="APPFL Bio Suite network logo">
    <div><div class="bio-brand">APPFL <span>Bio Suite</span></div>
    <div class="bio-brand-note">FEDERATION NETWORK</div></div>`;
  document.querySelector('.bio-suite-logo').src = branding.suite_logo || '';
  if (hiveLogo) {
    hiveLogo.classList.add('bio-hive-logo'); hiveLogo.alt = 'HiveWatch';
    const credit = document.createElement('div'); credit.className = 'bio-hive-credit';
    credit.append(hiveLogo); document.querySelector('.header-right').prepend(credit);
  }
  for (const id of ['hdr-round', 'hdr-acc', 'hdr-loss']) $(id).closest('.stat').hidden = true;
  $('hdr-clients').nextElementSibling.textContent = 'Sites shown';
  document.querySelector('.header-stats').insertAdjacentHTML('beforeend', `
    <div class="stat"><div class="stat-val" id="bio-experiment-count">0</div><div class="stat-label">Experiments</div></div>
    <div class="stat"><div class="stat-val" id="bio-country-count">0</div><div class="stat-label">Countries</div></div>`);
  $('theme-btn').textContent = '◐';
  $('theme-btn').setAttribute('aria-label', 'Toggle color theme');
  $('theme-btn').title = 'Toggle color theme';
  $('play-btn').setAttribute('aria-label', 'Play or pause run');
  $('speed-select').setAttribute('aria-label', 'Playback speed');

  const sidebar = document.querySelector('.sidebar');
  const runs = document.querySelector('.runs-panel');
  const roundSummary = sidebar.querySelector('.sidebar-section');
  const siteHeading = sidebar.querySelectorAll('.sidebar-section')[1];
  const logPanel = document.querySelector('.log-panel');
  const clientList = $('client-list');
  siteHeading.remove();
  sidebar.insertAdjacentHTML('afterbegin', `<div class="bio-tabs" role="tablist" aria-label="Explore federation">
    <button id="bio-experiments-tab" role="tab" aria-selected="true" aria-controls="bio-experiments-panel">Experiments</button>
    <button id="bio-results-tab" role="tab" aria-selected="false" aria-controls="bio-results-panel" tabindex="-1">Results</button>
    <button id="bio-runs-tab" role="tab" aria-selected="false" aria-controls="bio-runs-panel" tabindex="-1">Runs <span>↗</span></button></div>
    <section id="bio-experiments-panel" role="tabpanel" aria-labelledby="bio-experiments-tab">
      <div class="bio-panel-intro"><h1>Explore the network</h1><p>See where each experiment happens.</p></div>
      <fieldset id="bio-filters"><legend class="bio-sr-only">Filter sites by experiment</legend></fieldset>
      <div class="bio-sites-heading"><span>PARTICIPATING SITES</span><span id="bio-site-count">0</span></div>
      <p id="bio-filter-note" role="status"></p>
    </section>
    <section id="bio-results-panel" role="tabpanel" aria-labelledby="bio-results-tab" hidden></section>
    <section id="bio-runs-panel" role="tabpanel" aria-labelledby="bio-runs-tab" hidden>
      <div class="bio-panel-intro"><h1>Runs & playback</h1><p>Inspect live activity or replay a saved run.</p></div>
      <button id="bio-network-return" class="bio-text-button">← Return to federation network</button>
    </section>`);
  $('bio-experiments-panel').append(clientList);
  $('bio-runs-panel').append(runs, roundSummary, logPanel);

  const area = document.querySelector('.map-area');
  const resultsWorkspace = document.createElement('div'); resultsWorkspace.id = 'bio-results-workspace';
  resultsWorkspace.hidden = true; resultsWorkspace.setAttribute('aria-label', 'Experiment results');
  area.append(resultsWorkspace);
  const resultsView = new BioResults($('bio-results-panel'), resultsWorkspace);
  area.insertAdjacentHTML('afterbegin', `<div class="bio-map-toolbar">
    <div><div class="bio-eyebrow">CONNECTED SCIENCE</div><h2 id="bio-map-title">Federation footprint</h2></div>
    <div class="bio-view-switch" role="group" aria-label="Map view">
      <button id="bio-flat" aria-pressed="true">▱ <span>Flat map</span></button>
      <button id="bio-globe" aria-pressed="false">◎ <span>Globe</span></button>
    </div></div>
    <div id="bio-globe-container" hidden></div>
    <div id="bio-map-message" role="status" hidden></div>
    <div id="bio-site-detail" hidden></div>
    <div class="bio-map-footer"><div class="bio-legend"><span><i class="bio-dot"></i> Site</span>
      <span><i class="bio-dot bio-hub"></i> Coordinator</span></div>
      <span id="bio-geography-note">Institution locations · See site details for sources</span></div>
    <div id="bio-globe-controls" hidden><button id="bio-spin" aria-pressed="true">Pause rotation</button>
      <button id="bio-zoom-in" aria-label="Zoom globe in">+</button>
      <button id="bio-zoom-out" aria-label="Zoom globe out">−</button>
      <button id="bio-reset" aria-label="Reset globe view">↺</button>
      <span>Drag to rotate · Scroll to zoom</span></div>`);
  $('debug').hidden = true;

  function setTab(tab) {
    ui.tab = tab;
    document.body.classList.toggle('bio-show-runs', tab === 'runs');
    document.body.classList.toggle('bio-show-results', tab === 'results');
    resultsWorkspace.hidden = tab !== 'results';
    ui.globe?.setVisible(ui.view === 'globe' && tab !== 'results');
    if (tab === 'results') resultsView.open(ui.selected?.size === 1 ? [...ui.selected][0] : null);
    for (const name of ['experiments', 'results', 'runs']) {
      const active = tab === name;
      $(`bio-${name}-tab`).setAttribute('aria-selected', String(active));
      $(`bio-${name}-tab`).tabIndex = active ? 0 : -1;
      $(`bio-${name}-panel`).hidden = !active;
    }
    document.querySelector('.playback-bar').hidden = tab !== 'runs';
    $('debug').hidden = tab !== 'runs' || ui.source;
    syncLayers();
    requestAnimationFrame(() => map.invalidateSize());
  }
  const tabs = ['experiments', 'results', 'runs'];
  for (const tab of tabs) {
    $(`bio-${tab}-tab`).onclick = () => setTab(tab);
    $(`bio-${tab}-tab`).onkeydown = event => {
      if (!['ArrowLeft', 'ArrowRight', 'Home', 'End'].includes(event.key)) return;
      event.preventDefault();
      const next = event.key === 'Home' ? 'experiments' : event.key === 'End' ? 'runs'
        : tabs[(tabs.indexOf(tab) + (event.key === 'ArrowRight' ? 1 : 2)) % tabs.length];
      setTab(next); $(`bio-${next}-tab`).focus();
    };
  }
  setTab('experiments');
  initDraggableLogPanel();

  function discoverExperiments() {
    const clients = [...Object.values(state.clients), ...pb.rounds.flatMap(r => r.clients || [])];
    ui.experiments = [...new Set(clients.flatMap(memberships))].sort((a, b) => {
      const order = ['fine-mapping', 'gwas', 'flamby-heart-disease'];
      return (order.indexOf(a) < 0 ? 99 : order.indexOf(a))
        - (order.indexOf(b) < 0 ? 99 : order.indexOf(b)) || a.localeCompare(b);
    });
    const signature = JSON.stringify(ui.experiments);
    if ($('bio-filters').dataset.signature !== signature) {
      $('bio-filters').dataset.signature = signature;
      $('bio-filters').innerHTML = '<legend class="bio-sr-only">Filter sites by experiment</legend>';
      for (const name of [null, ...ui.experiments]) {
        const row = document.createElement('label'); row.className = 'bio-filter';
        const input = document.createElement('input'); input.type = 'checkbox';
        input.dataset.experiment = name || '*';
        const text = document.createElement('span'); text.textContent = name === null ? 'All experiments' : label(name);
        const count = document.createElement('span'); count.className = 'bio-filter-count';
        row.append(input, text, count); $('bio-filters').append(row);
        input.onchange = () => {
          if (name === null) ui.selected = input.checked ? null : new Set();
          else {
            if (ui.selected === null) ui.selected = new Set();
            if (input.checked) ui.selected.add(name); else ui.selected.delete(name);
          }
          closeDetail(); renderSidebar();
        };
      }
    }
    for (const input of $('bio-filters').querySelectorAll('input')) {
      const name = input.dataset.experiment;
      input.checked = name === '*' ? ui.selected === null : !!ui.selected?.has(name);
      input.closest('label').classList.toggle('selected', input.checked);
      input.closest('label').querySelector('.bio-filter-count').textContent = Object.values(state.clients)
        .filter(c => name === '*' || memberships(c).includes(name)).length;
    }
  }

  function cardHtml(id, c) {
    return `<button class="bio-site-card${ui.focus === id ? ' selected' : ''}" data-site="${esc(id)}">
      <span class="bio-card-top"><span class="bio-site-name">${esc(c.institution || id)}</span>
        <span class="bio-dot" style="background:${color(c)}" title="${esc(statusLabel(c))}"></span></span>
      <span class="bio-site-location">${esc([c.city, c.country].filter(Boolean).join(', ') || 'Location not declared')}</span>
      <span class="bio-tags">${memberships(c).map(e => `<span>${esc(label(e))}</span>`).join('')}</span>
      ${c.partnership_stage ? `<span class="bio-partner-stage" style="--stage-color:${color(c)}">${esc(c.partnership_stage)}</span>` : ''}
      <span class="bio-card-bottom"><span>${c.num_samples == null ? 'Samples not declared' : number(c.num_samples) + ' declared samples'}</span>
        <span>${placed(c) ? '↗' : 'Unplaced'}</span></span></button>`;
  }
  function detailsHtml(id, c) {
    const experimentRows = memberships(c).map(e => `<div class="bio-detail-row"><span>${esc(label(e))}</span>
      <strong>${esc(c[e] || 'Participating')}</strong></div>`).join('');
    return `<div class="popup-title">${esc(c.institution || id)}</div>
      <p class="bio-detail-location">${esc([c.city, c.country].filter(Boolean).join(', ') || 'Location not declared')}</p>
      ${experimentRows}<div class="bio-detail-row"><span>Samples across experiments</span><strong>${c.num_samples == null ? 'Not declared' : number(c.num_samples)}</strong></div>
      <div class="bio-detail-row"><span>${c.partnership_stage ? 'Project readiness' : 'Status'}</span><strong>${esc(statusLabel(c))}</strong></div>
      ${Array.isArray(c.contacts) ? c.contacts.map(contact => `<div class="bio-detail-row"><span>${esc(contact.name)}</span>
        <a href="mailto:${esc(encodeURI(contact.email))}">${esc(contact.email)}</a></div>`).join('') : ''}
      ${c.notes ? `<div class="bio-detail-row"><span>Notes</span><strong>${esc(c.notes)}</strong></div>` : ''}
      ${c.location_basis ? `<div class="bio-detail-row"><span>Map location</span><strong>${esc(c.location_basis)}</strong>
        ${/^https?:\/\//.test(c.location_source || '') ? `<a href="${esc(c.location_source)}" target="_blank" rel="noopener noreferrer">Location source ↗</a>` : ''}</div>` : ''}
      <details class="bio-metadata"><summary>Site details & metrics</summary><div class="client-metrics">${renderClientMetrics(c)}</div></details>`;
  }
  function closeDetail() { ui.focus = null; $('bio-site-detail').hidden = true; map.closePopup(); }
  focusClient = id => {
    const c = state.clients[id];
    if (!c || !matches(c)) return;
    ui.focus = id;
    if (placed(c)) {
      if (ui.view === 'globe') ui.globe.focus(Number(c.lat), Number(c.lng));
      else { map.flyTo([c.lat, c.lng], 5, { duration: 0.8 }); markers[id]?.openPopup(); }
    }
    if (ui.view === 'globe' || !placed(c)) {
      const panel = $('bio-site-detail'); panel.hidden = false;
      panel.innerHTML = '<button class="bio-detail-close" aria-label="Close site details">×</button>' + detailsHtml(id, c);
      panel.querySelector('button').onclick = () => { closeDetail(); renderSidebar(); };
    }
    renderSidebar();
  };
  clientList.onclick = event => {
    const card = event.target.closest('[data-site]'); if (card) focusClient(card.dataset.site);
  };

  // Keep complete HiveWatch state for playback; filtering changes only visible layers.
  function syncLayers() {
    for (const [id, c] of Object.entries(state.clients)) {
      const show = matches(c) && placed(c) && ui.view === 'flat' && ui.tab !== 'results';
      const toggle = (layer, enabled) => {
        if (!layer) return;
        if (enabled && !map.hasLayer(layer)) layer.addTo(map);
        else if (!enabled && map.hasLayer(layer)) layer.remove();
      };
      toggle(markers[id], show);
      toggle(lines[id], show && !!ui.server);
      toggle(packets[id]?.uplink, show && !!ui.server && ui.runId !== 'network');
      toggle(packets[id]?.downlink, show && !!ui.server && ui.runId !== 'network');
    }
    if (ui.view === 'globe' || ui.runId === 'network' || ui.tab === 'results') stopNetworkAnimation();
    if (ui.globe) ui.globe.setData(visible().map(([id, c]) => ({ ...c, client_id: id })), ui.server);
  }

  renderSidebar = () => {
    discoverExperiments();
    const entries = visible();
    const total = Object.keys(state.clients).length;
    $('hdr-clients').textContent = entries.length;
    $('bio-experiment-count').textContent = ui.experiments.length;
    $('bio-country-count').textContent = new Set(entries.map(([, c]) => c.country).filter(Boolean)).size;
    $('bio-site-count').textContent = `${entries.length} / ${total}`;
    $('sb-round').textContent = state.round || '—';
    $('sb-participants').textContent = total;
    $('sb-acc').textContent = fmtAcc(state.globalAcc);
    $('sb-loss').textContent = fmtLoss(state.globalLoss);
    const unplaced = entries.filter(([, c]) => !placed(c)).length;
    $('bio-filter-note').textContent = ui.selected === null ? 'All sites and partner institutions'
      : ui.selected.size ? 'Sites in any selected experiment' : 'Select an experiment to show its sites.';
    if (unplaced) $('bio-filter-note').textContent += ` · ${unplaced} unplaced`;
    const focused = document.activeElement?.dataset.site;
    clientList.innerHTML = entries.map(([id, c]) => cardHtml(id, c)).join('')
      || `<div class="bio-empty"><span>◎</span><h3>${ui.loading ? 'Loading the network…' : 'No sites to show'}</h3>
        <p>${esc(ui.error || (total ? 'Choose All experiments or another experiment.' : 'Build a network map or select a run.'))}</p></div>`;
    if (focused) [...clientList.querySelectorAll('[data-site]')].find(b => b.dataset.site === focused)?.focus();
    $('bio-map-message').hidden = !!entries.length || ui.loading;
    $('bio-map-message').textContent = ui.error || (total ? 'No sites match this selection' : 'No federation data loaded');
    syncLayers();
    const signature = `${ui.resultsRevision}|${ui.experiments.join(',')}`;
    if (ui.resultsSignature !== signature) {
      ui.resultsSignature = signature; resultsView.setData(ui.resultGroups, ui.experiments);
    }
  };

  // HiveWatch 0.2.1 coerces null to zero and guesses an undeclared coordinator.
  // Both conflict with this suite's declared-geography contract.
  toFiniteNumber = coordinate;
  inferServerFromClient = () => {};
  resolveServerLocation = async () => false;
  applyServerMetadata = server => {
    ui.server = server && placed(server) ? { ...server, lat: Number(server.lat), lng: Number(server.lng) } : null;
    if (serverMarker) { serverMarker.remove(); serverMarker = null; }
    serverResolvedFromMetadata = !!ui.server;
    if (ui.server) {
      SERVER_LL = [ui.server.lat, ui.server.lng];
      serverMarker = L.marker(SERVER_LL, { icon: L.divIcon({ className: 'bio-leaflet-hub',
        html: '<span></span>', iconSize: [20, 20], iconAnchor: [-6, 26] }), zIndexOffset: 1000 })
        .bindPopup(`<div class="popup-title">Coordinator</div><p>${esc(server.org || 'Coordinator')}</p>
        <p>${esc(server.city || '')}</p>`).addTo(map);
      for (const [id, line] of Object.entries(lines)) {
        const c = state.clients[id]; if (c && placed(c)) line.setLatLngs([[c.lat, c.lng], SERVER_LL]);
      }
    }
    syncLayers();
    return !!ui.server;
  };
  for (const dictionary of [state.clients, markers, lines, packets]) Object.setPrototypeOf(dictionary, null);
  addLog = (message, type = 'update') => {
    const entry = document.createElement('div'); entry.className = `log-entry ${type}`;
    const time = document.createElement('span'); time.className = 'ts';
    time.textContent = new Date().toLocaleTimeString('en', { hour12: false });
    const text = document.createElement('span'); text.className = 'msg'; text.textContent = message;
    entry.append(time, text); $('log-entries').prepend(entry);
    if ($('log-entries').children.length > 50) $('log-entries').lastElementChild.remove();
    $('log-count').textContent = `${++state.logCount} events`;
  };
  buildPacketTooltip = (id, direction) => {
    const c = state.clients[id] || {};
    const rows = { Round: state.round || '—', 'Grad norm': fmtNumber(c.gradient_norm, 2),
      'Grad magnitude': fmtNumber(c.gradient_magnitude, 4), Accuracy: fmtAcc(c.local_accuracy),
      Loss: fmtLoss(c.local_loss), Samples: c.num_samples ?? '—' };
    return `<div class="packet-meta"><div class="packet-title">${esc(id)}</div>
      <div class="packet-subtitle">${direction === 'uplink' ? 'Client to Server' : 'Server to Client'}</div>
      ${Object.entries(rows).map(([key, value]) => `<div class="packet-row"><span class="packet-key">${esc(key)}</span>
        <span class="packet-value">${esc(value)}</span></div>`).join('')}</div>`;
  };
  upsertClient = (id, data) => {
    const c = state.clients[id] = { ...state.clients[id], ...data };
    c.lat = coordinate(c.lat); c.lng = coordinate(c.lng);
    if (!placed(c)) return;
    const col = color(c);
    const ll = [c.lat, c.lng];
    const icon = L.divIcon({ className: 'bio-leaflet-marker',
      html: `<span style="--site-color:${col}"></span>`, iconSize: [20, 20], iconAnchor: [10, 10] });
    if (markers[id]) markers[id].setLatLng(ll).setIcon(icon).setPopupContent(detailsHtml(id, c));
    else markers[id] = L.marker(ll, { icon }).bindPopup(detailsHtml(id, c)).addTo(map);
    const style = { color: col, weight: 1.4, opacity: 0.5, dashArray: '4 8' };
    if (lines[id]) lines[id].setLatLngs([ll, SERVER_LL]).setStyle(style);
    else lines[id] = L.polyline([ll, SERVER_LL], style).addTo(map);
    ensureClientPackets(id, getPacketClientLL(id, c), col);
  };
  // JSONL fallback includes canonical round_end clients, even with no client_update events.
  buildRounds = events => {
    const rounds = new Map();
    for (const event of events) {
      if (event.round == null) continue;
      if (!rounds.has(event.round)) rounds.set(event.round, { round: event.round, clientsById: new Map() });
      const round = rounds.get(event.round);
      for (const c of event.clients || []) round.clientsById.set(c.client_id,
        { ...round.clientsById.get(c.client_id), ...Object.fromEntries(Object.entries(c).filter(([, v]) => v != null)) });
      if (event.event_type === 'comm_failure' && event.client_id) round.clientsById.set(event.client_id,
        { ...round.clientsById.get(event.client_id), client_id: event.client_id, status: 'failed' });
      if (event.event_type === 'round_end') {
        const rm = event.round_metrics || {};
        for (const [target, source] of Object.entries({ globalAcc: 'global_accuracy', globalLoss: 'global_loss',
          duration: 'round_duration_sec', divergence: 'gradient_divergence' })) {
          if (rm[source] != null) round[target] = rm[source];
        }
      }
    }
    return [...rounds.values()].sort((a, b) => a.round - b.round).map(({ clientsById, ...rest }) =>
      ({ ...rest, clients: [...clientsById.values()] }));
  };

  let spinning = !matchMedia('(prefers-reduced-motion: reduce)').matches;
  function spinChanged(value) {
    spinning = value; $('bio-spin').textContent = value ? 'Pause rotation' : 'Resume rotation';
    $('bio-spin').setAttribute('aria-pressed', String(value));
  }
  function setView(view) {
    closeDetail();
    if (view === 'globe' && !ui.globe) {
      try {
        ui.globe = new BioGlobe($('bio-globe-container'), { onSelect: focusClient, onSpinChange: spinChanged });
        ui.globe.setSpinning(spinning);
      } catch (error) {
        $('bio-map-message').textContent = 'Globe unavailable. You can continue using the flat map.';
        $('bio-map-message').hidden = false; return;
      }
    }
    ui.view = view;
    $('map').hidden = view !== 'flat';
    $('bio-globe-container').hidden = view !== 'globe';
    $('bio-globe-controls').hidden = view !== 'globe';
    $('bio-flat').setAttribute('aria-pressed', String(view === 'flat'));
    $('bio-globe').setAttribute('aria-pressed', String(view === 'globe'));
    area.classList.toggle('bio-globe-view', view === 'globe');
    ui.globe?.setVisible(view === 'globe' && ui.tab !== 'results');
    if (view === 'flat') requestAnimationFrame(() => map.invalidateSize());
    renderSidebar();
  }
  $('bio-flat').onclick = () => setView('flat');
  $('bio-globe').onclick = () => setView('globe');
  $('bio-spin').onclick = () => { spinChanged(!spinning); ui.globe?.setSpinning(spinning); };
  $('bio-zoom-in').onclick = () => ui.globe?.zoom(0.15);
  $('bio-zoom-out').onclick = () => ui.globe?.zoom(-0.15);
  $('bio-reset').onclick = () => ui.globe?.reset();
  spinChanged(spinning);
  document.addEventListener('keydown', event => { if (event.key === 'Escape') closeDetail(); });

  async function json(url) {
    const response = await fetch(url);
    if (!response.ok) throw new Error(`Could not load map data (HTTP ${response.status}).`);
    return response.json();
  }
  function fitNetwork() {
    const points = visible().map(([, c]) => c).filter(placed).map(c => [Number(c.lat), Number(c.lng)]);
    if (ui.server) points.push([ui.server.lat, ui.server.lng]);
    if (points.length) map.fitBounds(points, { paddingTopLeft: [55, 110], paddingBottomRight: [55, 75], maxZoom: 4 });
  }
  async function loadSource({ metadataUrl, eventsUrl, runId, live = false }) {
    const token = ++ui.loadToken;
    ui.loading = true; ui.error = ''; ui.runId = runId || 'external'; ui.queued = [];
    ui.selected = null; ui.resultGroups = []; ui.resultsRevision++;
    closeDetail(); stopPlayback(); clearMap(); applyServerMetadata(null);
    pb.rounds = []; pb.events = []; pb.runId = ui.runId; pb.index = 0;
    setMode('replay'); renderSidebar();
    try {
      let metadata = null, events = [];
      if (metadataUrl) {
        try { metadata = await json(metadataUrl); }
        catch (error) { if (!eventsUrl) throw error; }
      }
      if (!metadata?.rounds?.length && eventsUrl) events = await json(eventsUrl);
      if (token !== ui.loadToken) return;
      ui.resultGroups = metadata?.results || []; ui.resultsRevision++;
      live = !ui.source && ui.runId !== 'network' && !(metadata?.finished_at || events.some(e => e.event_type === 'finished'));
      if (live) ui.liveRuns.add(ui.runId); else ui.liveRuns.delete(ui.runId);
      pb.events = events;
      pb.rounds = metadata?.rounds?.length ? metadata.rounds : buildRounds(events);
      applyServerMetadata(metadata?.server || events.find(e => e.event_type === 'server_metadata')?.server);
      pb.rounds.forEach(applyRound);
      pb.index = pb.rounds.length;
      state.totalRounds = pb.rounds.length;
      $('play-btn').disabled = !pb.rounds.length || live;
      $('playback-mode').textContent = live ? 'live' : `replay · ${ui.runId}`;
      $('bio-map-title').textContent = ui.runId === 'network' ? 'Federation footprint' : `Run · ${ui.runId}`;
      if (live) setMode('live');
      if (ui.source) $('status-text').textContent = 'Network snapshot';
      buildRoundTicks(); updateProgress(); fitNetwork();
    } catch (error) {
      if (token !== ui.loadToken) return;
      ui.error = error.message || 'Unable to load the selected map.';
      $('status-text').textContent = 'Data unavailable';
      addLog(ui.error, 'drop');
    } finally {
      if (token === ui.loadToken) {
        ui.loading = false;
        const queued = ui.queued; ui.queued = [];
        for (const event of queued) receive(event);
        renderSidebar();
      }
    }
  }
  loadRunFromSource = loadSource;
  selectRun = async (runId, isLive) => {
    if (ui.source) return;
    const path = `${SERVER}/runs/${encodeURIComponent(runId)}`;
    await loadSource({ metadataUrl: `${path}/metadata`, eventsUrl: `${path}/events`, runId, live: isLive });
    if (ui.runId !== runId) return;
    for (const button of $('runs-list').querySelectorAll('[data-run]'))
      button.classList.toggle('active', button.dataset.run === runId);
  };
  loadRuns = async () => {
    if (ui.source) {
      $('runs-list').innerHTML = '<p class="bio-run-note">This export contains the network snapshot. Open the local viewer to browse runs and live activity.</p>';
      return [];
    }
    try {
      const runs = await json(`${SERVER}/runs`);
      $('runs-list').innerHTML = '';
      for (const run of runs) {
        const button = document.createElement('button');
        button.className = `run-item${run.run_id === ui.runId ? ' active' : ''}`;
        button.dataset.run = run.run_id;
        const live = ui.liveRuns.has(run.run_id);
        button.innerHTML = `<span class="run-id-text">${esc(run.run_id === 'network' ? 'Federation network' : run.run_id)}</span>
          <span class="run-meta">${esc(run.algorithm || '')}${live ? ' · live' : ''}</span>`;
        button.onclick = () => selectRun(run.run_id, ui.liveRuns.has(run.run_id));
        $('runs-list').append(button);
      }
      if (!runs.length) $('runs-list').innerHTML = '<p class="bio-run-note">No runs yet. Build a network map to get started.</p>';
      return runs;
    } catch (_) {
      $('runs-list').innerHTML = '<p class="bio-run-note">Unable to load runs. Use refresh to retry.</p>';
      return [];
    }
  };
  $('bio-network-return').onclick = () => { if (!ui.source) selectRun('network', false); setTab('experiments'); };

  function receive(msg) {
    if (msg.event_type === 'init') {
      if (msg.run_id && msg.run_id !== 'network' && msg.mode !== 'static') ui.liveRuns.add(msg.run_id);
      loadRuns();
      return;
    }
    if (msg.event_type === 'finished') { ui.liveRuns.delete(msg.run_id); loadRuns(); }
    if (!msg.run_id || msg.run_id !== ui.runId) return;
    if (ui.loading) { ui.queued.push(msg); return; }
    // A standing network can be rebuilt while other experiment runs keep streaming.
    if (ui.runId === 'network' && msg.event_type === 'round_end') {
      clearMap(); pb.events = []; pb.rounds = [];
    } else if (state.mode !== 'live' && ui.runId !== 'network') return;
    // Seed existing replay snapshots so newly streamed rounds retain their history.
    if (!pb.events.length && pb.rounds.length) {
      pb.events = pb.rounds.map(r => ({ event_type: 'round_end', run_id: ui.runId, round: r.round,
        clients: r.clients, round_metrics: { global_accuracy: r.globalAcc, global_loss: r.globalLoss,
          round_duration_sec: r.duration, gradient_divergence: r.divergence } }));
    }
    if (msg.event_type === 'comm_failure' && msg.client_id) {
      msg = { ...msg, event_type: 'client_update', clients: [{ ...state.clients[msg.client_id],
        client_id: msg.client_id, status: 'failed' }] };
    }
    handleLiveEvent(msg);
    pb.index = pb.rounds.length; updateProgress(); renderSidebar();
    if (msg.event_type === 'finished' && ui.runId !== 'network') {
      setMode('replay'); $('play-btn').disabled = !pb.rounds.length;
      $('playback-mode').textContent = `replay · ${ui.runId}`;
    }
  }
  connect = () => {
    const stream = new EventSource(SSE_URL);
    stream.onopen = () => {
      dbg('connected'); $('status-text').textContent = 'Connected';
      $('status-dot').classList.add('live'); loadRuns();
    };
    stream.onmessage = event => { try { receive(JSON.parse(event.data)); } catch (error) { console.error(error); } };
    stream.onerror = () => {
      dbg('disconnected'); $('status-text').textContent = 'Reconnecting…'; $('status-dot').classList.remove('live');
      // EventSource retries automatically without creating duplicate subscriptions.
    };
    window.addEventListener('pagehide', () => stream.close(), { once: true });
  };
  // Kept small and read-only for browser regression tests and downstream integrations.
  window.bioWatch = { memberships, placed, getState: () => ({ view: ui.view, tab: ui.tab,
    selected: ui.selected === null ? null : [...ui.selected], runId: ui.runId,
    visibleSites: visible().map(([id]) => id), experiments: [...ui.experiments] }) };
  renderSidebar();
  if (ui.source) {
    loadRuns();
    loadSource({ metadataUrl: METADATA_URL, eventsUrl: EVENTS_URL,
      runId: AUTO_RUN_ID || (METADATA_URL?.includes('network.map.json') ? 'network' : 'external') });
  } else {
    connect();
    loadRuns().then(runs => {
      const id = AUTO_RUN_ID || (runs.some(r => r.run_id === 'network') ? 'network' : runs[0]?.run_id);
      if (id) selectRun(id, ui.liveRuns.has(id));
    });
  }
})();
