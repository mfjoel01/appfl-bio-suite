/* Experiment results library: curated figures, sortable tables and sandboxed reports. */
(() => {
  'use strict';
  const labels = { 'fine-mapping': 'Fine-mapping', gwas: 'GWAS',
    'flamby-heart-disease': 'FLamby · Heart disease', caidf: 'CAIDF', cpg: 'CpG', tbd: 'Project TBD' };
  const title = id => labels[id] || id;
  const make = (tag, className, text) => {
    const element = document.createElement(tag);
    if (className) element.className = className;
    if (text !== undefined) element.textContent = text;
    return element;
  };
  const safeAsset = value => typeof value === 'string'
    && /^data:(image\/(png|jpeg|webp)|text\/(html|csv|tab-separated-values));base64,[a-zA-Z0-9+/=\s]+$/.test(value);
  class BioResults {
    constructor(sidebar, workspace) {
      this.sidebar = sidebar; this.workspace = workspace; this.groups = [];
      this.experiments = []; this.selected = ''; this.kind = 'all';
      const intro = make('div', 'bio-panel-intro');
      intro.append(make('h1', '', 'Experiment results'), make('p', '', 'Figures, tables and reports in one place.'));
      const label = make('label', 'bio-results-label', 'EXPERIMENT');
      this.selector = make('select', 'bio-result-select'); this.selector.id = 'bio-result-experiment';
      label.htmlFor = this.selector.id;
      this.selector.onchange = () => { this.selected = this.selector.value; this.kind = 'all'; this.render(); };
      this.navigation = make('nav', 'bio-results-nav'); this.navigation.setAttribute('aria-label', 'Result types');
      this.note = make('p', 'bio-results-note');
      sidebar.append(intro, label, this.selector, this.navigation, this.note);
    }
    setData(groups, experiments) {
      this.groups = Array.isArray(groups) ? groups : [];
      this.experiments = [...new Set([...experiments, ...this.groups.map(g => g.experiment)])];
      if (!this.experiments.includes(this.selected))
        this.selected = this.groups[0]?.experiment || this.experiments[0] || '';
      this.selector.replaceChildren();
      for (const id of this.experiments) {
        const option = make('option', '', title(id)); option.value = id;
        this.selector.append(option);
      }
      this.selector.value = this.selected;
      this.selector.disabled = !this.experiments.length;
      this.render();
    }
    open(experiment) {
      if (experiment && this.experiments.includes(experiment)) this.selected = experiment;
      this.selector.value = this.selected; this.render();
    }
    render() {
      const groups = this.groups.filter(g => g.experiment === this.selected);
      const all = groups.flatMap(g => g.artifacts || []);
      this.navigation.replaceChildren();
      for (const [kind, name] of [['all','All results'],['image','Figures'],['table','Tables'],['report','Reports']]) {
        const count = all.filter(a => kind === 'all' || a.kind === kind).length;
        const button = make('button', '', name);
        button.append(make('span', '', count)); button.setAttribute('aria-pressed', String(kind === this.kind));
        button.onclick = () => { this.kind = kind; this.render(); };
        this.navigation.append(button);
      }
      this.note.textContent = 'Only results selected for this viewer appear here. Site selection on the map does not alter these experiment-level results.';
      const header = make('div', 'bio-results-header');
      header.append(make('div', 'bio-eyebrow', 'EXPERIMENT RESULTS'),
        make('h1', '', this.selected ? title(this.selected) : 'Results library'),
        make('p', '', all.length ? `${all.length} published artifacts · Select a figure to expand it or download a table.`
          : 'No results have been added for this experiment yet.'));
      this.workspace.replaceChildren(header);
      if (!all.length) {
        const empty = make('div', 'bio-results-empty');
        empty.append(make('span', '', '▥'), make('h2', '', 'Results will appear here'),
          make('p', '', 'This experiment’s figures, tables and reports will be available when they are added to the viewer.'));
        this.workspace.append(empty); return;
      }
      let shown = 0;
      for (const group of groups) {
        const artifacts = (group.artifacts || []).filter(a => this.kind === 'all' || a.kind === this.kind);
        if (!artifacts.length) continue;
        shown += artifacts.length;
        const section = make('section', 'bio-result-group');
        section.append(make('h2', '', group.title), make('p', 'bio-group-description', group.description || ''));
        const grid = make('div', 'bio-result-grid');
        for (const artifact of artifacts) grid.append(this.artifact(artifact));
        section.append(grid); this.workspace.append(section);
      }
      if (!shown) this.workspace.append(make('p', 'bio-results-empty', 'No artifacts of this type have been added.'));
    }
    artifact(artifact) {
      const card = make('article', `bio-result-card bio-result-${artifact.kind}`);
      const header = make('div', 'bio-artifact-header');
      header.append(make('h3', '', artifact.title));
      if (safeAsset(artifact.download)) {
        const download = make('a', 'bio-download', 'Download ↓'); download.href = artifact.download;
        download.download = artifact.filename || 'result'; header.append(download);
      }
      card.append(header);
      if (artifact.description) card.append(make('p', 'bio-artifact-description', artifact.description));
      if (artifact.kind === 'image' && safeAsset(artifact.download)) {
        const button = make('button', 'bio-figure-button');
        button.setAttribute('aria-label', `Expand ${artifact.title}`);
        const image = make('img', 'bio-result-image'); image.src = artifact.download;
        image.alt = artifact.title; image.loading = 'lazy'; button.append(image);
        button.onclick = () => {
          const dialog = make('dialog', 'bio-figure-dialog');
          const close = make('button', 'bio-dialog-close', 'Close ×');
          close.onclick = () => dialog.close();
          const expanded = make('img', ''); expanded.src = artifact.download; expanded.alt = artifact.title;
          dialog.append(close, make('h2', '', artifact.title), expanded);
          dialog.onclose = () => dialog.remove(); dialog.onclick = e => { if (e.target === dialog) dialog.close(); };
          document.body.append(dialog); dialog.showModal();
        };
        card.append(button);
      } else if (artifact.kind === 'table') this.table(card, artifact);
      else if (artifact.kind === 'report' && safeAsset(artifact.download)) {
        const open = make('button', 'bio-open-report', 'Open report');
        open.onclick = () => {
          const frame = make('iframe', 'bio-report-frame'); frame.title = artifact.title;
          frame.setAttribute('sandbox', ''); frame.referrerPolicy = 'no-referrer'; frame.src = artifact.download;
          card.append(frame); open.remove();
        };
        card.append(open);
      }
      return card;
    }
    table(card, artifact) {
      const controls = make('div', 'bio-table-controls');
      const search = make('input', ''); search.type = 'search'; search.placeholder = 'Search preview rows…';
      search.setAttribute('aria-label', `Search ${artifact.title}`);
      const count = make('span', ''); count.setAttribute('role', 'status'); controls.append(search, count);
      const scroller = make('div', 'bio-table-scroll'); scroller.tabIndex = 0;
      scroller.setAttribute('role', 'region'); scroller.setAttribute('aria-label', `${artifact.title} table`);
      const table = make('table', ''); table.append(make('caption', 'bio-sr-only', artifact.title));
      const thead = make('thead', ''); const row = make('tr', ''); const tbody = make('tbody', '');
      let sort = null, direction = 1;
      const columns = Array.isArray(artifact.columns) ? artifact.columns : [];
      const rows = Array.isArray(artifact.rows) ? artifact.rows : [];
      const render = () => {
        const query = search.value.toLocaleLowerCase();
        const filtered = rows.filter(r => r.some(c => String(c).toLocaleLowerCase().includes(query)));
        if (sort !== null) filtered.sort((a, b) => {
          const left = String(a[sort]), right = String(b[sort]);
          const numeric = left.trim() && right.trim() && Number.isFinite(Number(left)) && Number.isFinite(Number(right));
          return direction * (numeric ? Number(left) - Number(right) : left.localeCompare(right));
        });
        tbody.replaceChildren();
        for (const cells of filtered) {
          const tr = make('tr', ''); for (const cell of cells) tr.append(make('td', '', cell)); tbody.append(tr);
        }
        count.textContent = `${filtered.length} shown · ${artifact.total_rows ?? rows.length} total rows`;
        [...row.children].forEach((th, index) => th.setAttribute('aria-sort', sort === index
          ? direction === 1 ? 'ascending' : 'descending' : 'none'));
      };
      columns.forEach((column, index) => {
        const th = make('th', ''); th.scope = 'col'; const button = make('button', '', column + ' ↕');
        button.onclick = () => { direction = sort === index ? -direction : 1; sort = index; render(); };
        th.append(button); row.append(th);
      });
      thead.append(row); table.append(thead, tbody); scroller.append(table); card.append(controls, scroller);
      if (artifact.total_rows > rows.length) card.append(make('p', 'bio-table-note',
        `Preview limited to the first ${rows.length} rows. Download includes all ${artifact.total_rows} rows; search and sorting apply to the preview.`));
      search.oninput = render; render();
    }
  }
  window.BioResults = BioResults;
})();
