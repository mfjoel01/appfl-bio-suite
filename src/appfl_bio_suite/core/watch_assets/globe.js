/* Interactive, offline orthographic globe for the Bio Suite watch viewer. */
(function () {
  'use strict';

  const TAU = Math.PI * 2;
  const COLORS = {
    idle: '#a1b7c5', registered: '#a1b7c5', training: '#56dfd3',
    active: '#56dfd3', connected: '#56dfd3', completed: '#76dca0',
    done: '#76dca0', error: '#ff8e91', failed: '#ff8e91', dropped: '#ff8e91', offline: '#8294a4',
  };
  const clamp = (value, low, high) => Math.max(low, Math.min(high, value));
  const coordinate = value => (typeof value === 'number' ||
    (typeof value === 'string' && value.trim() !== '')) && Number.isFinite(Number(value));
  function hasLocation(item) {
    return item && coordinate(item.lat) && coordinate(item.lng) &&
      Math.abs(Number(item.lat)) <= 90 && Math.abs(Number(item.lng)) <= 180;
  }

  class BioGlobe {
    constructor(container, options = {}) {
      if (!window.d3 || !window.BIO_LAND) throw new Error('Globe assets are unavailable.');
      this.container = container;
      this.options = options;
      this.clients = [];
      this.server = null;
      this.visible = false;
      this.rotation = [62, -28, 0];
      this.zoomFactor = 1;
      this.motionPreference = window.matchMedia('(prefers-reduced-motion: reduce)');
      this.theme = document.documentElement.dataset.theme || 'dark';
      this.spinning = !this.motionPreference.matches;
      this.frame = null;
      this.lastFrame = null;
      this.hitTargets = [];
      this.hovered = null;
      this.drag = null;
      this.width = 0;
      this.height = 0;
      this.canvas = document.createElement('canvas');
      this.canvas.className = 'bio-globe-canvas';
      this.canvas.tabIndex = 0;
      this.canvas.setAttribute('role', 'img');
      this.canvas.setAttribute('aria-label', 'Interactive globe. Drag or use arrow keys to rotate. Scroll or use plus and minus to zoom. Space pauses or resumes rotation. Home resets the view. Select a site from the site list for keyboard access.');
      Object.assign(this.canvas.style, {
        display: 'block', width: '100%', height: '100%', touchAction: 'none', cursor: 'grab',
      });
      this.context = this.canvas.getContext('2d');
      if (!this.context) throw new Error('Canvas is unavailable in this browser.');
      this.tooltip = document.createElement('div');
      this.tooltip.className = 'bio-globe-tooltip';
      this.tooltip.hidden = true;
      Object.assign(this.tooltip.style, {
        position: 'absolute', pointerEvents: 'none', zIndex: '5', maxWidth: '240px',
        padding: '10px 13px', border: '1px solid rgba(118,215,217,.28)',
        borderRadius: '10px', background: 'rgba(7,20,35,.96)', color: '#edf8fb',
        font: '12px/1.5 system-ui, sans-serif', boxShadow: '0 12px 40px #0006',
      });
      this.container.append(this.canvas, this.tooltip);
      this.themeObserver = new MutationObserver(() => this.setTheme(document.documentElement.dataset.theme));
      this.themeObserver.observe(document.documentElement, {attributes: true, attributeFilter: ['data-theme']});
      this.projection = d3.geoOrthographic().clipAngle(90).precision(.35);
      this.path = d3.geoPath(this.projection, this.context);
      this.graticule = d3.geoGraticule().step([20, 20])();
      this._bindEvents();
      this.resizeObserver = new ResizeObserver(() => this.resize());
      this.resizeObserver.observe(container);
      this.resize();
      this.setTheme(this.theme);
      if (this.options.onSpinChange) this.options.onSpinChange(this.spinning);
    }

    setData(clients, server) {
      this.clients = (Array.isArray(clients) ? clients : []).filter(hasLocation);
      this.server = hasLocation(server) ? server : null;
      this.routes = {
        type: 'MultiLineString',
        coordinates: this.server ? this.clients.map(client => [
          [Number(this.server.lng), Number(this.server.lat)],
          [Number(client.lng), Number(client.lat)],
        ]) : [],
      };
      this.hovered = null;
      this.tooltip.hidden = true;
      this.draw();
      this._schedule();
    }

    setTheme(theme) {
      this.theme = theme === 'light' ? 'light' : 'dark';
      const light = this.theme === 'light';
      this.tooltip.style.background = light ? 'rgba(255,255,255,.98)' : 'rgba(7,20,35,.96)';
      this.tooltip.style.color = light ? '#163f4d' : '#edf8fb';
      this.tooltip.style.borderColor = light ? '#b8d6dc' : 'rgba(118,215,217,.28)';
      this.tooltip.style.boxShadow = light ? '0 12px 40px #244e601c' : '0 12px 40px #0006';
      const description = this.tooltip.querySelector('div');
      if (description) description.style.color = light ? '#496b78' : '#a8c0d0';
      this.draw();
    }

    setVisible(visible) {
      this.visible = Boolean(visible);
      if (!this.visible) {
        this.hovered = null;
        this.tooltip.hidden = true;
      }
      this.resize();
      this._schedule();
    }

    focus(lat, lng) {
      if (!hasLocation({lat, lng})) return;
      this.setSpinning(false);
      this.rotation = [-Number(lng), -Number(lat), 0];
      this.hovered = null;
      this.tooltip.hidden = true;
      this.draw();
    }

    setSpinning(spinning) {
      this.spinning = Boolean(spinning);
      if (this.options.onSpinChange) this.options.onSpinChange(this.spinning);
      this._schedule();
    }

    zoom(delta) {
      if (!Number.isFinite(delta)) return;
      this.zoomFactor = clamp(this.zoomFactor * Math.exp(delta), .68, 2.6);
      this.hovered = null;
      this.tooltip.hidden = true;
      this.draw();
    }

    reset() {
      this.rotation = [62, -28, 0];
      this.zoomFactor = 1;
      this.hovered = null;
      this.tooltip.hidden = true;
      this.draw();
    }

    resize() {
      const rect = this.container.getBoundingClientRect();
      if (rect.width < 1 || rect.height < 1) return;
      this.width = rect.width;
      this.height = rect.height;
      const ratio = Math.min(window.devicePixelRatio || 1, 2);
      this.canvas.width = Math.round(this.width * ratio);
      this.canvas.height = Math.round(this.height * ratio);
      this.context.setTransform(ratio, 0, 0, ratio, 0, 0);
      this.draw();
    }

    _bindEvents() {
      const local = event => {
        const rect = this.canvas.getBoundingClientRect();
        return [event.clientX - rect.left, event.clientY - rect.top];
      };
      this.canvas.addEventListener('pointerdown', event => {
        if (event.button !== 0 || this.drag) return;
        this.setSpinning(false);
        this.canvas.focus({preventScroll: true});
        this.canvas.setPointerCapture(event.pointerId);
        this.drag = {id: event.pointerId, origin: local(event), rotation: this.rotation.slice(), moved: false};
        this.hovered = null;
        this.tooltip.hidden = true;
        this.canvas.style.cursor = 'grabbing';
      });
      this.canvas.addEventListener('pointermove', event => {
        const point = local(event);
        if (this.drag && this.drag.id === event.pointerId) {
          const dx = point[0] - this.drag.origin[0];
          const dy = point[1] - this.drag.origin[1];
          if (Math.hypot(dx, dy) > 4) this.drag.moved = true;
          const sensitivity = 70 / Math.max(this.radius || 1, 1);
          this.rotation = [this.drag.rotation[0] + dx * sensitivity,
            clamp(this.drag.rotation[1] - dy * sensitivity, -85, 85), 0];
          this.draw();
        } else if (!this.drag) {
          this.hovered = this._hit(point);
          this.canvas.style.cursor = this.hovered ? 'pointer' : 'grab';
          this._showTooltip(this.hovered, point);
          this.draw();
          this._schedule();
        }
      });
      const finishDrag = (event, cancelled = false) => {
        if (!this.drag || this.drag.id !== event.pointerId) return;
        const clicked = !cancelled && !this.drag.moved;
        this.drag = null;
        if (this.canvas.hasPointerCapture(event.pointerId)) this.canvas.releasePointerCapture(event.pointerId);
        this.canvas.style.cursor = 'grab';
        if (clicked) {
          const hit = this._hit(local(event));
          if (hit && !hit.server && this.options.onSelect) this.options.onSelect(hit.item.client_id);
        }
        this._schedule();
      };
      this.canvas.addEventListener('pointerup', event => finishDrag(event));
      this.canvas.addEventListener('pointercancel', event => finishDrag(event, true));
      this.canvas.addEventListener('lostpointercapture', event => finishDrag(event, true));
      this.canvas.addEventListener('pointerleave', () => {
        this.hovered = null;
        this.tooltip.hidden = true;
        this.draw();
        this._schedule();
      });
      this.canvas.addEventListener('wheel', event => {
        event.preventDefault();
        this.setSpinning(false);
        const units = event.deltaMode === 1 ? 16 : event.deltaMode === 2 ? this.height : 1;
        this.zoom(-clamp(event.deltaY * units, -160, 160) * .002);
      }, {passive: false});
      this.canvas.addEventListener('keydown', event => {
        const moves = {ArrowLeft: [-7, 0], ArrowRight: [7, 0], ArrowUp: [0, 7], ArrowDown: [0, -7]};
        if (moves[event.key]) {
          event.preventDefault();
          this.setSpinning(false);
          this.rotation[0] += moves[event.key][0];
          this.rotation[1] = clamp(this.rotation[1] + moves[event.key][1], -85, 85);
          this.draw();
        } else if (['+', '=', '-', '_', ' ', 'Home'].includes(event.key)) {
          event.preventDefault();
          if (event.key === ' ') this.setSpinning(!this.spinning);
          else if (event.key === 'Home') this.reset();
          else this.zoom(event.key === '-' || event.key === '_' ? -.18 : .18);
        }
      });
      this._visibilityChange = () => {
        if (!document.hidden) this.draw();
        this._schedule();
      };
      document.addEventListener('visibilitychange', this._visibilityChange);
      this._motionChange = event => { if (event.matches) this.setSpinning(false); };
      this.motionPreference.addEventListener('change', this._motionChange);
    }

    _schedule() {
      const animate = this.visible && !document.hidden && this.spinning && !this.drag && !this.hovered;
      if (!animate) {
        if (this.frame !== null) cancelAnimationFrame(this.frame);
        this.frame = null;
        this.lastFrame = null;
        return;
      }
      if (this.frame !== null) return;
      this.frame = requestAnimationFrame(time => {
        this.frame = null;
        const elapsed = this.lastFrame === null ? 0 : Math.min(time - this.lastFrame, 64);
        this.lastFrame = time;
        this.rotation[0] = (this.rotation[0] + elapsed * .0035) % 360;
        this.draw();
        this._schedule();
      });
    }

    _hit(point) {
      let closest = null;
      let distance = 15;
      for (const target of this.hitTargets) {
        const current = Math.hypot(point[0] - target.x, point[1] - target.y);
        if (current < distance) { closest = target; distance = current; }
      }
      return closest;
    }

    _showTooltip(target, point) {
      this.tooltip.hidden = !target;
      if (!target) return;
      const title = document.createElement('strong');
      title.textContent = target.item.institution || target.item.client_id || 'Coordinator';
      const description = document.createElement('div');
      description.style.color = this.theme === 'light' ? '#496b78' : '#a8c0d0';
      description.textContent = target.server ? 'Coordinator' :
        `${target.item.partnership_stage || target.item.status || 'Site'} · Select to inspect`;
      this.tooltip.replaceChildren(title, description);
      const left = clamp(point[0] + 17, 8, Math.max(8, this.width - this.tooltip.offsetWidth - 12));
      const top = clamp(point[1] - 16, 8, Math.max(8, this.height - this.tooltip.offsetHeight - 12));
      this.tooltip.style.left = `${left}px`;
      this.tooltip.style.top = `${top}px`;
    }

    draw() {
      if (!this.visible || !this.width || !this.height || document.hidden) return;
      const ctx = this.context;
      const light = this.theme === 'light';
      const w = this.width, h = this.height;
      const cx = w / 2, cy = h / 2;
      const radius = Math.min(w, h) * .405 * this.zoomFactor;
      this.radius = radius;
      this.projection.translate([cx, cy]).scale(radius).rotate(this.rotation);
      ctx.clearRect(0, 0, w, h);

      const sky = ctx.createRadialGradient(cx, cy, 0, cx, cy, Math.max(w, h) * .7);
      sky.addColorStop(0, light ? '#ffffff' : '#102639');
      sky.addColorStop(.6, light ? '#f0f7f9' : '#081726');
      sky.addColorStop(1, light ? '#e2edf2' : '#060f1b');
      ctx.fillStyle = sky;
      ctx.fillRect(0, 0, w, h);
      // Fixed deterministic stars remain still as the earth rotates.
      for (let i = 1; i <= 95; i++) {
        const x = ((i * 7919) % 997) / 997 * w;
        const y = ((i * 3571) % 991) / 991 * h;
        ctx.fillStyle = light ? 'rgba(82,139,160,.13)' : `rgba(171,210,228,${.11 + (i % 4) * .06})`;
        ctx.beginPath();
        ctx.arc(x, y, i % 7 === 0 ? 1 : .6, 0, TAU);
        ctx.fill();
      }

      const atmosphere = ctx.createRadialGradient(cx, cy, radius * .985, cx, cy, radius * 1.075);
      atmosphere.addColorStop(0, light ? 'rgba(56,165,188,.24)' : 'rgba(90,219,232,.29)');
      atmosphere.addColorStop(.28, 'rgba(58,174,207,.10)');
      atmosphere.addColorStop(1, 'rgba(35,123,167,0)');
      ctx.fillStyle = atmosphere;
      ctx.beginPath();
      ctx.arc(cx, cy, radius * 1.075, 0, TAU);
      ctx.fill();

      ctx.save();
      ctx.beginPath();
      ctx.arc(cx, cy, radius, 0, TAU);
      ctx.clip();
      const ocean = ctx.createRadialGradient(cx - radius * .36, cy - radius * .48, radius * .1, cx, cy, radius * 1.2);
      ocean.addColorStop(0, light ? '#f6fcff' : '#163c52');
      ocean.addColorStop(.65, light ? '#d8edf5' : '#0c273e');
      ocean.addColorStop(1, light ? '#bbdce7' : '#061321');
      ctx.fillStyle = ocean;
      ctx.fillRect(cx - radius, cy - radius, radius * 2, radius * 2);

      ctx.beginPath();
      this.path(this.graticule);
      ctx.strokeStyle = light ? 'rgba(63,119,141,.17)' : 'rgba(138,193,211,.115)';
      ctx.lineWidth = .65;
      ctx.stroke();

      const land = ctx.createLinearGradient(cx - radius, cy - radius, cx + radius, cy + radius);
      land.addColorStop(0, light ? '#bee4d9' : '#356b73');
      land.addColorStop(.5, light ? '#93c9bf' : '#25565f');
      land.addColorStop(1, light ? '#6fa99f' : '#163d4e');
      ctx.beginPath();
      this.path(window.BIO_LAND);
      ctx.fillStyle = land;
      ctx.fill();
      ctx.strokeStyle = light ? 'rgba(51,121,127,.5)' : 'rgba(121,205,206,.38)';
      ctx.lineWidth = .65;
      ctx.stroke();

      // Limb shading gives the orthographic projection a spherical appearance.
      const shade = ctx.createRadialGradient(cx - radius * .24, cy - radius * .2, radius * .25, cx, cy, radius);
      shade.addColorStop(0, 'rgba(0,6,17,0)');
      shade.addColorStop(.72, light ? 'rgba(23,68,87,.01)' : 'rgba(0,6,17,.03)');
      shade.addColorStop(1, light ? 'rgba(23,68,87,.16)' : 'rgba(0,6,17,.59)');
      ctx.fillStyle = shade;
      ctx.fillRect(cx - radius, cy - radius, radius * 2, radius * 2);

      if (this.routes && this.routes.coordinates.length) {
        ctx.beginPath();
        this.path(this.routes);
        ctx.strokeStyle = light ? 'rgba(0,125,131,.55)' : 'rgba(99,224,220,.32)';
        ctx.lineWidth = 1.15;
        ctx.stroke();
      }
      ctx.restore();
      ctx.beginPath();
      ctx.arc(cx, cy, radius, 0, TAU);
      ctx.strokeStyle = light ? 'rgba(49,132,155,.45)' : 'rgba(131,229,235,.4)';
      ctx.lineWidth = 1;
      ctx.stroke();

      const center = this.projection.invert([cx, cy]);
      this.hitTargets = [];
      const visiblePoint = item => {
        const location = [Number(item.lng), Number(item.lat)];
        return d3.geoDistance(center, location) < Math.PI / 2 - .008 ? this.projection(location) : null;
      };
      const coordinator = this.server && visiblePoint(this.server);
      const points = this.clients.map(item => ({item, point: visiblePoint(item)})).filter(entry => entry.point);
      for (const entry of points) {
        const [x, y] = entry.point;
        this._drawMarker(x, y, entry.item, false);
      }
      if (coordinator) {
        const [x, y] = coordinator;
        // Separate a co-located coordinator from its institution's selectable marker.
        const overlaps = points.some(entry => Math.hypot(entry.point[0] - x, entry.point[1] - y) < 17);
        if (overlaps) {
          ctx.beginPath();
          ctx.moveTo(x + 3, y - 3);
          ctx.lineTo(x + 15, y - 15);
          ctx.strokeStyle = 'rgba(246,199,120,.6)';
          ctx.lineWidth = 1;
          ctx.stroke();
        }
        this._drawMarker(x + (overlaps ? 19 : 0), y - (overlaps ? 19 : 0), this.server, true);
      }
    }

    _drawMarker(x, y, item, server) {
      const ctx = this.context;
      const stage = String(item.partnership_stage || '').slice(0, 1);
      const stageColors = { '0': '#8d9eaa', '1': '#669fe0', '2': '#d89e42', '4': '#15b9a6', X: '#db717e' };
      const color = server ? (this.theme === 'light' ? '#bd7a17' : '#f5c77e')
        : stageColors[stage] || (this.theme === 'light' && item.status === 'active' ? '#088981' : COLORS[item.status] || '#66ddd7');
      const hovered = this.hovered && this.hovered.server === server &&
        (server || this.hovered.item.client_id === item.client_id);
      ctx.save();
      ctx.shadowColor = color;
      ctx.shadowBlur = hovered ? 20 : 12;
      ctx.beginPath();
      ctx.arc(x, y, server ? 9 : hovered ? 10 : 8, 0, TAU);
      ctx.fillStyle = server ? 'rgba(245,199,126,.12)' : 'rgba(84,221,217,.12)';
      ctx.fill();
      ctx.shadowBlur = 0;
      ctx.strokeStyle = color;
      ctx.globalAlpha = hovered ? .85 : .4;
      ctx.lineWidth = 1;
      ctx.stroke();
      ctx.globalAlpha = 1;
      ctx.beginPath();
      if (server) {
        ctx.moveTo(x, y - 5);
        ctx.lineTo(x + 5, y);
        ctx.lineTo(x, y + 5);
        ctx.lineTo(x - 5, y);
        ctx.closePath();
      } else ctx.arc(x, y, hovered ? 4.5 : 3.4, 0, TAU);
      ctx.fillStyle = color;
      ctx.fill();
      ctx.strokeStyle = '#e5ffff';
      ctx.globalAlpha = .75;
      ctx.lineWidth = .75;
      ctx.stroke();
      ctx.restore();
      this.hitTargets.push({x, y, item, server});
    }

    destroy() {
      this.setVisible(false);
      this.resizeObserver.disconnect();
      this.themeObserver.disconnect();
      document.removeEventListener('visibilitychange', this._visibilityChange);
      this.motionPreference.removeEventListener('change', this._motionChange);
      this.canvas.remove();
      this.tooltip.remove();
    }
  }

  window.BioGlobe = BioGlobe;
})();
