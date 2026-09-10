/* Run against any freshly exported map.html. See README.md for browser setup. */
const assert = require('node:assert/strict');
const fs = require('node:fs');
const { chromium } = require(process.env.PLAYWRIGHT_MODULE || 'playwright');
const viewerUrl = process.env.WATCH_VIEWER_URL || 'http://127.0.0.1:7081/current/map.html';
const site = (id, experiments, lat, lng) => ({client_id:id, institution:`Institute ${id}`,
  experiments, lat, lng, country:'Example country', city:'Example city', num_samples:100, status:'active'});
const sites = [site('north','fine-mapping, gwas',40,-75), site('south','fine-mapping',-20,15),
  site('east','gwas',25,55), site('unplaced','new-experiment',null,null), site('zero','gwas',0,0)];
const network = {run_id:'network',algorithm:'federation network',finished_at:'2026-01-01',
  server:{lat:40,lng:-75,org:'Coordinator'},rounds:[{round:1,clients:sites}]};
const live = {run_id:'running',algorithm:'gwas',finished_at:null,server:network.server,
  rounds:[{round:1,clients:[{...sites[0],experiments:undefined,experiment:'gwas'}]}]};
const events = [
 {event_type:'init',run_id:'recording',algorithm:'gwas'},
 {event_type:'server_metadata',run_id:'recording',server:network.server},
 {event_type:'client_update',run_id:'recording',round:1,clients:[sites[0]]},
 {event_type:'round_end',run_id:'recording',round:1,clients:[{client_id:'north',lat:null,lng:null,experiments:null}],round_metrics:{global_accuracy:.7}},
 {event_type:'comm_failure',run_id:'recording',round:1,client_id:'north'},
 {event_type:'round_end',run_id:'recording',round:2,clients:[sites[1]],round_metrics:{global_accuracy:.8}},
 {event_type:'finished',run_id:'recording'}
];
(async () => {
 const browser = await chromium.launch({headless:true,
   ...(process.env.CHROMIUM_EXECUTABLE ? {executablePath:process.env.CHROMIUM_EXECUTABLE} : {})});
 const page = await browser.newPage({viewport:{width:1440,height:1000}});
 const errors=[]; page.on('pageerror',e=>errors.push(e.message));
 // Exercise real viewer code and loaders with deterministic public fixture data.
 await page.addInitScript(() => {
   window.testStreams=[];
   window.EventSource=class {
     constructor(url){this.url=url;testStreams.push(this);setTimeout(()=>this.onopen?.(),0)}
     close(){}
     emit(data){this.onmessage?.({data:JSON.stringify(data)})}
   };
 });
 await page.route('**/network.map.json',route=>route.fulfill({json:network}));
 await page.route('**/fixture-events.json',route=>route.fulfill({json:events}));
 await page.route('**/runs',route=>route.fulfill({json:[
   {run_id:'network',algorithm:'federation network'},{run_id:'running',algorithm:'gwas'},
   {run_id:'recording',algorithm:'gwas'}]}));
 await page.route('**/runs/*/metadata',route=> {
   const id=route.request().url().split('/').at(-2);
   return id==='recording' ? route.fulfill({status:404,body:'absent'}) : route.fulfill({json:id==='network'?network:live});
 });
 await page.route('**/runs/*/events',route=>route.fulfill({json:events}));
 const visible=()=>page.evaluate(()=>bioWatch.getState().visibleSites);
 const checkbox=name=>page.locator(`input[data-experiment="${name}"]`);
 const assertVisible=async ids=>{await page.waitForFunction(expected=>JSON.stringify(bioWatch.getState().visibleSites)===JSON.stringify(expected),ids);
   assert.deepEqual(await visible(),ids);
   assert.equal(await page.locator('.bio-site-card').count(),ids.length);
   assert.equal(await page.locator('#hdr-clients').innerText(),String(ids.length));};
 await page.goto(viewerUrl+'?metadata_url=network.map.json');
 await assertVisible(sites.map(s=>s.client_id));
 assert.equal(await page.locator('#map').isVisible(),true);
 assert.equal(await page.locator('#bio-runs-panel').isVisible(),false);
 for(const id of ['hdr-round','hdr-acc','hdr-loss']) assert.equal(await page.locator('#'+id).isVisible(),false);
 assert.equal(await checkbox('*').isChecked(),true);
 assert.deepEqual(await page.evaluate(()=>Object.keys(markers).filter(id=>map.hasLayer(markers[id]))),['north','south','east','zero']);
 await checkbox('fine-mapping').check(); await assertVisible(['north','south']);
 assert.deepEqual(await page.evaluate(()=>Object.keys(lines).filter(id=>map.hasLayer(lines[id]))),['north','south']);
 await checkbox('gwas').check(); await assertVisible(['north','south','east','zero']);
 await checkbox('fine-mapping').uncheck(); await assertVisible(['north','east','zero']);
 await page.evaluate(()=>{const original=BioGlobe.prototype.setData;BioGlobe.prototype.setData=function(...args){window.testGlobe=this;return original.apply(this,args)}});
 await page.locator('#bio-globe').click();
 await page.waitForFunction(()=>window.testGlobe?.visible);
 assert.deepEqual(await page.evaluate(()=>testGlobe.clients.map(c=>c.client_id)),['north','east','zero']);
 assert.equal(await page.locator('#map').isVisible(),false);
 await checkbox('gwas').uncheck();await assertVisible([]);
 assert.equal(await page.evaluate(()=>testGlobe.clients.length),0);
 assert.equal(await page.locator('#bio-map-message').isVisible(),true);
 await checkbox('*').check(); await assertVisible(sites.map(s=>s.client_id));
 assert.equal(await page.evaluate(()=>testGlobe.clients.length),4);
 await page.locator('[data-site="north"]').click();
 assert.equal(await page.locator('#bio-site-detail').isVisible(),true);
 assert.equal(await page.evaluate(()=>testGlobe.spinning),false);
 assert.deepEqual(await page.evaluate(()=>testGlobe.rotation),[75,-40,0]);
 await page.keyboard.press('Escape');
 assert.equal(await page.locator('#bio-site-detail').isVisible(),false);
 await page.locator('#bio-spin').click();
 const rotation=await page.evaluate(()=>testGlobe.rotation[0]);await page.waitForTimeout(150);
 assert.notEqual(await page.evaluate(()=>testGlobe.rotation[0]),rotation);
 await page.locator('#bio-flat').click();
 assert.equal(await page.evaluate(()=>testGlobe.frame),null);
 await page.locator('[data-site="zero"]').click();
 assert.equal(await page.evaluate(()=>markers.zero.isPopupOpen()),true);
 await page.locator('[data-site="unplaced"]').click();
 assert.equal(await page.locator('#bio-site-detail').isVisible(),true);
 assert.equal(await page.evaluate(()=>markers.unplaced===undefined),true);
 await page.keyboard.press('Escape');
 await page.locator('#bio-globe').click();
 const canvas=page.locator('canvas');const box=await canvas.boundingBox();
 await page.mouse.move(box.x+box.width*.5,box.y+box.height*.5);
 await page.mouse.down();await page.mouse.move(box.x+box.width*.6,box.y+box.height*.55);await page.mouse.up();
 assert.equal(await page.evaluate(()=>testGlobe.spinning),false);
 const zoom=await page.evaluate(()=>testGlobe.zoomFactor);
 await page.locator('#bio-zoom-in').click();assert.ok(await page.evaluate(()=>testGlobe.zoomFactor)>zoom);
 await canvas.focus();await page.keyboard.press('ArrowRight');await page.keyboard.press('Home');
 assert.deepEqual(await page.evaluate(()=>testGlobe.rotation),[62,-28,0]);
 await page.emulateMedia({reducedMotion:'reduce'});
 await page.reload();await assertVisible(sites.map(s=>s.client_id));
 await page.locator('#bio-globe').click();assert.equal(await page.locator('#bio-spin').getAttribute('aria-pressed'),'false');
 // Small screens retain reachable controls and have no horizontal overflow.
 await page.setViewportSize({width:390,height:844});
 assert.equal(await page.evaluate(()=>document.documentElement.scrollWidth<=innerWidth),true);
 assert.equal(await page.locator('#bio-globe').isVisible(),true);
 await page.locator('#bio-runs-tab').click();
 assert.equal(await page.locator('.playback-bar').isVisible(),true);
 assert.equal(await page.locator('#runs-list').isVisible(),true);
 if(process.env.WATCH_SCREENSHOT_DIR){fs.mkdirSync(process.env.WATCH_SCREENSHOT_DIR,{recursive:true});await page.screenshot({path:process.env.WATCH_SCREENSHOT_DIR+'/mobile.png'})}
 await page.setViewportSize({width:1440,height:1000});
 // JSONL fallback renders immediately, retains round_end clients and replays all frames.
 await page.goto(viewerUrl+'?events_url=fixture-events.json');
 await assertVisible(['north','south']);
 assert.equal(await page.evaluate(()=>state.clients.north.lat),40);
 assert.equal(await page.evaluate(()=>state.clients.north.status),'failed');
 await page.locator('#bio-runs-tab').click();await page.locator('#play-btn').click();
 await page.waitForFunction(()=>pb.index===1);assert.equal(await page.evaluate(()=>state.round),1);
 await page.waitForFunction(()=>pb.index===2);assert.equal(await page.evaluate(()=>state.round),2);
 // An unfinished run selected after its init was missed still receives updates.
 await page.goto(viewerUrl);await assertVisible(sites.map(s=>s.client_id));
 await page.locator('#bio-runs-tab').click();await page.locator('[data-run="running"]').click();
 await page.waitForFunction(()=>state.mode==='live');
 assert.deepEqual(await page.evaluate(()=>bioWatch.getState().experiments),['gwas']);
 await page.evaluate(()=>testStreams[0].emit({event_type:'client_update',run_id:'running',round:2,
   clients:[{client_id:'north',num_samples:222}]}));
 assert.equal(await page.evaluate(()=>state.clients.north.num_samples),222);
 assert.equal(await page.evaluate(()=>pb.rounds.length),2);
 await page.evaluate(()=>testStreams[0].emit({event_type:'comm_failure',run_id:'running',round:2,client_id:'north'}));
 assert.equal(await page.evaluate(()=>state.clients.north.status),'failed');
 await page.evaluate(()=>testStreams[0].emit({event_type:'finished',run_id:'running'}));
 assert.equal(await page.evaluate(()=>state.mode),'replay');
 assert.equal(await page.locator('#play-btn').isEnabled(),true);
 await page.locator('#bio-network-return').click();await assertVisible(sites.map(s=>s.client_id));
 // Unrelated live runs must never clear the network or change a chosen filter.
 await checkbox('fine-mapping').check();
 await page.evaluate(()=>{testStreams[0].emit({event_type:'init',run_id:'unrelated',mode:'live'});
   testStreams[0].emit({event_type:'round_end',run_id:'unrelated',round:10,clients:[{client_id:'intruder',lat:0,lng:0}]})});
 await assertVisible(['north','south']);
 assert.equal(await page.evaluate(()=>bioWatch.getState().runId),'network');
 // Log/packet/card/popup strings remain text, including an already-open popup.
 await page.evaluate(()=>{
   const id='<img src=x onerror="window.injected=true">';
   upsertClient(id,{client_id:id,institution:id,lat:5,lng:5,experiments:'gwas'});
   markers[id].openPopup();upsertClient(id,{city:id});addLog(id);
   document.getElementById('bio-site-detail').innerHTML=buildPacketTooltip(id,'uplink');
 });
 await page.waitForTimeout(100);assert.equal(await page.evaluate(()=>window.injected),undefined);
 assert.deepEqual(errors,[]);
 console.log('PASS: default flat/all, OR/none filters, markers/arcs/counts, globe interactions, null/zero coordinates, mobile, reduced motion, JSONL replay, late live selection, SSE isolation, escaped metadata');
 await browser.close();
})().catch(error=>{console.error(error);process.exit(1)});
