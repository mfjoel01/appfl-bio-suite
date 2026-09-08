"""TRS and TES: the tool pin that makes a result attributable, and the task that runs it.

These two are tested together because they are one mechanism: TRS says which container
image and which descriptor, TES puts exactly those into a task document. A test that
checked either alone would miss the join, which is the part that can silently be wrong --
a task that runs an image no pin names is a task whose results nobody can attribute.

The TES service itself is stubbed with a small HTTP server implementing the parts of the
API this suite calls. That is deliberate rather than a compromise: what needs testing is
that this client speaks TES correctly, and standing up Funnel to learn that would test
Funnel.
"""

from __future__ import annotations

import json
import threading
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer

import pytest

from appfl_bio_suite.core.ga4gh.tes import (
    TERMINAL_STATES,
    TesClient,
    TesError,
    build_site_task,
    write_task,
)
from appfl_bio_suite.core.ga4gh.trs import (
    TrsClient,
    TrsError,
    build_tool,
    default_tool_id,
    descriptor_checksum,
    verify_pin,
    write_registry,
)
from appfl_bio_suite.core.ga4gh.trs import serve as serve_trs

# ---------------------------------------------------------------------------
# TRS
# ---------------------------------------------------------------------------


def test_the_tool_describes_the_installed_site_stage():
    tool, files = build_tool()
    assert tool.versions and tool.versions[0].descriptor_type == ["CWL"]
    # The two modules APPFL actually ships. Their checksums are what makes the pin cover
    # the computation rather than a wrapper around it.
    assert "src/fine_mapping/dataset.py" in files
    assert "src/fine_mapping/trainer.py" in files
    assert "fine-mapping-site-stage.cwl" in files
    assert "Containerfile" in files


def test_the_descriptor_wraps_the_real_command():
    _, files = build_tool()
    cwl = files["fine-mapping-site-stage.cwl"]
    assert "baseCommand: [appfl-bio-suite, site-stage]" in cwl
    assert "--data-dir" in cwl and "--config" in cwl


def test_the_shipped_sources_in_the_tool_are_the_installed_ones():
    from appfl_bio_suite.core.experiments import get_spec

    _, files = build_tool()
    spec = get_spec("fine-mapping")
    for module in spec.shipped_modules:
        on_disk = (spec.package_path / f"{module}.py").read_text(encoding="utf-8")
        assert files[f"src/{spec.package}/{module}.py"] == on_disk


def test_a_pin_verifies_against_the_install(tmp_path):
    _, pin = write_registry(tmp_path)
    assert verify_pin(pin) == []


def test_a_tampered_pin_fails_and_says_what_to_do(tmp_path):
    _, pin = write_registry(tmp_path)
    bad = pin.model_copy(update={"descriptor_checksum": "0" * 64})
    problems = verify_pin(bad)
    assert problems and "does not match its pin" in problems[0]
    assert "re-publish and re-pin" in problems[0] or "changed the code" in problems[0]


def test_the_checksum_covers_paths_as_well_as_content():
    """Otherwise a file added and another removed can leave the digest unchanged."""
    base = {"a.cwl": "x", "b.txt": "y"}
    renamed = {"a.cwl": "x", "c.txt": "y"}
    assert descriptor_checksum(base) != descriptor_checksum(renamed)


def test_editing_a_shipped_module_would_change_the_pin(tmp_path):
    """The property the whole pin rests on, stated as a test rather than as a comment."""
    _, files = build_tool()
    edited = dict(files)
    edited["src/fine_mapping/trainer.py"] += "\n# an innocuous comment\n"
    assert descriptor_checksum(files) != descriptor_checksum(edited)


def test_the_tool_id_is_derived_from_the_repository_not_hardcoded():
    """A fork's runs must not be attributed to the upstream repository."""
    from appfl_bio_suite import REPO_URL
    from appfl_bio_suite.core.ga4gh.trs import repository_slug

    assert default_tool_id("fine-mapping", repo="github.com/someone/their-fork").startswith(
        "#workflow/github.com/someone/their-fork"
    )
    # And the default follows the one constant that already names this install's repo.
    assert repository_slug() in REPO_URL
    expected = f"#workflow/{repository_slug()}/fine-mapping-site-stage"
    assert default_tool_id("fine-mapping") == expected


def test_the_generated_tree_is_servable_and_conformant(tmp_path):
    _, pin = write_registry(
        tmp_path, image="registry.example.org/site-stage:0.1.0",
        registry_url="https://trs.example.org",
    )
    httpd = serve_trs(tmp_path, "127.0.0.1", 0)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    base = f"http://127.0.0.1:{httpd.server_address[1]}"
    try:
        client = TrsClient(base)
        assert client.service_info()["type"]["artifact"] == "trs"
        tool = client.tool(pin.id)
        assert tool.id == pin.id
        version = client.version(pin.id, pin.version)
        assert version.id == pin.version
        descriptor = client.descriptor(pin.id, pin.version)
        assert "cwlVersion" in (descriptor.content or "")
        paths = {f.path for f in client.files(pin.id, pin.version)}
        assert "fine-mapping-site-stage.cwl" in paths
        assert client.image_for(pin.id, pin.version) == "registry.example.org/site-stage:0.1.0"
    finally:
        httpd.shutdown()


def test_an_image_digest_is_pinned_as_a_digest(tmp_path):
    """A tag can be repushed; a digest is the only durable image identity."""
    digest = "a" * 64
    write_registry(
        tmp_path, image="registry.example.org/y:1.0", image_digest=digest, registry_url="http://x"
    )
    tool, _ = build_tool(image="registry.example.org/y:1.0", image_digest=digest)
    assert tool.versions[0].images[0].digest == digest


def test_publishing_without_an_image_registers_none(tmp_path):
    """No plausible-looking image is invented. A tool version that names an image nobody
    pushed fails at the executor, minutes into a run, instead of at the pin."""
    tool, files = build_tool()
    assert tool.versions[0].images == []
    assert "IMAGE-NOT-PUBLISHED" in files["fine-mapping-site-stage.cwl"]
    _, pin = write_registry(tmp_path)
    assert pin.image is None


def test_a_version_without_an_image_says_so_rather_than_guessing(tmp_path):
    """A TES executor cannot invent an image, so resolving one must fail loudly."""
    _, pin = write_registry(tmp_path)
    httpd = serve_trs(tmp_path, "127.0.0.1", 0)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    try:
        client = TrsClient(f"http://127.0.0.1:{httpd.server_address[1]}")
        with pytest.raises(TrsError, match="nothing for a TES executor to run"):
            client.image_for(pin.id, pin.version)
    finally:
        httpd.shutdown()


def test_dockstore_registration_is_emitted(tmp_path):
    write_registry(tmp_path)
    dockstore = (tmp_path / ".dockstore.yml").read_text(encoding="utf-8")
    assert "subclass: CWL" in dockstore
    assert "primaryDescriptorPath" in dockstore


# ---------------------------------------------------------------------------
# TES
# ---------------------------------------------------------------------------


def a_task(**overrides):
    kwargs = {
        "client_id": "anl",
        "experiment": "fine-mapping",
        "image": "ghcr.io/x/y@sha256:" + "a" * 64,
        "bundle_url": "drs://drs.test.org/" + "b" * 64,
        "outputs_url": "file:///scratch/out/anl",
        "run_config": {"client_id": "anl", "train_configs": {"locus_n_shards": 4}},
    }
    kwargs.update(overrides)
    return build_site_task(**kwargs)


def test_the_task_stages_the_bundle_by_drs_uri():
    task = a_task()
    bundle = next(i for i in task.inputs if i.name == "bundle")
    assert bundle.url.startswith("drs://")
    assert bundle.type == "DIRECTORY"


def test_the_run_config_is_inlined_not_staged():
    """A few hundred bytes that differ per site; an object store entry would be one more
    thing to keep in step with federation.yaml."""
    task = a_task()
    config = next(i for i in task.inputs if i.name == "run-config")
    assert config.url is None
    assert json.loads(config.content)["train_configs"]["locus_n_shards"] == 4


def test_the_executor_runs_the_trs_pinned_image_and_the_real_command():
    task = a_task()
    executor = task.executors[0]
    assert executor.image.startswith("ghcr.io/x/y@sha256:")
    assert executor.command[:3] == ["appfl-bio-suite", "site-stage", "fine-mapping"]
    # The same thread caps the Globus Compute workers get.
    assert executor.env["OMP_NUM_THREADS"] == "1"


def test_provenance_travels_in_the_tags():
    task = a_task(tags={"trs_id": "#workflow/x", "duo": "permitted"})
    assert task.tags["client_id"] == "anl"
    assert task.tags["trs_id"] == "#workflow/x"


def test_server_assigned_fields_are_stripped_from_a_create():
    task = a_task()
    task.id = "assigned-by-someone-else"
    task.state = "COMPLETE"
    body = task.request_body()
    assert "id" not in body and "state" not in body


def test_write_task_writes_the_submittable_document(tmp_path):
    path = write_task(a_task(), tmp_path / "anl.task.json")
    body = json.loads(path.read_text(encoding="utf-8"))
    assert body["executors"][0]["command"][0] == "appfl-bio-suite"
    assert "id" not in body


# -- the client, against a stub service --------------------------------------


class _StubTES(BaseHTTPRequestHandler):
    """The parts of TES 1.1 this client calls. States advance on each poll."""

    tasks: dict = {}
    sequence = ["QUEUED", "INITIALIZING", "RUNNING", "COMPLETE"]

    def _send(self, code, payload):
        body = json.dumps(payload).encode("utf-8")
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self):  # noqa: N802
        path = self.path.split("?")[0]
        if path.endswith("/service-info"):
            self._send(200, {"id": "stub", "name": "stub TES",
                             "type": {"group": "org.ga4gh", "artifact": "tes", "version": "1.1.0"}})
            return
        task_id = path.rsplit("/", 1)[-1]
        task = _StubTES.tasks.get(task_id)
        if task is None:
            self._send(404, {"message": "no such task"})
            return
        index = min(task["polls"], len(_StubTES.sequence) - 1)
        task["polls"] += 1
        self._send(200, {"id": task_id, "state": _StubTES.sequence[index],
                         "executors": task["executors"], "logs": []})

    def do_POST(self):  # noqa: N802
        length = int(self.headers.get("Content-Length", 0))
        body = json.loads(self.rfile.read(length) or b"{}")
        if self.path.endswith(":cancel"):
            self._send(200, {})
            return
        task_id = f"task-{len(_StubTES.tasks) + 1}"
        _StubTES.tasks[task_id] = {"polls": 0, "executors": body.get("executors", [])}
        self._send(200, {"id": task_id})

    def log_message(self, *args):
        return


@pytest.fixture
def tes_service():
    _StubTES.tasks = {}
    httpd = ThreadingHTTPServer(("127.0.0.1", 0), _StubTES)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    yield f"http://127.0.0.1:{httpd.server_address[1]}"
    httpd.shutdown()


def test_service_info_round_trip(tes_service):
    info = TesClient(tes_service).service_info()
    assert info["type"]["artifact"] == "tes"


def test_submit_and_wait_to_completion(tes_service):
    client = TesClient(tes_service)
    task_id = client.create_task(a_task())
    seen = []
    finished = client.wait(task_id, poll_seconds=0, on_state=lambda _i, s: seen.append(s))
    assert finished.state == "COMPLETE"
    assert "RUNNING" in seen
    assert finished.state in TERMINAL_STATES


def test_a_wait_that_times_out_cancels_the_task(tes_service, monkeypatch):
    """A driver that gave up while a site kept computing would burn its allocation."""
    monkeypatch.setattr(_StubTES, "sequence", ["RUNNING"])
    client = TesClient(tes_service)
    task_id = client.create_task(a_task())
    with pytest.raises(TesError, match="still RUNNING"):
        client.wait(task_id, poll_seconds=0, timeout_seconds=-1)


def test_an_unreachable_service_names_the_setting():
    with pytest.raises(TesError, match="ga4gh.tes.url"):
        TesClient("http://127.0.0.1:1").service_info()


def test_a_missing_task_is_an_error_carrying_the_body(tes_service):
    with pytest.raises(TesError, match="404"):
        TesClient(tes_service).get_task("nope")


# -- the payload wire format the TES path depends on -------------------------


def test_the_npz_payload_round_trips_exactly(tmp_path):
    """The TES path's wire format. Two things can break here and nothing else notices.

    The trainer's payload keys carry structure -- ``G::L0000::EUR`` -- and ``np.savez``
    uses them as member names inside a zip. If they were mangled, the aggregator's
    decoder would look for blocks that are present under other names and report a site
    that "returned no genotype blocks". And the manifest is bytes, not an array: it has
    to come back as the same uint8 tensor the APPFL path delivers, or the decoder would
    have to know which driver produced its input.
    """
    import numpy as np
    import torch

    from appfl_bio_suite.experiments.fine_mapping.site_stage import load_payload
    from appfl_bio_suite.experiments.fine_mapping.trainer import KEY_MANIFEST

    arrays = {
        "G::L0000::EUR": np.arange(9, dtype=np.float64).reshape(3, 3),
        "u::L0000::EUR": np.array([1.0, 2.0, 3.0]),
        "c::L0000_arch_rep0::EUR": np.array([0.5, -0.25, 0.125]),
        "s::L0000_arch_rep0::EUR": np.array([10.0, 20.0, 30.0]),
        # float32 is what a site sends when the Gram is provably exact in it.
        "G::L0001::AFR": np.arange(4, dtype=np.float32).reshape(2, 2),
    }
    manifest = {
        "client_id": "anl",
        "geno_blocks": [],
        "ga4gh": {"data_use": {"outcome": "permitted"}},
    }

    np.savez(tmp_path / "aggregates.npz", **arrays)
    (tmp_path / "aggregates.manifest.json").write_text(json.dumps(manifest), encoding="utf-8")

    payload = load_payload(tmp_path / "aggregates.npz")

    assert set(payload) == set(arrays) | {KEY_MANIFEST}
    for key, expected in arrays.items():
        assert torch.is_tensor(payload[key])
        np.testing.assert_array_equal(payload[key].numpy(), expected)
        assert payload[key].numpy().dtype == expected.dtype

    blob = payload[KEY_MANIFEST]
    assert blob.dtype is torch.uint8
    assert json.loads(bytes(blob.numpy().tobytes()).decode("utf-8")) == manifest
