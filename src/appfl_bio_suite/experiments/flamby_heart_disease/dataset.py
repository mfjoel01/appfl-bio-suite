"""FLamby dataset loader. SHIPPED TO WORKERS -- read the import rule before editing.

=============================================================================
THIS FILE'S SOURCE IS SENT OVER THE WIRE AND EXECUTED ON A PARTNER'S CLUSTER.
=============================================================================

APPFL's Globus Compute communicator reads this file on the coordinating driver, inlines
its text, and ships it to the worker::

    with open(client_config.data_configs.dataset_path) as file:
        client_config.data_configs.dataset_source = file.read()
    del client_config.data_configs.dataset_path

Two consequences, both load-bearing:

1. **It may only import what the partner has installed.** In particular it must NOT
   import from ``appfl_bio_suite`` -- a partner installs the suite only if we tell them
   to, and the point of keeping this self-contained is that they need not. Standard
   library plus ``flamby`` is the whole budget here.

   ``tests/test_shipped_modules.py`` enforces this by walking the AST. It is a test
   rather than a note because the failure mode -- a ``ModuleNotFoundError`` on a remote
   worker, minutes into a run, after a scheduler queue wait -- is expensive and looks
   like an infrastructure problem rather than an import problem.

   This is not hypothetical. The sibling GWAS experiment originally imported two helper
   modules from its own directory, which is why every partner's endpoint config needed a
   ``PYTHONPATH`` line, and why a missing one was the single largest source of
   partner-side breakage on that project.

2. **The path is resolved on the DRIVER, never on the partner's cluster.** A partner is
   never asked for a ``dataset_path``. An earlier version of the setup guide did ask, and
   the value went unused.

The heavy import is inside the function on purpose: this module is read as text by the
driver, which does not need ``flamby`` installed to ship it.
"""


def get_dataset(dataset: str, num_clients: int, client_id: int, **kwargs):
    """Return ``(train_dataset, test_dataset)`` for one site.

    Train is this site's own center (``pooled=False``); test is the **pooled** test set,
    so every site validates against identical held-out data. That is the benchmark's own
    protocol and it is what makes per-round numbers comparable across sites -- but it
    does mean every site holds the full test set, which is worth stating explicitly in a
    write-up so it is not mistaken for a leak.

    Args:
        dataset: FLamby dataset name. Only ``HeartDisease`` is wired up here.
        num_clients: The dataset's natural split size, which the loader asserts against.
            NOT the number of participating sites -- three sites drawing from a 4-center
            split is fine.
        client_id: Which center this site trains on.
    """
    if dataset != "HeartDisease":
        raise NotImplementedError(
            f"this experiment package handles 'HeartDisease', not {dataset!r}.\n"
            "FLamby ships six other datasets. Adding one means a new experiment package "
            "with its own model, loss and metric -- see docs/experiments/ for the "
            "four documents each experiment needs."
        )

    # Fed-Heart-Disease has four natural centers: the collection sites of the original
    # UCI study (Cleveland, Hungary, Switzerland, VA Long Beach). The heterogeneity is
    # real rather than a synthetic IID slicing, which is the reason for using it.
    if num_clients > 4:
        raise ValueError(
            f"Fed-Heart-Disease has 4 centers, but num_clients={num_clients}.\n"
            "Past four participating sites, either two sites share a center (defensible, "
            "but it must be disclosed) or the experiment moves to a dataset with more "
            "centers -- Fed-TCGA-BRCA has six."
        )
    if not 0 <= client_id < num_clients:
        raise ValueError(
            f"client_id={client_id} is outside [0, {num_clients}). Each site must be "
            "assigned a distinct center; check the `center` values in federation.yaml."
        )

    try:
        from flamby.datasets.fed_heart_disease import FedHeartDisease
    except ImportError as exc:
        raise ImportError(
            "flamby is not importable on this worker.\n"
            "\n"
            "Install it into the SAME environment the endpoint's worker_init activates, "
            "from a GitHub checkout:\n"
            "    git clone https://github.com/owkin/FLamby.git\n"
            "\n"
            "Do NOT run FLamby's `make install` (it builds a separate conda environment "
            "the endpoint never activates) and do NOT run `pip install -e .` without "
            "reading the setup guide first -- FLamby's setup.py hooks egg_info to shell "
            "out to an external package index, which can move pinned packages.\n"
            "\n"
            "See docs/partner/experiments/flamby-heart-disease.md."
        ) from exc

    test_dataset = FedHeartDisease(train=False, pooled=True)
    train_dataset = FedHeartDisease(train=True, center=client_id, pooled=False)
    return train_dataset, test_dataset
