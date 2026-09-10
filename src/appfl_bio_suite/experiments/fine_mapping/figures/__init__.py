"""Figure generation for the fine-mapping experiment. COORDINATOR-SIDE ONLY.

Never shipped to a worker, so these modules may import freely and may depend on
matplotlib -- which is deliberately absent from a partner's dependency surface.

    figstyle.py     the one visual system: validated palette, rcParams, helpers
    paper_plots.py  the SuSiEx-paper figure analogues, from the results table
    federation.py   what the paper has no analogue for: the cost of federating
    eda.py          the simulated package itself, before any fine-mapping runs
    detail.py       harvest per-credible-set / per-variant records from SuSiEx output

``plotting.py`` beside this package is the aggregator's in-process entry point and
stays where it is; it now delegates to ``paper_plots``.
"""

__all__ = ["figstyle", "paper_plots", "federation", "eda", "detail"]
