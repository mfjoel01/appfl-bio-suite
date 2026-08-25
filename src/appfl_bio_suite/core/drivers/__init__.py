"""Run drivers: the process that actually executes a federation.

Two, deliberately:

``globus_compute_driver``
    Production. Dispatches to remote endpoints across institutions.

``serial_driver``
    Loopback. Coordinator and all sites in one process, no Globus Compute and no
    scheduler. This is what lets someone validate an install before recruiting a single
    partner, and it is the only federated path CI can exercise. Built on APPFL's own
    non-Globus runner rather than a mock transport, so what it exercises is the real
    aggregation and training path.
"""
