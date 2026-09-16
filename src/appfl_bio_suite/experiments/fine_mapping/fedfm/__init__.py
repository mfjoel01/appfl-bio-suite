"""Maintained scientific fork of the imported FedFM simulation modules.

UPSTREAM.json records the original content hashes and explicitly identifies modules
changed after the September 2026 scientific audit. Equality is tested against the
analysis protocol and retained results, not inferred from byte identity with a legacy
source tree. Site code remains duplicated for shipment and tested for exact parity.
"""

__all__ = [
    "utils",
    "sampling",
    "locus_selection",
    "phenotype_sim",
    "qc",
    "validation",
    "fine_mapping",
    "fed_fine_mapping",
]
