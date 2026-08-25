"""FLamby Fed-Heart-Disease: multi-round FedAvg over a public tabular benchmark.

Partners obtain the data themselves -- it is public, ~40 KB, and has no registration
gate -- so there is no simulation stage and nothing for the coordinator to distribute.

`dataset.py` and `model.py` are shipped to workers as source text and must stay
self-contained. See dataset.py's module docstring for what that means.
"""
