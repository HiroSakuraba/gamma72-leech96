# Global exact Leech96 minimum-norm investigation

**No extremal 96-dimensional lattice is proved here.** This bundle adds a globally exhaustive *method*, with exact radius-pruning and sharded checkpoint/restart, to the existing Z/16 lattice candidate. Only one of 21 sign-normalized top-coordinate partitions is fully excluded by the new search.

Start with [GLOBAL_EXACT_SVP_2026-10-13.md](GLOBAL_EXACT_SVP_2026-10-13.md). Run `python test_global.py`. The input matrices in `z16_route/` come from the preceding construction experiment and are all rechecked for exact basis equivalence.

Everything in `global_state_110k.json` is explicitly **INCOMPLETE**. No GitHub merge or mathematical minimum claim should be made on its basis.
