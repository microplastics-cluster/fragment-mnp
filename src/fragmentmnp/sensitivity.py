import copy
import itertools
from typing import Any, Dict, Iterable, List

import numpy as np
import pandas as pd

from .fragmentmnp import FragmentMNP


def _set_nested(d: dict, path: str, value):
    keys = path.split(".")
    cur = d
    for key in keys[:-1]:
        if key.isdigit():
            key = int(key)
        cur = cur[key]
    last = keys[-1]
    if last.isdigit():
        last = int(last)
    cur[last] = value


def _get_additive_summary(out, additive=0):
    ts = out.get_additive_timeseries(additive)
    part = ts["particulate_total"]
    aq = ts["aqueous"]
    total = ts["total"]

    initial_total = float(total[0])
    final_total = float(total[-1])
    frac_released = float(aq[-1] / initial_total) if initial_total > 0 else np.nan

    def t_at_fraction(series, frac):
        target = frac * initial_total
        idx = np.where(series >= target)[0]
        if len(idx) == 0:
            return np.nan
        i = idx[0]
        return float(out.t[i])

    return {
        "initial_total": initial_total,
        "final_total": final_total,
        "final_particulate": float(part[-1]),
        "final_aqueous": float(aq[-1]),
        "fraction_released": frac_released,
        "t10_release": t_at_fraction(aq, 0.10),
        "t50_release": t_at_fraction(aq, 0.50),
        "t90_release": t_at_fraction(aq, 0.90),
    }


def run_sensitivity_grid(
    base_config: Dict[str, Any],
    base_data: Dict[str, Any],
    parameter_grid: Dict[str, Iterable[Any]],
    additive=0,
) -> pd.DataFrame:
    """
    parameter_grid keys are dotted paths into config/data:
      config.dt
      data.k_frag
      data.additives.0.pools.0.release.params.D_p
      data.additives.0.pools.0.initial_concs
    """
    keys = list(parameter_grid.keys())
    values_product = list(itertools.product(*(parameter_grid[k] for k in keys)))

    records: List[Dict[str, Any]] = []

    for run_id, values in enumerate(values_product):
        cfg = copy.deepcopy(base_config)
        dat = copy.deepcopy(base_data)

        for key, value in zip(keys, values):
            if key.startswith("config."):
                _set_nested(cfg, key[len("config."):], value)
            elif key.startswith("data."):
                _set_nested(dat, key[len("data."):], value)
            else:
                raise ValueError(f"Grid key must start with 'config.' or 'data.': {key}")

        out = FragmentMNP(cfg, dat).run()
        summary = _get_additive_summary(out, additive=additive)

        record = {"run_id": run_id}
        for key, value in zip(keys, values):
            record[key] = value
        record.update(summary)
        records.append(record)

    return pd.DataFrame.from_records(records)