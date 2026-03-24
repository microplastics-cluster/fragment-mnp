import matplotlib.pyplot as plt
import numpy as np


def plot_additive_fate_dashboard(out, additive=0):
    """
    One additive:
    - particulate total
    - released total
    - named medium pools
    """
    additive_name = out.additive_names[additive]
    ts = out.get_additive_timeseries(additive)

    fig, ax = plt.subplots()
    ax.plot(out.t, ts["particulate_total"], label="Particulate total")
    ax.plot(out.t, ts["aqueous"], label="Released total")
    ax.plot(out.t, ts["total"], "--", label="Tracked total")

    if getattr(out, "medium_pool_names", None) is not None:
        for i, pool_name in enumerate(out.medium_pool_names):
            if str(pool_name).startswith(f"{additive_name}:"):
                short = pool_name.split(":", 1)[1]
                ax.plot(out.t, out.c_medium_pool_species[i], label=f"Medium: {short}")

    ax.set_xlabel("Time")
    ax.set_ylabel("Mass concentration")
    ax.set_title(f"Additive fate dashboard: {additive_name}")
    ax.legend()
    return fig, ax


def plot_parent_product_dynamics(out, additive=0, parent_pool="dissolved_parent", product_pool="transformed_product"):
    additive_name = out.additive_names[additive]
    parent_name = f"{additive_name}:{parent_pool}"
    product_name = f"{additive_name}:{product_pool}"

    ip = out.get_medium_pool_index(parent_name)
    it = out.get_medium_pool_index(product_name)

    fig, ax = plt.subplots()
    ax.plot(out.t, out.c_medium_pool_species[ip], label=parent_pool)
    ax.plot(out.t, out.c_medium_pool_species[it], label=product_pool)
    ax.set_xlabel("Time")
    ax.set_ylabel("Mass concentration")
    ax.set_title(f"Parent-product dynamics: {additive_name}")
    ax.legend()
    return fig, ax


def plot_size_class_release_contributions(out, additive=0):
    ts = out.get_additive_timeseries(additive)
    initial = ts["particulate_by_size"][:, 0]
    final = ts["particulate_by_size"][:, -1]
    released = initial - final

    fig, ax = plt.subplots()
    ax.bar(np.arange(len(released)), released)
    ax.set_xlabel("Size class index")
    ax.set_ylabel("Released mass contribution")
    ax.set_title(f"Size-class release contributions: {out.additive_names[additive]}")
    return fig, ax