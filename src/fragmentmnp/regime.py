import numpy as np
import matplotlib.pyplot as plt


def classify_release_regime(D_p, K_pw, k_frag, radius_m, D_w=None):
    """
    Simple first-pass qualitative regime classifier.
    """
    tau_diff = (radius_m ** 2) / max(D_p, 1e-30)
    tau_frag = 1.0 / max(k_frag, 1e-30)

    ratio = tau_diff / tau_frag

    if ratio > 10:
        label = "fragmentation-dominated"
    elif ratio < 0.1:
        label = "diffusion-dominated"
    else:
        label = "intermediate"

    return {
        "tau_diff": tau_diff,
        "tau_frag": tau_frag,
        "ratio_diff_to_frag": ratio,
        "label": label,
    }


def add_regime_columns(df, D_p_col, k_frag_col, radius_col, K_pw_col=None, D_w_col=None):
    out = df.copy()
    labels = []
    tau_diff = []
    tau_frag = []
    ratio = []

    for _, row in out.iterrows():
        res = classify_release_regime(
            D_p=row[D_p_col],
            K_pw=row[K_pw_col] if K_pw_col else None,
            k_frag=row[k_frag_col],
            radius_m=row[radius_col],
            D_w=row[D_w_col] if D_w_col else None,
        )
        labels.append(res["label"])
        tau_diff.append(res["tau_diff"])
        tau_frag.append(res["tau_frag"])
        ratio.append(res["ratio_diff_to_frag"])

    out["tau_diff"] = tau_diff
    out["tau_frag"] = tau_frag
    out["ratio_diff_to_frag"] = ratio
    out["regime_label"] = labels
    return out


def plot_regime_map(df, x="tau_frag", y="tau_diff", label_col="regime_label"):
    fig, ax = plt.subplots()
    for label in sorted(df[label_col].unique()):
        sub = df[df[label_col] == label]
        ax.scatter(sub[x], sub[y], label=label)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.legend()
    ax.set_title("Release regime map")
    return fig, ax


def plot_release_metric_map(df, x, y, metric="fraction_released"):
    fig, ax = plt.subplots()
    sc = ax.scatter(df[x], df[y], c=df[metric])
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(x)
    ax.set_ylabel(y)
    ax.set_title(metric)
    plt.colorbar(sc, ax=ax, label=metric)
    return fig, ax