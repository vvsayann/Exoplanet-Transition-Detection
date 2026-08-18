

import numpy as np
import matplotlib.pyplot as plt
import lightkurve as lk

TARGET_NAME = "Kepler-10"
MISSION = "Kepler"
FLATTEN_WINDOW = 401
BLS_PERIOD_RANGE = np.arange(0.5, 2.0, 0.0005)
PUBLISHED_PERIOD_DAYS = 0.837

FIGURE_DIR = "figures"


def download_light_curve(target_name, mission):

    search_result = lk.search_lightcurve(target_name, mission=mission)
    print(search_result)
    lc = search_result[0].download()
    return lc


def clean_and_flatten(lc, window_length=FLATTEN_WINDOW):
    lc_clean = lc.remove_nans().remove_outliers(sigma=5)
    lc_flat = lc_clean.flatten(window_length=window_length)

    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    lc_clean.plot(ax=axes[0])
    axes[0].set_title("Before flattening")
    lc_flat.plot(ax=axes[1])
    axes[1].set_title("After flattening")
    plt.tight_layout()
    plt.savefig(f"{FIGURE_DIR}/before_after_flatten.png", dpi=150)
    plt.show()

    return lc_flat


def run_bls_search(lc_flat, period_range=BLS_PERIOD_RANGE):
    periodogram = lc_flat.to_periodogram(
        method="bls", period=period_range, frequency_factor=500
    )
    periodogram.plot()
    plt.savefig(f"{FIGURE_DIR}/bls_periodogram.png", dpi=150)
    plt.show()

    best_period = periodogram.period_at_max_power
    best_t0 = periodogram.transit_time_at_max_power
    best_duration = periodogram.duration_at_max_power

    print(f"Best period:   {best_period}")
    print(f"Best t0:       {best_t0}")
    print(f"Best duration: {best_duration}")

    return best_period, best_t0, best_duration


def fold_and_verify(lc_flat, best_period, best_t0, published_period=None):
    folded = lc_flat.fold(period=best_period, epoch_time=best_t0)

    folded.plot()
    plt.savefig(f"{FIGURE_DIR}/folded_transit.png", dpi=150)
    plt.show()

    folded.scatter()
    plt.xlim(-0.1, 0.1)
    plt.title("Zoomed near transit (phase 0)")
    plt.show()

    if published_period is not None:
        diff = abs(best_period.value - published_period)
        pct_diff = 100 * diff / published_period
        print(f"Published period: {published_period} d")
        print(f"Recovered period: {best_period.value:.6f} d")
        print(f"Difference: {pct_diff:.4f}%")

    return folded


if __name__ == "__main__":
    import os
    os.makedirs(FIGURE_DIR, exist_ok=True)

    lc = download_light_curve(TARGET_NAME, MISSION)
    lc_flat = clean_and_flatten(lc)
    best_period, best_t0, best_duration = run_bls_search(lc_flat)
    folded = fold_and_verify(lc_flat, best_period, best_t0, PUBLISHED_PERIOD_DAYS)