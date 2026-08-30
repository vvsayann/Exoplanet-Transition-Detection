#made on 26-08-2026 by ayann

import os
import numpy as np
import matplotlib.pyplot as plt
import lightkurve as lk

MISSION = "Kepler"
FLATTEN_WINDOW = 401

FIGURE_DIR = "figures"


TARGETS = [
    {"name": "Kepler-10", "period_range": np.arange(0.5, 2.0, 0.0005),
     "published_period": 0.837, "star_radius_solar": 1.065},
    {"name": "Kepler-8", "period_range": np.arange(3.4, 3.6, 0.0002),
     "published_period": 3.523, "star_radius_solar": 1.486},
    {"name": "Kepler-7", "period_range": np.arange(4.7, 5.0, 0.0002),
     "published_period": 4.886, "star_radius_solar": 1.966},
]


def download_light_curve(target_name, mission):
    search_result = lk.search_lightcurve(target_name, mission=mission, exptime=1800)
    print(search_result)
    lc = search_result[0].download()
    return lc


def clean_and_flatten(lc, window_length=FLATTEN_WINDOW, save_prefix=""):
    lc_clean = lc.remove_nans().remove_outliers(sigma=5)
    lc_flat = lc_clean.flatten(window_length=window_length)

    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
    lc_clean.plot(ax=axes[0])
    axes[0].set_title("Before flattening")
    lc_flat.plot(ax=axes[1])
    axes[1].set_title("After flattening")
    plt.tight_layout()
    plt.savefig(f"{FIGURE_DIR}/{save_prefix}before_after_flatten.png", dpi=150)
    plt.show()

    return lc_flat


def run_bls_search(lc_flat, period_range, save_prefix=""):
    periodogram = lc_flat.to_periodogram(
        method="bls", period=period_range, frequency_factor=500
    )
    periodogram.plot()
    plt.savefig(f"{FIGURE_DIR}/{save_prefix}bls_periodogram.png", dpi=150)
    plt.show()

    best_period = periodogram.period_at_max_power
    best_t0 = periodogram.transit_time_at_max_power
    best_duration = periodogram.duration_at_max_power

    print(f"Best period:   {best_period}")
    print(f"Best t0:       {best_t0}")
    print(f"Best duration: {best_duration}")

    return best_period, best_t0, best_duration


def fold_and_verify(lc_flat, best_period, best_t0, published_period=None, save_prefix=""):
    folded = lc_flat.fold(period=best_period, epoch_time=best_t0)

    folded.plot()
    plt.savefig(f"{FIGURE_DIR}/{save_prefix}folded_transit.png", dpi=150)
    plt.show()

    folded.scatter()
    plt.xlim(-0.1, 0.1)
    plt.title("Zoomed near transit (phase 0)")
    plt.savefig(f"{FIGURE_DIR}/{save_prefix}folded_transit_zoom.png", dpi=150)
    plt.show()

    if published_period is not None:
        diff = abs(best_period.value - published_period)
        pct_diff = 100 * diff / published_period
        print(f"Published period: {published_period} d")
        print(f"Recovered period: {best_period.value:.6f} d")
        print(f"Difference: {pct_diff:.4f}%")

    return folded



def estimate_planet_radius(folded, star_radius_solar):
    """Estimate planet radius from transit depth: depth = (R_planet/R_star)^2."""
    in_transit = folded.flux[np.abs(folded.phase.value) < 0.02]
    depth = 1.0 - np.nanmin(in_transit.value)

    planet_radius_solar = star_radius_solar * np.sqrt(depth)
    planet_radius_earth = planet_radius_solar * 109.2

    print(f"Transit depth: {depth:.6f} ({depth * 1e6:.1f} ppm)")
    print(f"Estimated planet radius: {planet_radius_earth:.2f} Earth radii")

    return planet_radius_earth


def run_full_pipeline(target_name, mission, period_range,
                       published_period=None, star_radius_solar=None,
                       flatten_window=FLATTEN_WINDOW):

    print(f"\n{'=' * 50}\nRunning pipeline for: {target_name}\n{'=' * 50}")
    save_prefix = f"{target_name.replace(' ', '_')}_"

    lc = download_light_curve(target_name, mission)
    lc_flat = clean_and_flatten(lc, window_length=flatten_window, save_prefix=save_prefix)
    best_period, best_t0, best_duration = run_bls_search(lc_flat, period_range, save_prefix=save_prefix)
    folded = fold_and_verify(lc_flat, best_period, best_t0, published_period, save_prefix=save_prefix)

    radius = None
    if star_radius_solar is not None:
        radius = estimate_planet_radius(folded, star_radius_solar)

    return {
        "target": target_name,
        "period": best_period.value,
        "duration": best_duration.value,
        "radius_earth": radius,
    }


if __name__ == "__main__":
    os.makedirs(FIGURE_DIR, exist_ok=True)

    results = []
    for t in TARGETS:
        result = run_full_pipeline(
            t["name"],
            MISSION,
            t["period_range"],
            published_period=t["published_period"],
            star_radius_solar=t["star_radius_solar"],
        )
        results.append(result)

    print(f"\n{'=' * 50}\nSummary\n{'=' * 50}")
    for r in results:
        print(r)