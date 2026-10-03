"""Automotive Dual-Clutch Transmission (DCT) Concentric Shafts Fatigue Design.

Mechanical Engineering Machine Design I - Sharif University of Technology
Authors: Zahra Khazaei & Arshia Edalatpajouh | Instructor: Dr. Keyvan Mohammadi

This module implements the Westinghouse (ASME-Elliptic/Goodman) fatigue design
criteria for both solid (inner) and hollow (outer) concentric transmission shafts,
incorporating Marin endurance modification factors, stress concentration factors,
static reaction loads, shear force and bending moment distributions, and keyway sizing.
"""

import math
import os
import matplotlib.pyplot as plt
import numpy as np


class ShaftDesignCalculator:
    """Fatigue analysis and diameter calculation engine for transmission shafts."""

    def __init__(
        self,
        engine_torque: float = 390.0,    # Maximum engine torque (N*m)
        engine_power_hp: float = 254.0,  # Engine rated power (hp)
        pressure_angle_deg: float = 20.0,# Gear tooth pressure angle
        safety_factor: float = 3.0,      # Design factor of safety against fatigue
        reliability: float = 0.9999,     # Target reliability (99.99%)
    ):
        self.t_engine = engine_torque
        self.p_hp = engine_power_hp
        self.phi = math.radians(pressure_angle_deg)
        self.fs = safety_factor
        self.reliability = reliability

        # Material Properties: AISI 1045 Hot Rolled (Shafts)
        self.s_ut_shaft = 570.0  # MPa
        self.s_y_shaft = 310.0   # MPa

        # Material Properties: AISI 1020 Cold Drawn (Keys)
        self.s_ut_key = 1020.0   # MPa
        self.s_y_key = 848.0     # MPa

        # Gear Pitch Radii (m)
        self.r_g1 = 0.150 / 2.0  # D = 150 mm -> r = 75 mm (0.075 m)
        self.r_g3 = 0.250 / 2.0  # D = 250 mm -> r = 125 mm (0.125 m)

    def calculate_marin_endurance_limit(self, d_initial: float = 50.0) -> dict:
        """Calculate the modified endurance limit Se using Marin factors."""
        # Uncorrected rotary beam endurance limit
        s_e_prime = 0.504 * self.s_ut_shaft  # 285.0 MPa

        # Surface condition factor ka (Ground finish: a = 1.58, b = -0.085)
        k_a = 1.58 * (self.s_ut_shaft ** -0.085)  # ~ 0.921

        # Size factor kb (for 51 mm < d <= 254 mm or d ~ 50 mm)
        k_b = 1.24 * (d_initial ** -0.107)        # ~ 0.816 - 0.847

        # Loading factor kc (Bending)
        k_c = 1.0

        # Temperature factor kd (T < 200 °C)
        k_d = 1.0

        # Reliability factor ke (for 99.99% reliability -> z = 3.719)
        k_e = 0.702

        # Miscellaneous factor kf
        k_f = 1.0

        # Modified endurance limit
        s_e = k_a * k_b * k_c * k_d * k_e * k_f * s_e_prime  # ~ 150.39 MPa

        return {
            "s_e_prime": s_e_prime,
            "k_a": k_a, "k_b": k_b, "k_c": k_c,
            "k_d": k_d, "k_e": k_e, "k_f": k_f,
            "s_e": s_e
        }

    @staticmethod
    def westinghouse_solid_diameter(
        fs: float, kf: float, kfs: float,
        mm: float, ma: float, tm: float, ta: float,
        sy: float, se: float
    ) -> float:
        """Compute required solid shaft diameter using the Westinghouse fatigue equation.

        Formula (d in mm, moments and torques in N*m converted with 1000):
        d = [ (32 * 1000 * FS / pi) * sqrt( kf^2 * (Mm/Sy + Ma/Se)^2 + kfs^2 * (Tm/Sy + Ta/Se)^2 ) ]^(1/3)
        """
        bending_term = (mm / sy) + (ma / se)
        term_1 = (kf ** 2) * (bending_term ** 2)
        torsion_term = (tm / sy) + (ta / se)
        term_2 = (kfs ** 2) * (torsion_term ** 2)
        inner_sqrt = math.sqrt(term_1 + term_2)
        coefficient = (32.0 * 1000.0 * fs) / math.pi
        d = (coefficient * inner_sqrt) ** (1.0 / 3.0)
        return d

    @staticmethod
    def westinghouse_hollow_diameter(
        fs: float, kf: float, kfs: float,
        mm: float, ma: float, tm: float, ta: float,
        sy: float, se: float, k_ratio: float
    ) -> float:
        """Compute outer diameter do for hollow shaft given diameter ratio k = di / do.

        do = [ (32 * 1000 * FS / (pi * (1 - k^4))) * sqrt( ... ) ]^(1/3)
        """
        bending_term = (mm / sy) + (ma / se)
        term_1 = (kf ** 2) * (bending_term ** 2)
        torsion_term = (tm / sy) + (ta / se)
        term_2 = (kfs ** 2) * (torsion_term ** 2)
        inner_sqrt = math.sqrt(term_1 + term_2)
        coefficient = (32.0 * 1000.0 * fs) / (math.pi * (1.0 - k_ratio ** 4))
        do = (coefficient * inner_sqrt) ** (1.0 / 3.0)
        return do

    def solve_shaft1_gear1_mode(self) -> dict:
        """Perform static and fatigue evaluation for Shaft 1 under 1st Gear engagement."""
        # Tangential force Wt and radial force Wr on Gear 1 (r = 0.075 m)
        wt = self.t_engine / self.r_g1  # 390 / 0.075 = 5200 N
        wr = wt * math.tan(self.phi)     # 5200 * tan(20°) = 1892.6 N

        # Reaction forces on Bearings B4 (x = 0) and B1 (x = 0.675 m), Gear 1 at x = 0.630 m
        # In vertical xy plane (Wr):
        # sum(M_B4) = 0 => Wr * 0.630 - B1_xy * 0.675 = 0
        b1_xy = (wr * 0.630) / 0.675  # 1766.46 N
        b4_xy = wr - b1_xy            # 126.17 N

        # In horizontal xz plane (Wt):
        b1_xz = (wt * 0.630) / 0.675  # 4853.33 N
        b4_xz = wt - b1_xz            # 346.67 N

        # Critical stepped sections (x in mm):
        # 1: B4 seat (x=10)
        # 2: Spline entry (x=40)
        # 3: Bearing B3 seat (x=170)
        # 4: Bearing B2 seat (x=410)
        # 5: Gear G1 shoulder (x=615)
        # 6: Gear G1 keyway (x=630)
        # 7: Bearing B1 seat (x=665)
        sections = [
            {"id": 1, "name": "Bearing B4 seat", "x": 10, "kf": 1.975, "kfs": 1.595, "has_torque": False},
            {"id": 2, "name": "Clutch Spline",   "x": 40, "kf": 1.975, "kfs": 1.595, "has_torque": True},
            {"id": 3, "name": "Bearing B3 seat", "x": 170, "kf": 1.975, "kfs": 1.595, "has_torque": True},
            {"id": 4, "name": "Bearing B2 seat", "x": 410, "kf": 1.975, "kfs": 1.595, "has_torque": True},
            {"id": 5, "name": "Gear G1 shoulder","x": 615, "kf": 1.975, "kfs": 1.595, "has_torque": True},
            {"id": 6, "name": "Gear G1 keyway",  "x": 630, "kf": 1.300, "kfs": 1.300, "has_torque": False},
            {"id": 7, "name": "Bearing B1 seat", "x": 665, "kf": 1.975, "kfs": 1.595, "has_torque": False},
        ]

        se_info = self.calculate_marin_endurance_limit(50.0)
        se = se_info["s_e"]

        results = []
        for sec in sections:
            x_m = sec["x"] / 1000.0
            # Macaulay bending moments
            m_xy = b4_xy * x_m - (wr * (x_m - 0.630) if x_m > 0.630 else 0.0)
            m_xz = b4_xz * x_m - (wt * (x_m - 0.630) if x_m > 0.630 else 0.0)
            ma = math.sqrt(m_xy**2 + m_xz**2)
            mm = 0.0  # Fully reversed rotating bending

            # Torsion
            if sec["has_torque"]:
                t_max = self.t_engine
                ta = t_max / 2.0  # 195 N*m
                tm = t_max / 2.0  # 195 N*m
            else:
                ta = 0.0
                tm = 0.0

            d_calc = self.westinghouse_solid_diameter(
                self.fs, sec["kf"], sec["kfs"],
                mm, ma, tm, ta, self.s_y_shaft, se
            )

            results.append({
                "id": sec["id"], "name": sec["name"], "x_mm": sec["x"],
                "m_xy": m_xy, "m_xz": m_xz, "ma": ma, "ta": ta, "tm": tm,
                "d_calc": d_calc
            })

        return {
            "wt": wt, "wr": wr,
            "b1_xy": b1_xy, "b4_xy": b4_xy,
            "b1_xz": b1_xz, "b4_xz": b4_xz,
            "sections": results
        }

    def solve_shaft2_hollow_iterative(self, di: float = 86.0) -> list:
        """Solve hollow shaft (Shaft 2) outer diameter using iterative trial-and-error for k."""
        # Representative section under critical bending Ma = 204.39 N*m, Ta = Tm = 195 N*m
        ma = 204.39
        mm = 0.0
        ta = 195.0
        tm = 195.0
        kf = 1.3
        kfs = 1.3
        se = 138.22  # MPa (using kb ~ 0.75 for hollow pipe)

        k_trials = [0.90, 0.92, 0.94, 0.96, 0.97, 0.974, 0.9743, 0.98]
        convergence_log = []

        for k in k_trials:
            do_calc = self.westinghouse_hollow_diameter(
                self.fs, kf, kfs, mm, ma, tm, ta, self.s_y_shaft, se, k
            )
            k_real = di / do_calc
            diff = abs(k - k_real)
            convergence_log.append({
                "k_assumed": k,
                "do_calculated": do_calc,
                "k_real": k_real,
                "error": diff
            })

        return convergence_log


def generate_shaft_plots(output_dir: str):
    """Generate visual engineering plots for bending moments and shaft diameters."""
    os.makedirs(output_dir, exist_ok=True)
    calc = ShaftDesignCalculator()
    s1_res = calc.solve_shaft1_gear1_mode()

    # 1. Bending Moment Diagram along Shaft 1
    x_span = np.linspace(0, 0.675, 500)
    b4_xy = s1_res["b4_xy"]
    b4_xz = s1_res["b4_xz"]
    wr = s1_res["wr"]
    wt = s1_res["wt"]

    m_xy = np.where(x_span <= 0.630, b4_xy * x_span, b4_xy * x_span - wr * (x_span - 0.630))
    m_xz = np.where(x_span <= 0.630, b4_xz * x_span, b4_xz * x_span - wt * (x_span - 0.630))
    m_total = np.sqrt(m_xy**2 + m_xz**2)

    plt.figure(figsize=(8, 4.5), dpi=300)
    plt.plot(x_span * 1000, m_xy, "b--", label="Vertical Plane $M_{xy}$")
    plt.plot(x_span * 1000, m_xz, "g-.", label="Horizontal Plane $M_{xz}$")
    plt.plot(x_span * 1000, m_total, "r-", lw=2, label=r"Total Resultant $M_a = \sqrt{M_{xy}^2 + M_{xz}^2}$")

    plt.xlabel("Axial Position along Shaft 1 (mm)", fontsize=11, fontweight="bold")
    plt.ylabel("Bending Moment (N·m)", fontsize=11, fontweight="bold")
    plt.title("Shaft 1 Bending Moment Distributions (1st Gear Mode)", fontsize=12, fontweight="bold")
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.axvline(630, color="gray", linestyle=":", label="Gear 1 Location (x = 630 mm)")
    plt.legend(frameon=True)
    plt.tight_layout()

    out_p1 = os.path.join(output_dir, "shaft1_bending_moments_plot.png")
    plt.savefig(out_p1)
    plt.close()
    print(f"Saved: {out_p1}")

    # 2. Calculated vs Standardized Diameters Comparison
    sec_names = [s["name"] for s in s1_res["sections"]]
    d_calculated = [s["d_calc"] for s in s1_res["sections"]]
    d_standard = [45, 48, 55, 55, 52, 52, 45]  # Standardized manufactured sizes

    x_indices = np.arange(len(sec_names))
    width = 0.35

    plt.figure(figsize=(9, 4.5), dpi=300)
    plt.bar(x_indices - width/2, d_calculated, width, label="Calculated Westinghouse Min. Dia", color="royalblue")
    plt.bar(x_indices + width/2, d_standard, width, label="Selected Standard Diameter", color="forestgreen")

    plt.ylabel("Diameter (mm)", fontsize=11, fontweight="bold")
    plt.title("Shaft 1 Diameter Sizing: Calculated Minimum vs Selected Standard", fontsize=12, fontweight="bold")
    plt.xticks(x_indices, sec_names, rotation=25, ha="right", fontsize=9, fontweight="bold")
    plt.grid(True, axis="y", linestyle="--", alpha=0.6)
    plt.legend(frameon=True)
    plt.tight_layout()

    out_p2 = os.path.join(output_dir, "shaft1_diameter_sizing_comparison.png")
    plt.savefig(out_p2)
    plt.close()
    print(f"Saved: {out_p2}")


if __name__ == "__main__":
    calc = ShaftDesignCalculator()
    res1 = calc.solve_shaft1_gear1_mode()
    print("=" * 68)
    print("SHAFT 1 FATIGUE SIZING RESULTS (1ST GEAR ENGAGEMENT)")
    print("=" * 68)
    for s in res1["sections"]:
        print(f"Section {s['id']} ({s['name']:<18}) @ x = {s['x_mm']:3d} mm : Ma = {s['ma']:6.2f} N*m | d_min = {s['d_calc']:5.2f} mm")
    print("=" * 68)

    print("\nSHAFT 2 (HOLLOW) ITERATIVE CONVERGENCE SAMPLE:")
    iter_log = calc.solve_shaft2_hollow_iterative(86.0)
    for row in iter_log:
        print(f"k_assumed: {row['k_assumed']:.4f} -> do_calc: {row['do_calculated']:6.2f} mm -> k_real: {row['k_real']:.4f} | Error: {row['error']:.4f}")

    figures_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "docs", "images")
    generate_shaft_plots(figures_dir)
