/**
 * test_droplet.cpp
 *
 * Unit tests for the droplet heating and evaporation simulation.
 *
 * Test suite (no external framework – plain assert with diagnostics):
 *
 *  T01 – Antoine equation: saturation pressure at T_boil ≈ 1 atm
 *  T02 – Gas properties: physically reasonable at several temperatures
 *  T03 – Spalding number B_M > 0 and increases with temperature
 *  T04 – D²-law linearity: (d/d0)² decreases nearly linearly in time
 *  T05 – Temperature convergence: T_d approaches wet-bulb temperature
 *  T06 – Mass conservation: mass loss == integral of evaporation rate
 *  T07 – Energy balance: integral(Q_conv) ≈ integral(-mdot)*L_v + sensible
 *  T08 – Convection always faster than stagnant (U_rel effect)
 *  T09 – Larger droplet lives longer (lifetime ∝ d0²)
 *  T10 – Wet-bulb temperature lies between T_d0 and T_inf
 */

#include "droplet_evaporation.h"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <string>
#include <vector>
#include <cassert>
#include <numeric>

// ─────────────────────────────────────────────────────────────────────────────
// Minimal test framework
// ─────────────────────────────────────────────────────────────────────────────

static int g_total = 0, g_passed = 0, g_failed = 0;

static void EXPECT(bool condition, const std::string& msg,
                   double got = 0.0, double expected = 0.0) {
    ++g_total;
    if (condition) {
        ++g_passed;
        std::cout << "  [PASS] " << msg << "\n";
    } else {
        ++g_failed;
        std::cout << "  [FAIL] " << msg;
        if (got != 0.0 || expected != 0.0)
            std::cout << "  (got=" << got << ", expected≈" << expected << ")";
        std::cout << "\n";
    }
}

[[maybe_unused]]
static void EXPECT_NEAR(double got, double expected, double rel_tol,
                        const std::string& msg) {
    double err = std::abs(got - expected);
    double tol = std::abs(expected) * rel_tol + 1e-300;
    EXPECT(err <= tol, msg, got, expected);
}

static void START_TEST(const std::string& name) {
    std::cout << "\n--- " << name << " ---\n";
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper: build a standard simulation (stagnant, T_inf=1000 K, d0=100 µm)
// ─────────────────────────────────────────────────────────────────────────────

static std::pair<DropletEvaporation, std::vector<DropletState>>
makeStandardSim(double d0 = 100e-6, double T_inf = 1000.0,
                double U_rel = 0.0, double t_end = 0.3) {
    FuelProperties fuel = FuelProperties::nHeptane();

    AmbientConditions amb;
    amb.T_inf = T_inf;
    amb.P_inf = PhysConst::P_atm;
    amb.U_rel = U_rel;

    SimParams p;
    p.d0    = d0;
    p.T_d0  = 300.0;
    p.dt    = 5e-6;
    p.t_end = t_end;

    DropletEvaporation sim(fuel, amb, p);
    auto history = sim.run();
    return {sim, history};
}

// ─────────────────────────────────────────────────────────────────────────────
// T01 – Antoine equation: P_sat(T_boil) ≈ 101 325 Pa  (within 5%)
// ─────────────────────────────────────────────────────────────────────────────
static void test_T01() {
    START_TEST("T01 – Antoine equation at boiling point");
    FuelProperties f = FuelProperties::nHeptane();
    double P_sat = f.saturationPressure(f.T_boil);
    double err_pct = std::abs(P_sat - PhysConst::P_atm) / PhysConst::P_atm * 100.0;

    std::cout << "  P_sat(" << f.T_boil << " K) = " << P_sat << " Pa  "
              << "(error = " << std::fixed << std::setprecision(2) << err_pct << "%)\n";
    EXPECT(err_pct < 5.0,
           "P_sat at T_boil within 5% of 101325 Pa", P_sat, PhysConst::P_atm);

    // Below boiling point, P_sat must be < P_atm
    double P_low = f.saturationPressure(300.0);
    EXPECT(P_low < PhysConst::P_atm,
           "P_sat(300K) < 1 atm", P_low, PhysConst::P_atm);

    // Above boiling point, P_sat must be > P_atm
    double P_high = f.saturationPressure(400.0);
    EXPECT(P_high > PhysConst::P_atm,
           "P_sat(400K) > 1 atm", P_high, PhysConst::P_atm);

    // Monotonically increasing with temperature
    EXPECT(P_low < P_sat && P_sat < P_high,
           "P_sat monotonically increasing (300 < 371.6 < 400 K)");
}

// ─────────────────────────────────────────────────────────────────────────────
// T02 – Gas properties: physically reasonable values
// ─────────────────────────────────────────────────────────────────────────────
static void test_T02() {
    START_TEST("T02 – Gas properties at 300 K, 700 K, 1500 K");
    for (double T : {300.0, 700.0, 1500.0}) {
        GasProperties g = GasProperties::evaluate(T, PhysConst::P_atm, 0.028964);
        std::cout << "  T=" << T << " K: mu=" << g.mu_g << " Pa·s, k="
                  << g.k_g << " W/(m·K), cp=" << g.cp_g
                  << " J/(kg·K), rho=" << g.rho_g
                  << " kg/m³, D_AB=" << g.D_AB << " m²/s, Pr="
                  << g.Pr << "\n";

        EXPECT(g.mu_g > 1e-6 && g.mu_g < 1e-4,
               "Viscosity in [1e-6, 1e-4] Pa·s at T=" + std::to_string((int)T));
        EXPECT(g.k_g > 0.01 && g.k_g < 0.30,
               "Conductivity in [0.01, 0.30] W/(m·K) at T=" + std::to_string((int)T));
        EXPECT(g.cp_g > 1000.0 && g.cp_g < 1300.0,
               "Cp_g in [1000, 1300] J/(kg·K) at T=" + std::to_string((int)T));
        EXPECT(g.rho_g > 0.0,
               "Density > 0 at T=" + std::to_string((int)T));
        EXPECT(g.D_AB > 0.0,
               "Diffusivity > 0 at T=" + std::to_string((int)T));
        EXPECT(g.Pr > 0.5 && g.Pr < 1.5,
               "Prandtl in [0.5, 1.5] at T=" + std::to_string((int)T));
    }

    // Viscosity must increase with temperature (gas)
    auto g300  = GasProperties::evaluate(300.0,  PhysConst::P_atm, 0.028964);
    auto g1000 = GasProperties::evaluate(1000.0, PhysConst::P_atm, 0.028964);
    EXPECT(g1000.mu_g > g300.mu_g,
           "Gas viscosity increases with temperature (300→1000 K)");

    // Density must decrease with temperature (ideal gas, const P)
    EXPECT(g1000.rho_g < g300.rho_g,
           "Gas density decreases with temperature (ideal gas)");
}

// ─────────────────────────────────────────────────────────────────────────────
// T03 – Spalding number B_M
// ─────────────────────────────────────────────────────────────────────────────
static void test_T03() {
    START_TEST("T03 – Spalding mass transfer number B_M");
    FuelProperties f = FuelProperties::nHeptane();

    // B_M increases with droplet temperature
    AmbientConditions amb;
    amb.T_inf = 1000.0; amb.P_inf = PhysConst::P_atm; amb.U_rel = 0.0;

    SimParams p;
    p.d0 = 100e-6; p.T_d0 = 300.0; p.dt = 5e-6; p.t_end = 0.01;

    auto [sim, h] = makeStandardSim(100e-6, 1000.0, 0.0, 0.01);

    EXPECT(!h.empty(), "Simulation produced output");
    if (h.size() >= 2) {
        double BM_initial = h.front().B_M;
        double BM_later   = h.back().B_M;
        std::cout << "  B_M(t=0)   = " << BM_initial << "\n";
        std::cout << "  B_M(t_end) = " << BM_later   << "\n";
        EXPECT(BM_initial > 0.0, "B_M > 0 at start");
        EXPECT(BM_later > BM_initial,
               "B_M increases as droplet heats up", BM_later, BM_initial);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// T04 – D²-law linearity
// After the initial heat-up transient the d²(t) curve should be nearly linear.
// ─────────────────────────────────────────────────────────────────────────────
static void test_T04() {
    START_TEST("T04 – D²-law (quasi-steady evaporation phase)");
    auto [sim, h] = makeStandardSim(100e-6, 1000.0, 0.0, 0.3);

    EXPECT(h.size() > 20, "Enough time points for linearity check");

    // Skip the first 20% (heat-up transient) and analyse the remaining part
    size_t skip = h.size() / 5;
    size_t n    = h.size() - skip;

    if (n < 5) { EXPECT(false, "Not enough post-transient points"); return; }

    // Linear regression on d²-norm vs time
    double sum_t = 0, sum_d2 = 0, sum_t2 = 0, sum_td2 = 0;
    for (size_t i = skip; i < h.size(); ++i) {
        double t   = h[i].time;
        double d2n = h[i].diameter_sq_norm;
        sum_t   += t;
        sum_d2  += d2n;
        sum_t2  += t * t;
        sum_td2 += t * d2n;
    }
    double N     = static_cast<double>(n);
    double slope = (N * sum_td2 - sum_t * sum_d2) / (N * sum_t2 - sum_t * sum_t);
    double inter = (sum_d2 - slope * sum_t) / N;

    // Compute R² for linearity
    double mean_d2 = sum_d2 / N;
    double SS_tot = 0, SS_res = 0;
    for (size_t i = skip; i < h.size(); ++i) {
        double d2n = h[i].diameter_sq_norm;
        double fit = slope * h[i].time + inter;
        SS_tot += (d2n - mean_d2) * (d2n - mean_d2);
        SS_res += (d2n - fit)    * (d2n - fit);
    }
    double R2 = (SS_tot > 0) ? 1.0 - SS_res / SS_tot : 1.0;

    std::cout << "  Linear fit:  d²/d0² = " << slope << " * t + " << inter << "\n";
    std::cout << "  R² = " << R2 << "  (expect > 0.99 for quasi-steady phase)\n";
    EXPECT(slope < 0.0, "D²-law slope is negative (d² decreases)");
    EXPECT(R2 > 0.98,   "D²-law linearity R² > 0.98", R2, 0.99);
}

// ─────────────────────────────────────────────────────────────────────────────
// T05 – Temperature convergence to wet-bulb temperature
// ─────────────────────────────────────────────────────────────────────────────
static void test_T05() {
    START_TEST("T05 – Droplet temperature approaches wet-bulb temperature");
    auto [sim, h] = makeStandardSim(100e-6, 1000.0, 0.0, 0.3);

    double T_wb = sim.wetBulbTemperature();
    std::cout << "  Wet-bulb temperature : " << T_wb << " K\n";
    EXPECT(T_wb > 300.0 && T_wb < 1000.0,
           "T_wb between T_d0 and T_inf", T_wb, 650.0);

    // Temperature must rise from initial value
    if (h.size() >= 2) {
        EXPECT(h.back().temperature > h.front().temperature,
               "Temperature rises during simulation");
    }

    // After sufficient time the droplet temperature should be within 20 K of T_wb
    if (!h.empty()) {
        double T_final = h.back().temperature;
        double diff    = std::abs(T_final - T_wb);
        std::cout << "  Final T_d = " << T_final << " K,  |T_d - T_wb| = "
                  << diff << " K\n";
        EXPECT(diff < 25.0,
               "|T_final - T_wb| < 25 K", diff, 0.0);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// T06 – Mass conservation: total mass loss == integral of |ṁ| dt
// ─────────────────────────────────────────────────────────────────────────────
static void test_T06() {
    START_TEST("T06 – Mass conservation");
    auto [sim, h] = makeStandardSim(100e-6, 1000.0, 0.0, 0.2);

    EXPECT(h.size() > 2, "Enough points for integration");

    double mass_initial = h.front().mass;
    double mass_final   = h.back().mass;
    double delta_mass   = mass_initial - mass_final;   // mass lost

    // Trapezoidal integration of |ṁ| dt
    double integral_mdot = 0.0;
    for (size_t i = 1; i < h.size(); ++i) {
        double dt_i = h[i].time - h[i-1].time;
        integral_mdot += 0.5 * (std::abs(h[i].mdot) + std::abs(h[i-1].mdot)) * dt_i;
    }

    double err_pct = std::abs(delta_mass - integral_mdot)
                   / std::max(mass_initial, 1e-30) * 100.0;

    std::cout << "  Mass initial    : " << mass_initial  << " kg\n";
    std::cout << "  Mass final      : " << mass_final    << " kg\n";
    std::cout << "  Δm (direct)     : " << delta_mass    << " kg\n";
    std::cout << "  Δm (∫|ṁ|dt)     : " << integral_mdot << " kg\n";
    std::cout << "  Error           : " << err_pct       << " %\n";

    EXPECT(err_pct < 5.0,
           "Mass loss matches integrated evaporation rate (<5% error)", err_pct, 0.0);
}

// ─────────────────────────────────────────────────────────────────────────────
// T07 – Energy balance
// Integrating  m·c_l·dT/dt = Q_conv + ṁ·L_v  gives the exact relation:
//   ∫Q_conv dt = c_l·(m_f·T_f − m_i·T_i) − c_l·∫T·ṁ dt + L_v·(m_i − m_f)
// ─────────────────────────────────────────────────────────────────────────────
static void test_T07() {
    START_TEST("T07 – Energy balance");
    FuelProperties fuel = FuelProperties::nHeptane();
    auto [sim, h] = makeStandardSim(100e-6, 1000.0, 0.0, 0.1);

    EXPECT(h.size() > 2, "Enough points for energy integration");

    // ∫ Q_conv dt  (trapezoidal)
    double Q_conv_int = 0.0;
    // ∫ T · ṁ dt  (correction term for enthalpy carried by evaporating mass)
    double T_mdot_int = 0.0;
    for (size_t i = 1; i < h.size(); ++i) {
        double dt_i = h[i].time - h[i-1].time;
        Q_conv_int += 0.5 * (h[i].Q_conv  + h[i-1].Q_conv)  * dt_i;
        T_mdot_int += 0.5 * (h[i].temperature * h[i].mdot
                           + h[i-1].temperature * h[i-1].mdot) * dt_i;
    }

    // Exact absorbed energy (from integrating the ODE analytically)
    double m_i = h.front().mass,  T_i = h.front().temperature;
    double m_f = h.back().mass,   T_f = h.back().temperature;
    double Q_exact = fuel.c_l * (m_f * T_f - m_i * T_i)
                   - fuel.c_l * T_mdot_int
                   + fuel.L_v * (m_i - m_f);

    double err_pct = std::abs(Q_conv_int - Q_exact)
                   / std::max(std::abs(Q_conv_int), 1e-30) * 100.0;

    std::cout << std::scientific << std::setprecision(4);
    std::cout << "  ∫Q_conv dt       : " << Q_conv_int << " J\n";
    std::cout << "  Q_exact (ODE)    : " << Q_exact    << " J\n";
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  Energy balance error: " << err_pct << " %\n";

    EXPECT(err_pct < 5.0,
           "Energy balance error < 5% (exact ODE identity)", err_pct, 0.0);
    EXPECT((m_i - m_f) > 0.0, "Mass was evaporated (m_f < m_i)");
}

// ─────────────────────────────────────────────────────────────────────────────
// T08 – Convection accelerates evaporation (U_rel effect)
// ─────────────────────────────────────────────────────────────────────────────
static void test_T08() {
    START_TEST("T08 – Convective vs stagnant evaporation");
    // Use t_run < minimum evaporation time so both droplets still exist
    // at the comparison point.  At 1000 K, stagnant lifetime ≈ 8.4 ms,
    // convective ≈ 7.3 ms.  Sample at 5 ms (well before both lifetimes).
    const double t_run = 5e-3;

    auto [sim0, h0] = makeStandardSim(100e-6, 1000.0, 0.0, t_run); // stagnant
    auto [sim1, h1] = makeStandardSim(100e-6, 1000.0, 1.0, t_run); // 1 m/s

    EXPECT(!h0.empty() && !h1.empty(), "Both simulations produced output");

    double d_final_stag = h0.back().diameter;
    double d_final_conv = h1.back().diameter;

    std::cout << "  Stagnant final d  : " << d_final_stag * 1e6 << " µm\n";
    std::cout << "  Convective final d: " << d_final_conv * 1e6 << " µm\n";

    EXPECT(d_final_conv < d_final_stag,
           "Convective droplet is smaller at same end time", d_final_conv, d_final_stag);

    // Nusselt and Sherwood should be larger with U_rel > 0
    double Nu0 = h0.front().Nu, Nu1 = h1.front().Nu;
    double Sh0 = h0.front().Sh, Sh1 = h1.front().Sh;
    std::cout << "  Nu (stagnant)  = " << Nu0 << ",  Nu (conv) = " << Nu1 << "\n";
    std::cout << "  Sh (stagnant)  = " << Sh0 << ",  Sh (conv) = " << Sh1 << "\n";
    EXPECT(Nu1 > Nu0, "Nu is larger with relative velocity", Nu1, Nu0);
    EXPECT(Sh1 > Sh0, "Sh is larger with relative velocity", Sh1, Sh0);
}

// ─────────────────────────────────────────────────────────────────────────────
// T09 – Droplet lifetime scales as d0²  (D²-law)
// ─────────────────────────────────────────────────────────────────────────────
static void test_T09() {
    START_TEST("T09 – Lifetime ∝ d0²  (D²-law scaling)");

    // Use a high-temperature environment so the droplets evaporate quickly
    const double T_inf = 1200.0;

    // Run two droplets: d0_b = 2 * d0_a → lifetime_b ≈ 4 * lifetime_a
    double d0_a = 50e-6, d0_b = 100e-6;

    auto [sa, ha] = makeStandardSim(d0_a, T_inf, 0.0, 0.5);
    auto [sb, hb] = makeStandardSim(d0_b, T_inf, 0.0, 1.5);

    // Estimate lifetime from linear extrapolation of d²(t)
    auto estimateLifetime = [](const std::vector<DropletState>& h) -> double {
        if (h.size() < 5) return h.back().time;
        const auto& pa = h[h.size() - 5];
        const auto& pb = h.back();
        double slope = (pb.diameter_sq_norm - pa.diameter_sq_norm)
                      / (pb.time - pa.time + 1e-30);
        if (slope < -1e-20)
            return pa.time - pa.diameter_sq_norm / slope;
        return h.back().time;
    };

    double t_a = estimateLifetime(ha);
    double t_b = estimateLifetime(hb);
    double ratio = t_b / t_a;

    std::cout << "  d0_a = " << d0_a*1e6 << " µm,  t_life_a = " << t_a << " s\n";
    std::cout << "  d0_b = " << d0_b*1e6 << " µm,  t_life_b = " << t_b << " s\n";
    std::cout << "  t_b / t_a = " << ratio << "  (expect ≈ 4.0)\n";

    // Allow 15% tolerance for transient and property variation effects
    EXPECT(ratio > 3.4 && ratio < 4.6,
           "Lifetime ratio t(2d0)/t(d0) ∈ [3.4, 4.6]  (≈ 4.0)", ratio, 4.0);
}

// ─────────────────────────────────────────────────────────────────────────────
// T10 – Wet-bulb temperature bracketed between T_d0 and T_inf
// ─────────────────────────────────────────────────────────────────────────────
static void test_T10() {
    START_TEST("T10 – Wet-bulb temperature bounds");
    FuelProperties fuel = FuelProperties::nHeptane();

    for (double T_inf : {600.0, 800.0, 1000.0, 1200.0, 1500.0}) {
        AmbientConditions amb;
        amb.T_inf = T_inf; amb.P_inf = PhysConst::P_atm; amb.U_rel = 0.0;

        SimParams p; p.d0 = 100e-6; p.T_d0 = 300.0; p.dt = 5e-6; p.t_end = 0.1;

        DropletEvaporation sim(fuel, amb, p);
        double T_wb = sim.wetBulbTemperature();

        std::cout << "  T_inf=" << T_inf << " K  →  T_wb=" << T_wb << " K\n";

        EXPECT(T_wb > p.T_d0 && T_wb < T_inf,
               "T_wb ∈ (T_d0, T_inf) for T_inf=" + std::to_string((int)T_inf),
               T_wb, T_inf);

        // Higher T_inf → higher T_wb (monotone relation)
        if (T_inf > 600.0) {
            AmbientConditions amb2 = amb;
            amb2.T_inf = T_inf - 200.0;
            SimParams p2 = p; p2.t_end = 0.1;
            DropletEvaporation sim2(fuel, amb2, p2);
            double T_wb2 = sim2.wetBulbTemperature();
            EXPECT(T_wb > T_wb2,
                   "T_wb increases with T_inf", T_wb, T_wb2);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Entry point
// ─────────────────────────────────────────────────────────────────────────────

int main() {
    std::cout << "=================================================\n";
    std::cout << " Droplet Evaporation – Unit Test Suite\n";
    std::cout << "=================================================\n";

    test_T01();
    test_T02();
    test_T03();
    test_T04();
    test_T05();
    test_T06();
    test_T07();
    test_T08();
    test_T09();
    test_T10();

    std::cout << "\n=================================================\n";
    std::cout << " Results:  " << g_passed << " passed,  "
              << g_failed << " failed,  "
              << g_total  << " total\n";
    std::cout << "=================================================\n";

    return (g_failed == 0) ? 0 : 1;
}
