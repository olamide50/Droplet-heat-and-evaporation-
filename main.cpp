/**
 * main.cpp
 *
 * Driver program for the droplet heating and evaporation simulation.
 * Runs three demonstration cases:
 *   Case 1 – Stagnant environment (U_rel = 0), hot air at 1000 K
 *   Case 2 – Convective environment (U_rel = 1 m/s), hot air at 1000 K
 *   Case 3 – Parametric study: varying ambient temperature (600–1500 K)
 *
 * Output is written to CSV files for post-processing / plotting.
 */

#include "droplet_evaporation.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

static void writeCSV(const std::string& filename,
                     const std::vector<DropletState>& history,
                     double /*d0*/) {
    std::ofstream out(filename);
    if (!out.is_open()) {
        std::cerr << "ERROR: cannot open " << filename << "\n";
        return;
    }
    out << std::fixed << std::setprecision(8);
    out << "time_s,diameter_m,diameter_um,temperature_K,"
           "mass_kg,d2_norm,Re,Nu,Sh,B_M,mdot_kg_s,Q_conv_W,K_evap_m2_s\n";
    for (const auto& s : history) {
        out << s.time        << ','
            << s.diameter    << ','
            << s.diameter * 1e6 << ','   // µm
            << s.temperature << ','
            << s.mass        << ','
            << s.diameter_sq_norm << ','
            << s.Re          << ','
            << s.Nu          << ','
            << s.Sh          << ','
            << s.B_M         << ','
            << s.mdot        << ','
            << s.Q_conv      << ','
            << s.K_evap      << '\n';
    }
    std::cout << "  Wrote " << history.size() << " rows to " << filename << "\n";
}

static void printSummary(const std::string& label,
                         const std::vector<DropletState>& h,
                         double K_analytical,
                         double T_wb) {
    if (h.empty()) return;
    const auto& first = h.front();
    const auto& last  = h.back();

    std::cout << "\n=== " << label << " ===\n";
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "  Initial diameter : " << first.diameter * 1e6    << " µm\n";
    std::cout << "  Initial temp     : " << first.temperature        << " K\n";
    std::cout << "  Wet-bulb temp    : " << T_wb                     << " K\n";
    std::cout << "  Analytical K_D2  : " << K_analytical * 1e6      << " mm²/s  ("
              << std::scientific << std::setprecision(3) << K_analytical
              << std::fixed << std::setprecision(4) << " m²/s)\n";

    double t_evap = last.time;
    double d_final = last.diameter * 1e6;
    std::cout << "  End time         : " << t_evap                   << " s\n";
    std::cout << "  Final diameter   : " << d_final                  << " µm\n";
    std::cout << "  Final temp       : " << last.temperature         << " K\n";

    // Estimate numerical lifetime from D²-law intercept
    // Find time where d²/d0² ≈ 0 by linear fit over last third of history
    if (h.size() > 10) {
        const auto& p1 = h[h.size() * 2 / 3];
        const auto& p2 = h.back();
        double slope = (p2.diameter_sq_norm - p1.diameter_sq_norm)
                      / (p2.time - p1.time);
        if (slope < 0.0) {
            double t_complete = p1.time - p1.diameter_sq_norm / slope;
            std::cout << "  Extrapolated lifetime : " << t_complete << " s\n";
        }
    }
    // Compare K at midpoint
    const auto& mid = h[h.size() / 2];
    std::cout << "  Numerical K_D2 (mid): " << mid.K_evap * 1e6     << " mm²/s\n";
    std::cout << "  B_M (initial)       : " << first.B_M            << "\n";
    std::cout << "  Nu  (initial)       : " << first.Nu             << "\n";
    std::cout << "  Sh  (initial)       : " << first.Sh             << "\n";
    std::cout << "  Peak Q_conv [W]     : ";
    double peak_Q = 0.0;
    for (const auto& s : h) peak_Q = std::max(peak_Q, std::abs(s.Q_conv));
    std::cout << peak_Q << "\n";
}

// ─────────────────────────────────────────────────────────────────────────────
// Main
// ─────────────────────────────────────────────────────────────────────────────

int main() {
    std::cout << "=================================================\n";
    std::cout << " Droplet Heating and Evaporation Simulation\n";
    std::cout << " Fuel: n-Heptane (C7H16) in air\n";
    std::cout << "=================================================\n";

    FuelProperties fuel = FuelProperties::nHeptane();
    std::cout << "\nFuel : " << fuel.name << "\n";
    std::cout << "  MW  = " << fuel.MW_F   * 1000 << " g/mol\n";
    std::cout << "  rho = " << fuel.rho_l         << " kg/m³\n";
    std::cout << "  L_v = " << fuel.L_v   / 1000  << " kJ/kg\n";
    std::cout << "  T_b = " << fuel.T_boil         << " K\n";

    // ─── Case 1: Stagnant, T_inf = 1000 K ─────────────────────────────────
    {
        AmbientConditions amb;
        amb.T_inf = 1000.0;
        amb.P_inf = PhysConst::P_atm;
        amb.U_rel = 0.0;

        SimParams p;
        p.d0          = 100e-6;  // 100 µm
        p.T_d0        = 300.0;   // 300 K
        p.dt          = 5e-6;    // 5 µs
        p.t_end       = 0.015;   // slightly beyond expected lifetime (~11 ms)
        p.save_stride = 5;       // save every 5 steps → ~500 points for smooth plots

        DropletEvaporation sim(fuel, amb, p);
        double K   = sim.D2LawConstant();
        double T_wb = sim.wetBulbTemperature();
        auto   h   = sim.run();

        printSummary("Case 1 – Stagnant (U=0 m/s), T_inf=1000 K", h, K, T_wb);
        writeCSV("results_case1.csv", h, p.d0);
    }

    // ─── Case 2: Convective flow, T_inf = 1000 K ──────────────────────────
    {
        AmbientConditions amb;
        amb.T_inf = 1000.0;
        amb.P_inf = PhysConst::P_atm;
        amb.U_rel = 1.0;    // 1 m/s relative velocity

        SimParams p;
        p.d0          = 100e-6;
        p.T_d0        = 300.0;
        p.dt          = 5e-6;
        p.t_end       = 0.015;
        p.save_stride = 5;

        DropletEvaporation sim(fuel, amb, p);
        double K   = sim.D2LawConstant();
        double T_wb = sim.wetBulbTemperature();
        auto   h   = sim.run();

        printSummary("Case 2 – Convective (U=1 m/s), T_inf=1000 K", h, K, T_wb);
        writeCSV("results_case2.csv", h, p.d0);
    }

    // ─── Case 3: Parametric – ambient temperature sweep ───────────────────
    {
        std::cout << "\n=== Case 3 – Parametric T_inf sweep (stagnant) ===\n";
        std::cout << std::fixed << std::setprecision(3);
        std::cout << std::setw(12) << "T_inf [K]"
                  << std::setw(14) << "T_wb [K]"
                  << std::setw(18) << "K_D2 [mm²/s]"
                  << std::setw(16) << "t_life [s]\n";
        std::cout << std::string(60, '-') << "\n";

        std::ofstream param_out("results_parametric.csv");
        param_out << "T_inf_K,T_wb_K,K_D2_m2s,t_evap_s,d_final_um\n";

        std::vector<double> T_list = {600, 700, 800, 900, 1000,
                                      1100, 1200, 1300, 1400, 1500};
        for (double T_inf : T_list) {
            AmbientConditions amb;
            amb.T_inf = T_inf;
            amb.P_inf = PhysConst::P_atm;
            amb.U_rel = 0.0;

            SimParams p;
            p.d0    = 100e-6;
            p.T_d0  = 300.0;
            p.dt    = 5e-6;
            p.t_end = 2.0;  // generous upper bound

            DropletEvaporation sim(fuel, amb, p);
            double K   = sim.D2LawConstant();
            double T_wb = sim.wetBulbTemperature();
            auto   h   = sim.run();

            const auto& last = h.back();

            // Estimate total lifetime via linear extrapolation of d²(t)
            double t_life = last.time;
            if (last.diameter > 1e-8 && h.size() > 5) {
                const auto& pa = h[h.size() - 5];
                const auto& pb = h.back();
                double slope = (pb.diameter_sq_norm - pa.diameter_sq_norm)
                              / (pb.time - pa.time);
                if (slope < -1e-20)
                    t_life = pa.time - pa.diameter_sq_norm / slope;
            }

            std::cout << std::setw(12) << T_inf
                      << std::setw(14) << T_wb
                      << std::setw(18) << K * 1e6
                      << std::setw(16) << t_life << "\n";

            param_out << std::fixed << std::setprecision(6)
                      << T_inf << ',' << T_wb << ',' << K << ','
                      << t_life << ',' << last.diameter * 1e6 << '\n';
        }
        std::cout << "  Results saved to results_parametric.csv\n";
    }

    // ─── Summary ──────────────────────────────────────────────────────────
    std::cout << "\n=================================================\n";
    std::cout << " Simulation complete.  Output files:\n";
    std::cout << "   results_case1.csv       – stagnant case\n";
    std::cout << "   results_case2.csv       – convective case\n";
    std::cout << "   results_parametric.csv  – T_inf parametric sweep\n";
    std::cout << "=================================================\n";
    return 0;
}
