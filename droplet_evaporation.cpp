/**
 * droplet_evaporation.cpp
 *
 * Implementation of the droplet heating and evaporation model.
 * See droplet_evaporation.h for full documentation.
 */

#include "droplet_evaporation.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

using namespace PhysConst;

// ─────────────────────────────────────────────────────────────────────────────
// FuelProperties
// ─────────────────────────────────────────────────────────────────────────────

FuelProperties FuelProperties::nHeptane() {
    FuelProperties f;
    f.name      = "n-Heptane (C7H16)";
    f.MW_F      = 0.100204;     // kg/mol
    f.rho_l     = 684.0;        // kg/m³  (at 298 K)
    f.c_l       = 2250.0;       // J/(kg·K)
    f.L_v       = 317000.0;     // J/kg   (at boiling point 371.6 K)
    f.T_boil    = 371.6;        // K
    f.T_crit    = 540.2;        // K
    f.cp_vapor  = 1600.0;       // J/(kg·K)  (approximate)
    // Antoine coefficients (NIST, T in °C, P in mmHg):
    f.A_Antoine = 6.89386;
    f.B_Antoine = 1264.02;
    f.C_Antoine = 216.640;
    return f;
}

double FuelProperties::saturationPressure(double T_K) const {
    // Clamp to physical range
    double T_C = T_K - 273.15;
    T_C = std::min(T_C, T_crit - 273.15 - 1.0);
    T_C = std::max(T_C, -100.0);

    double log10_P_mmHg = A_Antoine - B_Antoine / (C_Antoine + T_C);
    double P_mmHg       = std::pow(10.0, log10_P_mmHg);
    return P_mmHg * 133.322; // mmHg → Pa
}

// ─────────────────────────────────────────────────────────────────────────────
// GasProperties
// ─────────────────────────────────────────────────────────────────────────────

GasProperties GasProperties::evaluate(double T_film, double P_Pa, double MW_air) {
    GasProperties g;
    g.T_film = T_film;

    // --- Viscosity: Sutherland's law for air ---
    const double T_ref_S = 273.15, mu_ref = 1.716e-5, S_const = 110.4;
    g.mu_g = mu_ref * std::pow(T_film / T_ref_S, 1.5)
             * (T_ref_S + S_const) / (T_film + S_const);

    // --- Thermal conductivity: power law for air ---
    g.k_g = 0.0241 * std::pow(T_film / 273.15, 0.82);

    // --- Specific heat: 4th-order polynomial fit for air [J/(kg·K)] ---
    // Valid ~200 K – 2000 K  (Gordon & McBride coefficients, simplified)
    const double T = T_film;
    g.cp_g = 1030.5 - 0.19975 * T + 3.9734e-4 * T * T
             - 2.319e-7  * T * T * T + 4.940e-11 * T * T * T * T;
    g.cp_g = std::max(g.cp_g, 1000.0);  // physical lower bound

    // --- Density: ideal gas for air ---
    g.rho_g = P_Pa * MW_air / (R_universal * T_film);

    // --- Binary diffusivity: Chapman-Enskog power-law scaling ---
    // D_AB ~ T^1.75 / P  (Fuller et al., 1966)
    // Reference: D(n-heptane / air) ≈ 7.7e-6 m²/s at 300 K, 1 atm
    g.D_AB = 7.7e-6 * std::pow(T_film / 300.0, 1.75)
             / (P_Pa / P_atm);

    // --- Dimensionless numbers ---
    g.Pr = g.mu_g * g.cp_g / g.k_g;
    g.Sc = g.mu_g / (g.rho_g * g.D_AB);

    return g;
}

// ─────────────────────────────────────────────────────────────────────────────
// DropletEvaporation – constructor
// ─────────────────────────────────────────────────────────────────────────────

DropletEvaporation::DropletEvaporation(const FuelProperties&   fuel,
                                       const AmbientConditions& ambient,
                                       const SimParams&         params)
    : fuel_(fuel), ambient_(ambient), params_(params)
{
    if (params_.d0     <= 0.0) throw std::invalid_argument("d0 must be > 0");
    if (params_.dt     <= 0.0) throw std::invalid_argument("dt must be > 0");
    if (params_.t_end  <= 0.0) throw std::invalid_argument("t_end must be > 0");
    if (params_.T_d0   <= 0.0) throw std::invalid_argument("T_d0 must be > 0");
    if (ambient_.T_inf <= 0.0) throw std::invalid_argument("T_inf must be > 0");
}

// ─────────────────────────────────────────────────────────────────────────────
// Core: compute derivatives  dd/dt  and  dT/dt
// ─────────────────────────────────────────────────────────────────────────────

DropletEvaporation::State2
DropletEvaporation::derivatives(const State2& s, DropletState* diag) const {
    const double d   = s.d;
    const double T_d = s.T;

    // Guard: fully-evaporated droplet – return zero derivatives to stop integration
    if (d <= 0.0) {
        if (diag) {
            diag->Re = diag->Nu = diag->Sh = 0.0;
            diag->B_M = diag->mdot = diag->Q_conv = diag->K_evap = 0.0;
        }
        return {0.0, 0.0};
    }

    // Film-rule temperature (1/3 rule gives better accuracy, use simple 1/2)
    double T_film = 0.5 * (T_d + ambient_.T_inf);
    T_film = std::max(T_film, 200.0);

    GasProperties g = GasProperties::evaluate(T_film, ambient_.P_inf,
                                               ambient_.MW_air);

    // ── Reynolds number ──────────────────────────────────────────────────────
    double Re = (ambient_.U_rel > 0.0)
                ? g.rho_g * ambient_.U_rel * d / g.mu_g
                : 0.0;

    // ── Spalding mass transfer number  B_M ───────────────────────────────────
    double P_sat = fuel_.saturationPressure(T_d);
    // Mole fraction at surface
    double X_Fs = std::min(P_sat / ambient_.P_inf, 0.9999);
    // Convert to mass fraction (fuel / air binary mixture)
    double ratio = X_Fs * fuel_.MW_F / ((1.0 - X_Fs) * ambient_.MW_air);
    double Y_Fs  = ratio / (1.0 + ratio);
    Y_Fs = std::max(Y_Fs, ambient_.Y_F_inf);   // avoid B_M < 0

    double B_M = (Y_Fs - ambient_.Y_F_inf) / (1.0 - Y_Fs);
    B_M = std::max(B_M, 1e-12);

    // ── Ranz-Marshall  Nu₀ and Sh₀ ───────────────────────────────────────────
    double Re05 = std::sqrt(Re);
    double Nu0  = 2.0 + 0.6 * Re05 * std::pow(g.Pr, 1.0 / 3.0);
    double Sh0  = 2.0 + 0.6 * Re05 * std::pow(g.Sc, 1.0 / 3.0);

    // ── Spalding blowing correction (Abramzon-Sirignano) ─────────────────────
    // F_M corrects Sh for mass blowing at the surface
    double F_M = std::log(1.0 + B_M) / B_M;
    double Sh  = 2.0 + (Sh0 - 2.0) * F_M;

    // Heat-transfer Spalding number B_T (relates to B_M via Le number)
    // Simple approximation: B_T ≈ (1+B_M)^(cp_vapor/cp_g / Le) - 1
    // Here we use B_T = cp_vapor*(T_inf - T_d) / L_v  (1st-law estimate)
    double B_T  = fuel_.cp_vapor * (ambient_.T_inf - T_d) / fuel_.L_v;
    B_T = std::max(B_T, 1e-12);
    double F_T  = std::log(1.0 + B_T) / B_T;
    double Nu   = 2.0 + (Nu0 - 2.0) * F_T;

    // ── Evaporation rate  ṁ  [kg/s]  (negative = mass loss) ─────────────────
    double mdot = -PI * d * g.rho_g * g.D_AB * Sh * std::log(1.0 + B_M);

    // ── Convective heat transfer  Q_conv  [W] ────────────────────────────────
    double Q_conv = PI * d * g.k_g * Nu * (ambient_.T_inf - T_d);

    // ── Droplet mass and specific heat ───────────────────────────────────────
    double m_d  = fuel_.rho_l * PI / 6.0 * d * d * d;

    // ── Energy equation: m*c_l * dT/dt = Q_conv + ṁ*L_v ────────────────────
    // (mdot < 0, L_v > 0  →  mdot*L_v < 0 cools the droplet)
    double dTdt = (Q_conv + mdot * fuel_.L_v) / (m_d * fuel_.c_l);

    // ── Diameter ODE:  dm/dt = ρ_l * (π/2) * d² * (dd/dt) ──────────────────
    double dddt = mdot / (fuel_.rho_l * PI / 2.0 * d * d);

    // ── Optionally fill diagnostic record ────────────────────────────────────
    if (diag) {
        diag->Re     = Re;
        diag->Nu     = Nu;
        diag->Sh     = Sh;
        diag->B_M    = B_M;
        diag->mdot   = mdot;
        diag->Q_conv = Q_conv;
        diag->K_evap = -dddt * 2.0 * d;  // d(d²)/dt = 2d*(dd/dt)
    }

    return {dddt, dTdt};
}

// ─────────────────────────────────────────────────────────────────────────────
// Build a full DropletState record
// ─────────────────────────────────────────────────────────────────────────────

DropletState DropletEvaporation::makeRecord(double t, const State2& s) const {
    DropletState rec;
    rec.time             = t;
    rec.diameter         = s.d;
    rec.temperature      = s.T;
    rec.mass             = fuel_.rho_l * PI / 6.0 * s.d * s.d * s.d;
    rec.diameter_sq_norm = (s.d / params_.d0) * (s.d / params_.d0);
    derivatives(s, &rec);   // fills Re, Nu, Sh, B_M, mdot, Q_conv, K_evap
    return rec;
}

// ─────────────────────────────────────────────────────────────────────────────
// Analytical D²-law evaporation constant  K  [m²/s]
// ─────────────────────────────────────────────────────────────────────────────

double DropletEvaporation::D2LawConstant() const {
    // Quasi-steady D²-law constant evaluated at the wet-bulb temperature,
    // which is where the droplet spends most of its lifetime.
    // K = 8 * rho_g * D_AB * ln(1 + B_M) / rho_l
    double T_wb  = wetBulbTemperature();
    double T_film = 0.5 * (T_wb + ambient_.T_inf);
    GasProperties g = GasProperties::evaluate(T_film, ambient_.P_inf,
                                               ambient_.MW_air);
    double B_M_qs;
    {
        double P_sat = fuel_.saturationPressure(T_wb);
        double X_Fs  = std::min(P_sat / ambient_.P_inf, 0.9999);
        double ratio = X_Fs * fuel_.MW_F / ((1.0 - X_Fs) * ambient_.MW_air);
        double Y_Fs  = ratio / (1.0 + ratio);
        B_M_qs       = std::max((Y_Fs - ambient_.Y_F_inf) / (1.0 - Y_Fs), 1e-12);
    }
    return 8.0 * g.rho_g * g.D_AB * std::log(1.0 + B_M_qs) / fuel_.rho_l;
}

// ─────────────────────────────────────────────────────────────────────────────
// Wet-bulb temperature  (iterative energy balance)
// ─────────────────────────────────────────────────────────────────────────────

double DropletEvaporation::wetBulbTemperature() const {
    // Wet-bulb: dT/dt = 0  →  Q_conv + ṁ * L_v = 0
    // Solve:  Nu*(T_inf - T_wb) = Sh * (L_v / cp_gas) * ln(1+B_M)
    // Simple bisection on T_wb ∈ [T_d0, T_inf]
    const double tol = 1e-4;
    const int    max_iter = 200;

    auto residual = [&](double T_wb) -> double {
        State2 s{params_.d0, T_wb};
        DropletState diag;
        derivatives(s, &diag);
        // At steady state: dT/dt = 0  →  Q_conv + mdot*L_v = 0
        return diag.Q_conv + diag.mdot * fuel_.L_v;
    };

    double T_lo = params_.T_d0 + 1.0;
    // Physical upper bound: at 1 atm the droplet cannot exceed the boiling
    // point (higher T_inf → droplet boils, T_wb → T_boil).
    double T_hi = std::min(ambient_.T_inf - 1.0, fuel_.T_boil - 0.1);
    if (T_hi <= T_lo) return fuel_.T_boil;

    // Verify bracketing; if both endpoints have the same sign the droplet
    // is in the boiling regime for this ambient condition.
    if (residual(T_lo) * residual(T_hi) > 0.0) {
        return fuel_.T_boil;
    }

    double T_mid = 0.5 * (T_lo + T_hi);
    for (int i = 0; i < max_iter; ++i) {
        T_mid = 0.5 * (T_lo + T_hi);
        double r_mid = residual(T_mid);
        if (std::abs(T_hi - T_lo) < tol) break;
        if (r_mid * residual(T_lo) < 0.0)
            T_hi = T_mid;
        else
            T_lo = T_mid;
    }
    return T_mid;
}

// ─────────────────────────────────────────────────────────────────────────────
// Main simulation loop  (4th-order Runge-Kutta)
// ─────────────────────────────────────────────────────────────────────────────

std::vector<DropletState> DropletEvaporation::run() {
    const double dt      = params_.dt;
    const double t_end   = params_.t_end;
    const double d_min   = 1.0e-9;   // 1 nm – consider fully evaporated

    // ── Temperature freeze threshold ─────────────────────────────────────────
    // Near T_boil, dP_sat/dT is steep, so dB_M/dT_d is large.  This makes
    // the energy ODE stiff: effectively dT/dt ~ S/d² where S is a stiffness
    // coefficient (~65 W/(m·K) for n-heptane near 1000 K).  Explicit RK4 with
    // step dt is stable only when d > d_stiff = sqrt(S * dt / (rho_l*pi/6*c_l)).
    // We capture T at d_capture (slightly above d_stiff for safety) and freeze
    // it thereafter so the temperature report is physically meaningful.
    const double S_stiff   = 70.0;  // W/(m·K), conservative stiffness bound
    const double d_stiff   = std::sqrt(S_stiff * dt
                                       / (fuel_.rho_l * PhysConst::PI / 6.0
                                          * fuel_.c_l));
    // Also ensure the droplet has had time to heat up (at least 5% of d0 evaporated)
    // Use 1.5× safety factor; also ensure at least 10% of d0 has evaporated
    // so the heating transient is underway before we freeze T.
    const double d_capture = std::max(1.5 * d_stiff, 0.10 * params_.d0);
    double T_frozen    = params_.T_d0;
    bool   T_is_frozen = false;

    // Determine save stride
    int stride = params_.save_stride;
    if (stride <= 0) {
        // Default: ~1000 output points over the full run
        stride = std::max(1, static_cast<int>(t_end / dt / 1000));
    }

    std::vector<DropletState> history;
    history.reserve(1001);

    State2 s{params_.d0, params_.T_d0};
    double t = 0.0;

    // Save initial state
    history.push_back(makeRecord(t, s));

    // ── RK4 lambda ───────────────────────────────────────────────────────────
    // Intermediate diameters are clamped to >= 0 to prevent singularities;
    // derivatives() returns {0,0} when d <= 0, so the integration degrades
    // gracefully as the droplet approaches complete evaporation.
    auto rk4_step = [&](const State2& s_in) -> State2 {
        State2 k1 = derivatives(s_in);

        State2 s2{std::max(s_in.d + 0.5 * dt * k1.d, 0.0),
                  std::clamp(s_in.T + 0.5 * dt * k1.T, 200.0, ambient_.T_inf)};
        State2 k2 = derivatives(s2);

        State2 s3{std::max(s_in.d + 0.5 * dt * k2.d, 0.0),
                  std::clamp(s_in.T + 0.5 * dt * k2.T, 200.0, ambient_.T_inf)};
        State2 k3 = derivatives(s3);

        State2 s4{std::max(s_in.d + dt * k3.d, 0.0),
                  std::clamp(s_in.T + dt * k3.T, 200.0, ambient_.T_inf)};
        State2 k4 = derivatives(s4);

        return {
            s_in.d + dt / 6.0 * (k1.d + 2.0*k2.d + 2.0*k3.d + k4.d),
            s_in.T + dt / 6.0 * (k1.T + 2.0*k2.T + 2.0*k3.T + k4.T)
        };
    };

    // ── Integration loop ─────────────────────────────────────────────────────
    for (int step = 1; t + dt <= t_end + 0.5 * dt; ++step) {
        // Terminate BEFORE calling rk4_step on a fully-evaporated droplet
        if (!std::isfinite(s.d) || s.d < d_min) {
            if (history.empty() || history.back().time < t - 0.5 * dt)
                history.push_back(makeRecord(t, s));
            break;
        }

        State2 s_new = rk4_step(s);

        // Physical bounds
        s_new.d = std::max(s_new.d, 0.0);
        s_new.T = std::clamp(s_new.T, 200.0, ambient_.T_inf);

        // Capture quasi-steady temperature the first time d drops below
        // d_capture; freeze it thereafter to avoid the 1/d² singularity
        // in the energy equation that appears as d → 0.
        if (!T_is_frozen && s_new.d < d_capture) {
            T_frozen    = s.T;   // T at last stable size
            T_is_frozen = true;
        }
        if (T_is_frozen) s_new.T = T_frozen;

        t += dt;
        s  = s_new;

        // Record output
        if (step % stride == 0 && std::isfinite(s.d)) {
            history.push_back(makeRecord(t, s));
        }

        // Check full evaporation at end of step
        if (!std::isfinite(s.d) || s.d < d_min) {
            // Clamp to zero for a clean final record
            s.d = 0.0;
            history.push_back(makeRecord(t, s));
            break;
        }
    }

    return history;
}
