/**
 * droplet_evaporation.h
 *
 * Droplet Heating and Evaporation Simulation
 * ==========================================
 * Models a single liquid fuel droplet (n-heptane, C7H16) evaporating in a hot
 * quiescent or convective air environment.
 *
 * Physical Model
 * --------------
 * - Quasi-steady gas-phase assumption
 * - Uniform droplet temperature (infinite liquid thermal conductivity)
 * - Temperature-dependent gas-phase properties (Sutherland's law, power laws)
 * - Ranz-Marshall correlations for Nusselt and Sherwood numbers
 * - Spalding blowing correction factors
 * - Antoine equation for fuel vapor pressure
 * - 4th-order Runge-Kutta time integration
 *
 * Governing Equations
 * -------------------
 *  Mass:    dm/dt = m_dot = -pi * d * rho_g * D_AB * Sh * ln(1 + B_M)
 *  Energy:  m * c_l * dT_d/dt = Q_conv + m_dot * L_v
 *           Q_conv = pi * d * k_g * Nu * (T_inf - T_d)
 *  Diameter: d(d^2)/dt = -K,  K = 8*rho_g*D_AB*ln(1+B_M)/rho_l  (D^2-law)
 *
 * References
 * ----------
 *  [1] Ranz & Marshall, Chem. Eng. Prog., 48, 141-146 (1952)
 *  [2] Spalding, ARS Journal, 29, 828-835 (1959)
 *  [3] Law, Prog. Energy Combust. Sci., 8, 171-201 (1982)
 *  [4] Abramzon & Sirignano, Int. J. Heat Mass Transfer, 32, 1605-1618 (1989)
 */

#pragma once

#include <vector>
#include <string>
#include <functional>
#include <cmath>

// ─────────────────────────────────────────────────────────────────────────────
// Physical constants
// ─────────────────────────────────────────────────────────────────────────────
namespace PhysConst {
    constexpr double PI          = 3.14159265358979323846;
    constexpr double R_universal = 8.314;          // J/(mol·K)
    constexpr double P_atm       = 101325.0;       // Pa
}

// ─────────────────────────────────────────────────────────────────────────────
// Fuel / liquid-phase properties  (n-heptane, C7H16)
// ─────────────────────────────────────────────────────────────────────────────
struct FuelProperties {
    std::string name;
    double MW_F;        // Molecular weight [kg/mol]
    double rho_l;       // Liquid density [kg/m³]
    double c_l;         // Liquid specific heat [J/(kg·K)]
    double L_v;         // Latent heat of vaporization [J/kg]
    double T_boil;      // Normal boiling temperature [K]
    double T_crit;      // Critical temperature [K]
    double cp_vapor;    // Vapor specific heat [J/(kg·K)]
    // Antoine equation: log10(P_vap [mmHg]) = A - B / (C + T [°C])
    double A_Antoine;
    double B_Antoine;
    double C_Antoine;

    /** Returns saturation pressure [Pa] at temperature T_K [K]. */
    double saturationPressure(double T_K) const;

    /** Default: n-heptane */
    static FuelProperties nHeptane();
};

// ─────────────────────────────────────────────────────────────────────────────
// Ambient / gas-phase conditions
// ─────────────────────────────────────────────────────────────────────────────
struct AmbientConditions {
    double T_inf;       // Far-field gas temperature [K]
    double P_inf;       // Ambient pressure [Pa]
    double U_rel;       // Droplet–gas relative velocity [m/s]
    double Y_F_inf;     // Fuel mass fraction in far field (0 = pure air)
    double MW_air;      // Molecular weight of air [kg/mol]

    AmbientConditions()
        : T_inf(1000.0), P_inf(PhysConst::P_atm),
          U_rel(0.0),    Y_F_inf(0.0), MW_air(0.028964) {}
};

// ─────────────────────────────────────────────────────────────────────────────
// Instantaneous gas-phase transport properties (evaluated at film temperature)
// ─────────────────────────────────────────────────────────────────────────────
struct GasProperties {
    double T_film;  // Film temperature [K]
    double k_g;     // Thermal conductivity [W/(m·K)]
    double mu_g;    // Dynamic viscosity [Pa·s]
    double rho_g;   // Density [kg/m³]
    double cp_g;    // Specific heat [J/(kg·K)]
    double D_AB;    // Binary diffusivity fuel-air [m²/s]
    double Pr;      // Prandtl number
    double Sc;      // Schmidt number

    /**
     * Evaluate all gas properties at the given film temperature and pressure.
     * Uses Sutherland's law for viscosity, power-law for conductivity and
     * specific heat, and Chapman-Enskog scaling for diffusivity.
     */
    static GasProperties evaluate(double T_film, double P_Pa, double MW_air);
};

// ─────────────────────────────────────────────────────────────────────────────
// Time-resolved droplet state (output record)
// ─────────────────────────────────────────────────────────────────────────────
struct DropletState {
    double time;        // Simulation time [s]
    double diameter;    // Droplet diameter [m]
    double temperature; // Droplet temperature [K]
    double mass;        // Droplet mass [kg]
    double diameter_sq_norm; // (d/d0)²  – for D²-law plots
    double Re;          // Reynolds number
    double Nu;          // Nusselt number (with blowing correction)
    double Sh;          // Sherwood number (with blowing correction)
    double B_M;         // Spalding mass transfer number
    double mdot;        // Evaporation rate [kg/s]  (< 0 means mass loss)
    double Q_conv;      // Convective heat rate to droplet [W]
    double K_evap;      // Instantaneous D²-law evaporation constant [m²/s]
};

// ─────────────────────────────────────────────────────────────────────────────
// Simulation parameters
// ─────────────────────────────────────────────────────────────────────────────
struct SimParams {
    double d0;          // Initial diameter [m]
    double T_d0;        // Initial droplet temperature [K]
    double dt;          // Time step [s]
    double t_end;       // Maximum simulation time [s]
    int    save_stride; // Save state every N integration steps (0 = auto)

    SimParams()
        : d0(100e-6), T_d0(300.0), dt(1e-5), t_end(1.0), save_stride(0) {}
};

// ─────────────────────────────────────────────────────────────────────────────
// Main simulation class
// ─────────────────────────────────────────────────────────────────────────────
class DropletEvaporation {
public:
    DropletEvaporation(const FuelProperties&  fuel,
                       const AmbientConditions& ambient,
                       const SimParams&         params);

    /**
     * Run the simulation and return time-history of DropletState.
     * Integration stops when d < 1 nm (complete evaporation) or t > t_end.
     */
    std::vector<DropletState> run();

    /** Evaporation constant K [m²/s] predicted by the classical D²-law
     *  (evaluated at the initial conditions for a quick estimate). */
    double D2LawConstant() const;

    /** Wet-bulb temperature [K] (iterative solution). */
    double wetBulbTemperature() const;

    // Read-back accessors
    const FuelProperties&   fuel()    const { return fuel_;    }
    const AmbientConditions& ambient() const { return ambient_; }
    const SimParams&         params()  const { return params_;  }

private:
    FuelProperties   fuel_;
    AmbientConditions ambient_;
    SimParams        params_;

    // Internal 2-DOF state: [diameter, temperature]
    struct State2 { double d, T; };

    /**
     * Compute time derivatives (dd/dt, dT/dt) given current State2 and
     * optionally return the instantaneous diagnostics.
     */
    State2 derivatives(const State2& s,
                       DropletState* diag = nullptr) const;

    /** Build a DropletState record from a State2 plus diagnostics. */
    DropletState makeRecord(double t, const State2& s) const;
};
