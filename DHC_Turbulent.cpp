#include "DHC_Turbulent.h"
#include<iostream>
#include <fstream>
#include <chrono>
#include <cmath>
#include <vector>
#include <algorithm>
#include <ranges>

using namespace std;
using namespace std::chrono;

// Function to allocate matrices
vector<vector<double>> CreateMatrix(int rows, int cols, double init_val = 0.0) {
    return vector(rows, vector(cols, init_val));
}

// Base Class to implement multiple methods
class ConvectiveScheme {
public:
    virtual void ComputeRt(double rho, int Ny, int Nx, double dx, double lambda, double cp, vector<vector<double> > &u,
                           vector<vector<double> > &v, vector<vector<double> > &Rt,
                           vector<vector<double> > &T) = 0; // Pure virtual function

    virtual void ComputeRu(vector<vector<double> > &u, vector<vector<double> > &v, vector<vector<double> > &Ru,
                           double rho, double dx, double mu, int Nx, int Ny) = 0;

    virtual void ComputeRv(vector<vector<double> > &u, vector<vector<double> > &v, vector<vector<double> > &Rv,
                           double rho, double dx, double mu, int Nx, int Ny) = 0;

    virtual void Nusselt(int Nx, int Ny, double Tcold, double Thot, double L, double alpha, double dx, double &Nux_avg,
                         vector<vector<double> > &T1star, vector<vector<double> > &T1, vector<vector<double> > &qx,
                         vector<vector<double> > &u1star, vector<vector<double> > &u1, vector<vector<double> > &dTdx,
                         vector<double> &Nux, double dxstar, double dystar) = 0;

    virtual ~ConvectiveScheme() = default; // Virtual destructor
};

// Central Differencing Scheme (CDS) method
class CDSMethod : public ConvectiveScheme {

public:
    // Function to compute R(t) on each node of stagg-T mesh (same as stagg-P)
    void ComputeRt(double rho, int Ny, int Nx, double dx, double lambda, double cp, vector<vector<double> > &u,
                   vector<vector<double> > &v, vector<vector<double> > &Rt, vector<vector<double> > &T) override {
        for (int i = 1; i < Ny-1; i++) {
            for (int j = 1; j < Nx-1; j++) {
                Rt[i][j] = -rho * dx * cp * (u[i][j + 1] * 1.0 / 2 * (T[i][j] + T[i][j + 1]) - u[i][j] * 1.0 / 2 * (
                                            T[i][j] + T[i][j - 1]) + v[i][j] * 1.0 / 2 * (T[i][j] + T[i - 1][j]) - v[
                                            i + 1][j] * 1.0 / 2 * (T[i][j] + T[i + 1][j])) + lambda * (
                               T[i][j + 1] + T[i][j - 1] + T[i - 1][j] + T[i + 1][j] - 4 * T[i][j]);
            }
        }
    }

    // Function to compute R(u) on internal nodes of stagg-x mesh
    void ComputeRu(vector<vector<double> > &u, vector<vector<double> > &v, vector<vector<double> > &Ru, double rho,
                   double dx, double mu, int Nx, int Ny) override {
        for (int i = 1; i < Ny - 1; i++) {
            for (int j = 1; j < Nx; j++) {
                Ru[i][j] = -rho * dx * (1.0 / 2 * (v[i][j - 1] + v[i][j]) * 1.0 / 2 * (u[i][j] + u[i - 1][j])
                                        - 1.0 / 2 * (v[i + 1][j - 1] + v[i + 1][j]) * 1.0 / 2 * (u[i][j] + u[i + 1][j])
                                        + 1.0 / 2 * (u[i][j] + u[i][j + 1]) * 1.0 / 2 * (u[i][j] + u[i][j + 1])
                                        - 1.0 / 2 * (u[i][j] + u[i][j - 1]) * 1.0 / 2 * (u[i][j] + u[i][j - 1]))
                           + mu * (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1] - 4 * u[i][j]);
            }
        }
    }

    // Function to compute R(v) on internal nodes of stagg-y mesh
    void ComputeRv(vector<vector<double> > &u, vector<vector<double> > &v, vector<vector<double> > &Rv, double rho,
                   double dx, double mu, int Nx, int Ny) override {
        for (int i = 1; i < Ny; i++) {
            for (int j = 1; j < Nx - 1; j++) {
                Rv[i][j] = -rho * dx * (1.0 / 2 * (v[i][j] + v[i - 1][j]) * 1.0 / 2 * (v[i][j] + v[i - 1][j])
                                        - 1.0 / 2 * (v[i][j] + v[i + 1][j]) * 1.0 / 2 * (v[i][j] + v[i + 1][j])
                                        + 1.0 / 2 * (u[i - 1][j + 1] + u[i][j + 1]) * 1.0 / 2 * (v[i][j] + v[i][j + 1])
                                        - 1.0 / 2 * (u[i - 1][j] + u[i][j]) * 1.0 / 2 * (v[i][j] + v[i][j - 1]))
                           + mu * (v[i - 1][j] + v[i + 1][j] + v[i][j - 1] + v[i][j + 1] - 4 * v[i][j]);
            }
        }
    }

    void Nusselt(int Nx, int Ny, double Tcold, double Thot, double H, double alpha, double dx, double &Nux_avg,
             vector<vector<double>> &T1star, vector<vector<double>> &T_avg, vector<vector<double>> &qx,
             vector<vector<double>> &u1star, vector<vector<double>> &u_avg, vector<vector<double>> &dTdx,
             vector<double> &Nux, double dxstar, double dystar) override {

        // Adimensionalizzazione basata su H (altezza della cavità)
        dxstar = dx / H;
        dystar = dx / H;

        // Temperatura adimensionale
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                T1star[i][j] = (T_avg[i][j] - Tcold) / (Thot - Tcold);
            }
        }

        // Velocità orizzontale adimensionale: u* = u H / α
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx + 1; j++) {
                u1star[i][j] = u_avg[i][j] * H / alpha;
            }
        }

        // Derivata e flusso sulla parete calda (x = 0)
        for (int i = 0; i < Ny; ++i) {
            dTdx[i][0] = (T1star[i][1] - T1star[i][0]) / dxstar;
            double T_central = 0.5 * (T1star[i][1] + T1star[i][0]);
            qx[i][0] = u1star[i][0] * T_central - dTdx[i][0];
        }

        // Celle interne: schema CDS
        for (int i = 0; i < Ny; ++i) {
            for (int j = 1; j < Nx; ++j) {
                dTdx[i][j] = (T1star[i][j] - T1star[i][j - 1]) / dxstar;
                double T_central = 0.5 * (T1star[i][j] + T1star[i][j - 1]);
                qx[i][j] = u1star[i][j] * T_central - dTdx[i][j];
            }
        }

        // Integrazione verticale (in y*)
        for (int j = 0; j < Nx; ++j) {
            double sum = 0.0;
            for (int i = 0; i < Ny; ++i) {
                sum += qx[i][j] * dystar;
            }
            Nux[j] = sum;
        }

        // Nusselt medio lungo la parete
        Nux_avg = 0.0;
        for (int j = 0; j < Nx; ++j) {
            Nux_avg += Nux[j];
        }
        Nux_avg /= Nx;
    }
};

// Upwind Differencing Scheme (UDS) method
class UDSMethod : public ConvectiveScheme {

public:
    // Function to compute R(t) on each node of stagg-T mesh (same as stagg-P)
    void ComputeRt(double rho, int Ny, int Nx, double dx, double lambda, double cp, vector<vector<double> > &u,
                   vector<vector<double> > &v, vector<vector<double> > &Rt, vector<vector<double> > &T) override {
        for (int i = 1; i < Ny-1; i++) {
            for (int j = 1; j < Nx-1; j++) {
                double convW = u[i][j] > 0 ? u[i][j] * T[i][j - 1] : u[i][j] * T[i][j];

                double convE = u[i][j + 1] > 0 ? u[i][j + 1] * T[i][j] : u[i][j + 1] * T[i][j + 1];

                double convS = v[i + 1][j] > 0 ? v[i + 1][j] * T[i + 1][j] : v[i + 1][j] * T[i][j];

                double convN = v[i][j] > 0 ? v[i][j] * T[i][j] : v[i][j] * T[i - 1][j];

                Rt[i][j] = -rho * dx * cp * (convE - convW + convN - convS) +
                           lambda * (T[i][j + 1] + T[i][j - 1] + T[i - 1][j] + T[i + 1][j] - 4 * T[i][j]);
            }
        }
    }

    // Function to compute R(u) on internal nodes of stagg-x mesh
    void ComputeRu(vector<vector<double> > &u, vector<vector<double> > &v,
                   vector<vector<double> > &Ru, double rho, double dx, double mu, int Nx, int Ny) override {
        for (int i = 1; i < Ny - 1; i++) {
            for (int j = 1; j < Nx; j++) {
                double v_n = (v[i][j - 1] + v[i][j]) / 2;
                double conv_n = v_n > 0 ? v_n * u[i][j] : v_n * u[i - 1][j];

                double v_s = (v[i + 1][j - 1] + v[i + 1][j]) / 2;
                double conv_s = v_s > 0 ? v_s * u[i + 1][j] : v_s * u[i][j];

                double conv_w = u[i][j] > 0 ? u[i][j] * u[i][j - 1] : u[i][j] * u[i][j];

                double conv_e = u[i][j + 1] > 0 ? u[i][j + 1] * u[i][j] : u[i][j + 1] * u[i][j + 1];

                Ru[i][j] = -rho * dx * (conv_e - conv_w + conv_n - conv_s)
                           + mu * (u[i - 1][j] + u[i + 1][j] + u[i][j - 1] + u[i][j + 1] - 4 * u[i][j]);
            }
        }
    }

    //Function to compute R(v) on internal nodes of stagg-y mesh
    void ComputeRv(vector<vector<double> > &u, vector<vector<double> > &v, vector<vector<double> > &Rv, double rho,
                   double dx, double mu, int Nx, int Ny) override {
        for (int i = 1; i < Ny; i++) {
            for (int j = 1; j < Nx - 1; j++) {
                double conv_n = v[i][j] > 0 ? v[i][j] * v[i - 1][j] : v[i][j] * v[i][j];

                double conv_s = v[i + 1][j] > 0 ? v[i + 1][j] * v[i][j] : v[i + 1][j] * v[i + 1][j];

                double u_e = 0.5 * (u[i - 1][j + 1] + u[i][j + 1]);
                double conv_e = u_e > 0 ? u_e * v[i][j] : u_e * v[i][j + 1];

                double u_w = 0.5 * (u[i - 1][j] + u[i][j]);
                double conv_w = u_w > 0 ? u_w * v[i][j - 1] : u_w * v[i][j];

                Rv[i][j] = -rho * dx * (conv_e - conv_w + conv_n - conv_s)
                           + mu * (v[i - 1][j] + v[i + 1][j] + v[i][j - 1] + v[i][j + 1] - 4 * v[i][j]);
            }
        }
    }

    // Funzione per il calcolo del Nusselt medio (UDS) con adimensionalizzazione basata su H
    void Nusselt(int Nx, int Ny, double Tcold, double Thot, double H, double alpha, double dx, double &Nux_avg,
                 vector<vector<double> > &T1star, vector<vector<double> > &T_avg, vector<vector<double> > &qx,
                 vector<vector<double> > &u1star, vector<vector<double> > &u_avg, vector<vector<double> > &dTdx,
                 vector<double> &Nux, double dxstar, double dystar) override {
        // Adimensionalizzazione su H
        dxstar = dx / H;
        dystar = dx / H;

        // Temperatura adimensionale
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                T1star[i][j] = (T_avg[i][j] - Tcold) / (Thot - Tcold);
            }
        }

        // Velocità adimensionale u* = u * H / α
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx + 1; j++) {
                u1star[i][j] = u_avg[i][j] * H / alpha;
            }
        }

        // Calcolo flusso UDS (Upwind Differencing Scheme)
        for (int i = 0; i < Ny; ++i) {
            // Parete calda x = 0
            dTdx[i][0] = (T1star[i][1] - T1star[i][0]) / dxstar;
            double T_upwind_0 = (u1star[i][0] > 0.0) ? T1star[i][0] : T1star[i][1];
            qx[i][0] = u1star[i][0] * T_upwind_0 - dTdx[i][0];

            // Seconda colonna
            dTdx[i][1] = (T1star[i][1] - T1star[i][0]) / dxstar;
            double T_upwind_1 = (u1star[i][1] > 0.0) ? T1star[i][0] : T1star[i][1];
            qx[i][1] = u1star[i][1] * T_upwind_1 - dTdx[i][1];

            // Celle interne
            for (int j = 2; j < Nx; ++j) {
                if (u1star[i][j] > 0.0) {
                    dTdx[i][j] = (T1star[i][j - 1] - T1star[i][j - 2]) / dxstar;
                } else {
                    dTdx[i][j] = (T1star[i][j] - T1star[i][j - 1]) / dxstar;
                }
                double T_upwind = (u1star[i][j] > 0.0) ? T1star[i][j - 1] : T1star[i][j];
                qx[i][j] = u1star[i][j] * T_upwind - dTdx[i][j];
            }
        }

        // Integrazione verticale (in y*)
        for (int j = 0; j < Nx; ++j) {
            double sum = 0.0;
            for (int i = 0; i < Ny; ++i) {
                sum += qx[i][j] * dystar;
            }
            Nux[j] = sum;
        }

        // Nusselt medio
        Nux_avg = 0.0;
        for (int j = 0; j < Nx; ++j) {
            Nux_avg += Nux[j];
        }
        Nux_avg /= Nx;
    }
};

// Function to choose the method
ConvectiveScheme* createCalculator(const string &method) {
    if (method == "CDS") {
        return new CDSMethod();
    } else if (method == "UDS") {
        return new UDSMethod();
    } else {
        cerr << "Not a valid method!" << endl;
        return nullptr;
    }
}

// Function to apply boundary conditions
void BoundaryConditions(int Nx, int Ny, double Thot, double Tcold, vector<vector<double> > &Tn_1,
                        vector<vector<double> > &Tn, vector<vector<double> > &T1) {
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            // LEFT WALL (x = 0): Dirichlet BC
            if (j == 0) {
                Tn_1[i][j] = Thot;
                Tn[i][j] = Thot;
                T1[i][j] = Thot;
            }
            // RIGHT WALL (x = L): Dirichlet BC
            else if (j == Nx - 1) {
                Tn_1[i][j] = Tcold;
                Tn[i][j] = Tcold;
                T1[i][j] = Tcold;
            }
        }
    }
    // TOP WALL (y = L -> i = 0): Neumann dT/dz = 0 -> T[0][j] = T[1][j]
    for (int j = 0; j < Nx; j++) {
        Tn_1[0][j] = Tn_1[1][j];
        Tn[0][j] = Tn[1][j];
        T1[0][j] = T1[1][j];
    }
    // BOTTOM WALL (y = 0 -> i = Ny - 1): Neumann dT/dz = 0 -> T[Ny-1][j] = T[Ny-2][j]
    for (int j = 0; j < Nx; j++) {
        Tn_1[Ny - 1][j] = Tn_1[Ny - 2][j];
        Tn[Ny - 1][j] = Tn[Ny - 2][j];
        T1[Ny - 1][j] = T1[Ny - 2][j];
    }
}

// Function to solve the pressure field with Gauss-Seidel solver on all nodes of stagg-P mesh
void PoissonSolver(double maxResP, double maxIteP, double rho, double dx, double dt, int Nx, int Ny, double omega_P,
                   vector<vector<double> > &P1, vector<vector<double> > &vP, vector<vector<double> > &uP,
                   vector<vector<double> > &Pg) {
    double resP = maxResP + 1;
    int iteP = 0;
    while (resP > maxResP && iteP < maxIteP) {
        double maxDiffP = 0.0;
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                if (i == 0)
                    P1[i][j] = 1.0 / 1 * (P1[i+1][j] - rho * dx / dt * (vP[i][j] - vP[i+1][j] + uP[i][j+1] - uP[i][j]));
                else if (j == 0)
                    P1[i][j] = 1.0 / 1 * (P1[i][j+1] - rho * dx / dt * (vP[i][j] - vP[i+1][j] + uP[i][j+1] - uP[i][j]));
                else if (i == Ny - 1)
                    P1[i][j] = 1.0 / 1 * (P1[i-1][j] - rho * dx / dt * (vP[i][j] - vP[i+1][j] + uP[i][j+1] - uP[i][j]));
                else if (j == Nx - 1)
                    P1[i][j] = 1.0 / 1 * (P1[i][j-1] - rho * dx / dt * (vP[i][j] - vP[i+1][j] + uP[i][j+1] - uP[i][j]));
                else
                    P1[i][j] = 1.0 / 4 * (P1[i+1][j] + P1[i][j+1] + P1[i-1][j] + P1[i][j-1]
                        - rho * dx / dt * (vP[i][j] - vP[i+1][j] + uP[i][j+1] - uP[i][j]));
                double diffP = fabs(P1[i][j] - Pg[i][j]);
                if (diffP > maxDiffP) {
                    maxDiffP = diffP;
                }
            }
        }
        resP = maxDiffP;
        // Update P guess with over-relaxation
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                if (omega_P < 1e-8) {
                    Pg[i][j] = P1[i][j];
                } else {
                    Pg[i][j] = Pg[i][j] + omega_P * (P1[i][j] - Pg[i][j]);
                }
            }
        }
        iteP++;
    }
}

int main () {
    // Start timer
    auto start = high_resolution_clock::now();

    /////// Physical data ///////

    int Nx;
    int Ny;

    // Ask for mesh refinement
    cout << "Insert number of control volumes on x direction: ";
    cin >> Nx;
    Ny = Nx * 4;
    double Ra;
    double L = 1.0;
    double H = 4.0;
    double dx = L / Nx;
    double Pr = 0.71;
    double Thot = 1.0;
    double Tcold = 0.0;
    double Tinf = 0.5;
    double cp = 0.71;
    double lambda = 1.0;
    double mu = Pr * lambda / cp;
    double rho = 1.0;
    double beta = 1.0;

    // Ask for desired Rayleigh number
    cout << "Insert Rayleigh number: ";
    cin >> Ra;
    double g = Ra * mu * lambda / (cp * beta * pow(rho, 2) * (Thot - Tcold) * pow(H, 3));   //HO SCAMBIATO L CON H
    double alpha = lambda / (rho * cp);
    double Nux_avg = 0.0;
    double Vmax;
    double dxstar;
    double dystar;
    double Ek_old = 0.0;

    /////// Numerical data ///////

    int ite = 0;
    int sample_count = 0;
    double maxResP = 1e-3;
    double maxIte = 1e10;
    double maxResT = 1e-3;
    double maxResU = maxResT;
    double maxResV = maxResT;
    double res1 = maxResU + 1;
    double res2 = maxResV + 1;
    double res3 = maxResT + 1;
    double res3_old = res3;
    double t_count = 0.0;
    double t_count_avg = 0.0;
    double total_GS_time = 0.0;
    // Starting time for averaging
    double t_0 = 0.6;
    // Total duration of the simulation
    const double t_min = 2.5;
    double dt = 1e-6;
    double dtc, dtd, dtt;
    double omega_P = 1.8;
    bool switched = false;
    string method;

    // Ask for convective scheme
    cout << "Choose the method ('CDS', 'UDS', or 'AUTO'): ";
    cin >> method;
    ranges::transform(method, method.begin(), ::toupper);
    ConvectiveScheme* calculator = nullptr;
    ConvectiveScheme* calculator_UDS = nullptr;
    ConvectiveScheme* calculator_CDS = nullptr;
    if (method == "CDS" || method == "UDS") {
        calculator = createCalculator(method);
        if (!calculator) {
            cerr << "Error: Invalid method selected." << endl;
            return 1;
        }
        cout << "Running with fixed " << method << " scheme.\n";
    } else if (method == "AUTO") {
        calculator_UDS = createCalculator("UDS");
        calculator_CDS = createCalculator("CDS");
        if (!calculator_UDS || !calculator_CDS) {
            cerr << "Error: Could not create convection schemes." << endl;
            return 1;
        }
        calculator = calculator_UDS; // start with UDS
        cout << "AUTO mode: starting with UDS, will switch to CDS based on residuals.\n";
    } else {
        cerr << "Error: Method must be 'CDS', 'UDS', or 'AUTO'." << endl;
        return 1;
    }

    /////// Definition of matrices and vectors ///////

    // Pressure, Temperature, u and v
    auto P1    = CreateMatrix(Ny, Nx);
    auto Pg    = CreateMatrix(Ny, Nx);
    auto Tn_1 = CreateMatrix(Ny, Nx);
    auto Tn = CreateMatrix(Ny, Nx);
    auto T1 = CreateMatrix(Ny, Nx);
    auto un_1  = CreateMatrix(Ny, Nx + 1);
    auto un    = CreateMatrix(Ny, Nx + 1);
    auto u1    = CreateMatrix(Ny, Nx + 1);
    auto uP    = CreateMatrix(Ny, Nx + 1);
    auto vn_1  = CreateMatrix(Ny + 1, Nx);
    auto vn    = CreateMatrix(Ny + 1, Nx);
    auto v1    = CreateMatrix(Ny + 1, Nx);
    auto vP    = CreateMatrix(Ny + 1, Nx);

    // Convective-diffusive, Buoyancy term R() matrices
    auto Ru1    = CreateMatrix(Ny, Nx + 1);
    auto Run    = CreateMatrix(Ny, Nx + 1);
    auto Run_1  = CreateMatrix(Ny, Nx + 1);
    auto Rv1    = CreateMatrix(Ny + 1, Nx);
    auto Rvn    = CreateMatrix(Ny + 1, Nx);
    auto Rvn_1  = CreateMatrix(Ny + 1, Nx);
    auto Rt1 = CreateMatrix(Ny, Nx);
    auto Rtn = CreateMatrix(Ny, Nx);
    auto Rtn_1 = CreateMatrix(Ny, Nx);
    auto Rbn = CreateMatrix(Ny + 1, Nx);
    auto Rbn_1 = CreateMatrix(Ny + 1, Nx);

    // Nusselt number
    auto T1star = CreateMatrix(Ny, Nx);
    auto u1star = CreateMatrix(Ny, Nx + 1);
    auto v1star = CreateMatrix(Ny + 1, Nx);
    auto dTdx = CreateMatrix(Ny, Nx);
    auto qx = CreateMatrix(Ny, Nx);
    auto V = CreateMatrix(Ny, Nx);
    vector Nux(Ny, 0.0);
    vector u1star_L2(Ny, 0.0);
    vector v1star_L2(Ny, 0.0);

    // Averaged fields
    auto T_avg = CreateMatrix(Ny, Nx, 0.0);
    auto u_avg = CreateMatrix(Ny, Nx + 1, 0.0);
    auto v_avg = CreateMatrix(Ny + 1, Nx, 0.0);
    auto uu_avg = CreateMatrix(Ny, Nx + 1, 0.0);
    auto vv_avg = CreateMatrix(Ny + 1, Nx, 0.0);
    auto uv_avg = CreateMatrix(Ny, Nx, 0.0);
    auto psi = CreateMatrix(Ny, Nx, 0.0);

    // Kinetic Energy
    auto u_var = CreateMatrix(Ny, Nx + 1, 0.0);
    auto v_var = CreateMatrix(Ny + 1, Nx, 0.0);
    auto uv_fluct = CreateMatrix(Ny, Nx, 0.0);
    auto ke = CreateMatrix(Ny, Nx, 0.0);
    auto tke = CreateMatrix(Ny, Nx, 0.0);
    auto u_interp_sum = CreateMatrix(Ny, Nx, 0.0);
    auto uu_interp_sum = CreateMatrix(Ny, Nx, 0.0);
    auto v_interp_sum = CreateMatrix(Ny + 1, Nx, 0.0);
    auto vv_interp_sum = CreateMatrix(Ny + 1, Nx, 0.0);
    auto uv_interp_sum = CreateMatrix(Ny + 1, Nx, 0.0);

    /////// Probes analysis ///////

    int i_probe_bottom_left = Ny - 2;
    int j_probe_bottom_left = 1;
    int i_probe_bottom_right = Ny - 2;
    int j_probe_bottom_right = Nx - 2;
    int i_probe_top_left = 1;
    int j_probe_top_left = 1;
    int i_probe_top_right = 1;
    int j_probe_top_right = Nx - 2;

    /////// Initialization ///////

    // Initial velocity fields = 0.0
    for (int i = 1; i < Ny; i++) {
        for (int j = 1; j < Nx + 1; j++) {
            u1[i][j] = 0.0;
            un[i][j] = u1[i][j];
            un_1[i][j] = un[i][j];
        }
    }
    for (int i = 1; i < Ny + 1; i++) {
        for (int j = 1; j < Nx; j++) {
            v1[i][j] = 0.0;
            vn[i][j] = v1[i][j];
            vn_1[i][j] = vn[i][j];
        }
    }
    // Initial pressure and linearly distributed temperature field
    for (int i = 0; i < Ny; i++) {
        for (int j = Nx-1; j >= 0; --j) {
            Pg[i][j] = 0.0;
            T1[i][j] = 0;
            Tn[i][j] = T1[i][j];
            Tn_1[i][j] = Tn[i][j];
        }
    }

    // Apply boundary conditions
    BoundaryConditions(Nx, Ny, Thot, Tcold, Tn_1, Tn, T1);

    // Compute R()^n-1
    calculator->ComputeRu(un_1, vn_1, Run_1, rho, dx, mu, Nx, Ny);
    calculator->ComputeRv(un_1, vn_1, Rvn_1, rho, dx, mu, Nx, Ny);
    calculator->ComputeRt(rho, Ny, Nx, dx, lambda, cp, un_1, vn_1, Rtn_1, Tn_1);

    // Compute R()^n
    calculator->ComputeRu(un, vn, Run, rho, dx, mu, Nx, Ny);
    calculator->ComputeRv(un, vn, Rvn, rho, dx, mu, Nx, Ny);
    calculator->ComputeRt(rho, Ny, Nx, dx, lambda, cp, un, vn, Rtn, Tn);

    // Open Info file
    ofstream Infos ("Info_Ra_" + to_string(static_cast<int>(Ra)) + ".txt");
    Infos << "==== SIMULATION INFO ====" << endl;
    Infos << "RAYLEIGH = " << Ra << endl;
    Infos << "Mesh size = " << Nx << " x " << Ny << endl;
    Infos << "Convective scheme = " << method << endl;
    Infos << "Gauss-Seidel tollerance = " << maxResP << endl;
    Infos << "Relaxation Poisson = " << omega_P << endl;
    Infos << "Start averaging time  = " << t_0 << " s" << endl;
    Infos << "Averaging window      = " << t_min - t_0 << " s" << endl;
    Infos << "Total simulation time = " << t_min << " s" << endl;

    // Open probes and kinetic budget script files
    ofstream probes_out("ProbesHistory_Ra_" + to_string(static_cast<int>(Ra)) + ".txt");
    ofstream energy_out("EnergyTerms_Ra_" + to_string(static_cast<int>(Ra)) + ".txt");

    /////// Time loop ///////

    while (t_count <= t_min && ite < maxIte) {
        double maxDiff1 = 0.0;
        double maxDiff2 = 0.0;
        double maxDiff3 = 0.0;

        // Step 1a
        for (int i = 1; i < Ny - 1; i++) {
            for (int j = 1; j < Nx; j++) {
                uP[i][j] = un[i][j] + dt / (rho * pow(dx, 2)) * (3.0 / 2 * Run[i][j] - 1.0 / 2 * Run_1[i][j]);
            }
        }

        // Step 1b + Boussinesq approx.
        for (int i = 1; i < Ny; i++) {
            for (int j = 1; j < Nx - 1; j++) {
                Rbn[i][j] = rho * pow(dx, 2) * beta * (Tn[i][j] - Tinf) * g;
                Rbn_1[i][j] = rho * pow(dx, 2) * beta * (Tn_1[i][j] - Tinf) * g;
                vP[i][j] = vn[i][j] + dt / (rho * pow(dx, 2)) * (
                               3.0 / 2 * Rvn[i][j] - 1.0 / 2 * Rvn_1[i][j] + 3.0 / 2 * Rbn[i][j] - 1.0 / 2 * Rbn_1[i][
                                   j]);
            }
        }

        // Step 2
        if (t_count > t_0) {
            maxResP = 1e-4;
        }
        auto start_GS = std::chrono::high_resolution_clock::now();
        PoissonSolver(maxResP, maxIte, rho, dx, dt, Nx, Ny, omega_P, P1, vP, uP, Pg);
        auto end_GS = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> dt_GS = end_GS - start_GS;
        total_GS_time += dt_GS.count();

        // Step 3
        for (int i = 1; i < Ny - 1; i++) {
            for (int j = 1; j < Nx; j++) {
                u1[i][j] = uP[i][j] - dt / (rho * dx) * (P1[i][j] - P1[i][j - 1]);
                double diff1 = fabs(u1[i][j] - un[i][j]);
                if (diff1 > maxDiff1) {
                    maxDiff1 = diff1;
                }
            }
        }
        for (int i = 1; i < Ny; i++) {
            for (int j = 1; j < Nx - 1; j++) {
                v1[i][j] = vP[i][j] - dt / (rho * dx) * (P1[i - 1][j] - P1[i][j]);
                double diff2 = fabs(v1[i][j] - vn[i][j]);
                if (diff2 > maxDiff2) {
                    maxDiff2 = diff2;
                }
            }
        }

        // Step 4 - Temperature evaluation
        for (int i = 1; i < Ny - 1; i++) {
            for (int j = 1; j < Nx - 1; j++) {
                T1[i][j] = Tn[i][j] + dt / (pow(dx, 2) * rho * cp) * (3.0 / 2 * Rtn[i][j] - 1.0 / 2 * Rtn_1[i][j]);
                double diff3 = fabs(T1[i][j] - Tn[i][j]);
                if (diff3 > maxDiff3) {
                    maxDiff3 = diff3;
                }
            }
        }
        // Update temperature boundary conditions
        for (int j = 0; j < Nx; j++) {
            T1[0][j] = T1[1][j];
            T1[Ny - 1][j] = T1[Ny - 2][j];
        }

        // Step 5: Non-dimensional kinetic energy balance

        double Ek = 0.0;
        double error = 0.0;
        double R_total = 0.0;
        double R_diff = 0.0;
        double R_conv = 0.0;
        double R_press = 0.0;
        double R_buoy = 0.0;

        // Reference energy scale: rho * alpha^2
        double Ek_ref = rho * alpha * alpha;

        for (int i = 1; i < Ny-1; i++) {
            for (int j = 1; j < Nx-1; j++) {
                // Interpolate u and v to cell center (dimensional)
                double u_c = 0.5 * (u1[i][j] + u1[i][j + 1]);
                double v_c = 0.5 * (v1[i][j] + v1[i + 1][j]);

                // Kinetic energy (adimensionalized at the end)
                Ek += 0.5 * (u_c * u_c + v_c * v_c) * dx * dx;

                // Viscous diffusion term
                double lap_u = (u1[i][j + 1] + u1[i][j - 1] + u1[i - 1][j] + u1[i + 1][j] - 4.0 * u1[i][j]) / (dx * dx);
                double lap_v = (v1[i][j + 1] + v1[i][j - 1] + v1[i - 1][j] + v1[i + 1][j] - 4.0 * v1[i][j]) / (dx * dx);
                R_diff += mu * (u_c * lap_u + v_c * lap_v) * dx * dx;

                // Convective transport term
                double dudx = (u1[i][j + 1] - u1[i][j - 1]) / (2.0 * dx);
                double dudy = (u1[i - 1][j] - u1[i + 1][j]) / (2.0 * dx);
                double dvdx = (v1[i][j + 1] - v1[i][j - 1]) / (2.0 * dx);
                double dvdy = (v1[i - 1][j] - v1[i + 1][j]) / (2.0 * dx);
                double conv_u = u_c * dudx + v_c * dudy;
                double conv_v = u_c * dvdx + v_c * dvdy;
                R_conv -= rho * (u_c * conv_u + v_c * conv_v) * dx * dx;

                // Pressure work
                double dpdx = (P1[i][j + 1] - P1[i][j - 1]) / (2.0 * dx);
                double dpdy = (P1[i - 1][j] - P1[i + 1][j]) / (2.0 * dx);
                R_press -= (u_c * dpdx + v_c * dpdy) * dx * dx;

                // Buoyancy term
                double Fb = rho * beta * (T1[i][j] - Tinf) * g;
                R_buoy += Fb * v_c * dx * dx;
            }
        }

        // Final nondimensionalization
        Ek /= Ek_ref;
        R_diff /= Ek_ref;
        R_conv /= Ek_ref;
        R_press /= Ek_ref;
        R_buoy /= Ek_ref;

        // Time derivative (dimensionless time)
        double dt_star = dt * alpha / (H * H);
        double dEk_num;
        if (t_count > 0) {
            dEk_num = (Ek - Ek_old) / dt_star;
        } else {
            dEk_num = 0.0;
        }

        // Total residual and mismatch
        R_total = R_diff + R_conv + R_press + R_buoy;
        error = fabs(R_total - dEk_num);

        // Output (already nondimensional)
        energy_out << t_count << " "
                << Ek << " "
                << R_diff << " "
                << R_conv << " "
                << R_press << " "
                << R_buoy << " "
                << R_total << " "
                << dEk_num << " "
                << error << endl;

        // Update old energy
        Ek_old = Ek;

        // Update time counter and residual
        t_count += dt;
        res1 = maxDiff1;
        res2 = maxDiff2;
        res3 = maxDiff3;

        // Switch method criteria (when 'AUTO' is selected)
        if (method == "AUTO" && !switched && res3 < 1e-3 && fabs(res3 - res3_old) / res3_old < 1e-2) {
            calculator = calculator_CDS;
            switched = true;
            cout << "Switched from UDS to CDS at iteration " << ite << " (res3 = " << res3 << ")\n";
        }

        // Update temperature residual
        res3_old = res3;

        // Update u^n and u^n-1
        for (int i = 1; i < Ny - 1; i++) {
            for (int j = 1; j < Nx; j++) {
                un_1[i][j] = un[i][j];
                un[i][j] = u1[i][j];
            }
        }

        // Update v^n and v^n-1
        for (int i = 1; i < Ny; i++) {
            for (int j = 1; j < Nx - 1; j++) {
                vn_1[i][j] = vn[i][j];
                vn[i][j] = v1[i][j];
            }
        }

        // Update T^n and T^n-1
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                Tn_1[i][j] = Tn[i][j];
                Tn[i][j] = T1[i][j];
            }
        }

        // Compute R()^n+1
        calculator->ComputeRu(u1, v1, Ru1, rho, dx, mu, Nx, Ny);
        calculator->ComputeRv(u1, v1, Rv1, rho, dx, mu, Nx, Ny);
        calculator->ComputeRt(rho, Ny, Nx, dx, lambda, cp, u1, v1, Rt1, T1);

        // Update R()^n and R()^n-1
        for (int i = 1; i < Ny - 1; i++) {
            for (int j = 1; j < Nx; j++) {
                Run_1[i][j] = Run[i][j];
                Run[i][j] = Ru1[i][j];
            }
        }
        for (int i = 1; i < Ny; i++) {
            for (int j = 1; j < Nx - 1; j++) {
                Rvn_1[i][j] = Rvn[i][j];
                Rvn[i][j] = Rv1[i][j];
            }
        }
        for (int i = 1; i < Ny-1; i++) {
            for (int j = 1; j < Nx-1; j++) {
                Rtn_1[i][j] = Rtn[i][j];
                Rtn[i][j] = Rt1[i][j];
            }
        }

        // Print probes history
        probes_out << t_count << " "
                   << T1[i_probe_bottom_left][j_probe_bottom_left] << " " << u1[i_probe_bottom_left][j_probe_bottom_left] << " " << v1[i_probe_bottom_left][j_probe_bottom_left] << " "
                   << T1[i_probe_bottom_right][j_probe_bottom_right] << " " << u1[i_probe_bottom_right][j_probe_bottom_right] << " " << v1[i_probe_bottom_right][j_probe_bottom_right] << " "
                   << T1[i_probe_top_left][j_probe_top_left] << " " << u1[i_probe_top_left][j_probe_top_left] << " " << v1[i_probe_top_left][j_probe_top_left] << " "
                   << T1[i_probe_top_right][j_probe_top_right] << " " << u1[i_probe_top_right][j_probe_top_right] << " " << v1[i_probe_top_right][j_probe_top_right]
                   << endl;

        // Sum of the fields for averaging
        if (t_count > t_0) {
            for (int i = 0; i < Ny; i++) {
                for (int j = 0; j < Nx; j++) {
                    T_avg[i][j] += T1[i][j] * dt;

                    u_avg[i][j] += 0.5 * (u1[i][j] + u1[i][j + 1]) * dt;
                    v_avg[i][j] += 0.5 * (v1[i][j] + v1[i + 1][j]) * dt;

                    uu_avg[i][j] += pow(0.5 * (u1[i][j] + u1[i][j + 1]), 2) * dt;
                    vv_avg[i][j] += pow(0.5 * (v1[i][j] + v1[i + 1][j]), 2) * dt;
                    uv_avg[i][j] += 0.5 * (u1[i][j] + u1[i][j + 1]) * 0.5 * (v1[i][j] + v1[i + 1][j]) * dt;
                }
            }
            t_count_avg += dt;
            sample_count++;
        }

        //Update time step
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                V[i][j] = sqrt(pow(u1[i][j], 2) + pow(v1[i][j], 2));
            }
        }
        Vmax = std::ranges::max(std::views::join(V));
        dtc = 0.35 * dx * L / Vmax;
        dtd = 0.20 * pow(dx, 2) / (mu / rho);
        dtt = 0.20 * pow(dx, 2) / (lambda / (rho * cp));
        dt = min({dtc, dtd, dtt});

        // Go to next time step and print values of interest for debugging purposes
        ite++;
        if (ite % 100 == 0) {
            cout << "SUMMARY: ite = " << ite
                    << " | resT = " << res3
                    << " | resU = " << res1
                    << " | resV = " << res2
                    << " | t = " << t_count
                    << " | dt = " << dt
                    << endl;
        }
    }

    // Close probes history and kinetic budget files
    probes_out.close();
    energy_out.close();

    // Average Fields
    if (sample_count > 0) {
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                T_avg[i][j] /= t_count_avg;
                u_avg[i][j] /= t_count_avg;
                v_avg[i][j] /= t_count_avg;
                uu_avg[i][j] /= t_count_avg;
                vv_avg[i][j] /= t_count_avg;
                uv_avg[i][j] /= t_count_avg;
            }
        }
        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                u_var[i][j] = uu_avg[i][j] - pow(u_avg[i][j], 2);
                v_var[i][j] = vv_avg[i][j] - pow(v_avg[i][j], 2);
                uv_fluct[i][j] = uv_avg[i][j] - u_avg[i][j] * v_avg[i][j];
                ke[i][j] = 0.5 * (pow(u_avg[i][j],2) + pow(v_avg[i][j], 2));
                tke[i][j] = 0.5 * (u_var[i][j] + v_var[i][j]);
            }
        }
    }
    Infos << "Fields averaged at " << t_count << " seconds" << endl;

    // Stream Function Calculation
    psi[0][0] = 0.0;
    // Vertical integration in the first column
    for (int i = 1; i < Ny; ++i) {
        psi[i][0] = psi[i - 1][0] + dx * u_avg[i - 1][0];
    }
    // Horizontal integration in the rows
    for (int i = 0; i < Ny; ++i) {
        for (int j = 1; j < Nx; ++j) {
            psi[i][j] = psi[i][j - 1] - dx * v_avg[i][j - 1];
        }
    }

    /////// .txt files for plotting ///////

    // Averaged temperature distribution
    ofstream TemperatureDistribution("TemperatureDistribution_Ra_" + to_string(static_cast<int>(Ra)) + ".txt");
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            TemperatureDistribution << j * dx + dx / 2 << " " << H - i * dx - dx / 2 << " " << T_avg[i][j] << endl;
        }
        TemperatureDistribution << "\n";
    }

    // Averaged velocity u distribution
    ofstream VelocityUDistribution("VelocityUDistribution_Ra_" + to_string(static_cast<int>(Ra)) + ".txt");
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx + 1; j++) {
            VelocityUDistribution << j * dx + dx / 2 << " " << H - i * dx - dx / 2 << " " << u_avg[i][j] << endl;
        }
        VelocityUDistribution << "\n";
    }

    // Averaged velocity v distribution
    ofstream VelocityVDistribution("VelocityVDistribution_Ra_" + to_string(static_cast<int>(Ra)) + ".txt");
    for (int i = 0; i < Ny + 1; i++) {
        for (int j = 0; j < Nx; j++) {
            VelocityVDistribution << j * dx + dx / 2 << " " << H - i * dx - dx / 2 << " " << v_avg[i][j] << endl;
        }
        VelocityVDistribution << "\n";
    }

    // u(y) at L/2
    ofstream Uyl2("UyL2_Ra_" + to_string(static_cast<int>(Ra)) + ".txt");
    Uyl2 << 4 << " " << 0 << endl;
    for (int i = 0; i < Ny; i++) {
        Uyl2 << H - i * dx - dx / 2 << " " << (u_avg[i][Nx / 2] + u_avg[i][Nx / 2 + 1]) / 2 << endl;
    }
    Uyl2 << 0 << " " << 0 << endl;

    // v(x) at L/2
    ofstream Vxl2("VxL2_Ra_" + to_string(static_cast<int>(Ra)) + ".txt");
    Vxl2 << 0 << " " << 0 << endl;
    for (int j = 0; j < Nx; j++) {
        Vxl2 << j * dx + dx / 2 << " " << (v_avg[Ny / 2][j] + v_avg[Ny / 2 + 1][j]) / 2 << endl;
    }
    Vxl2 << 1 << " " << 0 << endl;

    // Stream Function
    ofstream StreamFunction("StreamFunction_Ra_" + to_string(static_cast<int>(Ra)) + ".txt");
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            double x = j * dx + dx / 2;
            double y = H - i * dx - dx / 2;
            StreamFunction << x << " " << y << " " << psi[i][j] << endl;
        }
        StreamFunction << "\n";
    }

    // Reynolds stresses, turbulent kinetic energy and mean kinetic energy
    ofstream Turbulent_kinetic("Turbulent_KE_" + to_string(static_cast<int>(Ra)) + ".txt");
    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            double x = j * dx + dx / 2;
            double y = H - i * dx - dx / 2;
            Turbulent_kinetic << x << " " << y << " " << u_var[i][j] << " " << v_var[i][j] << " " << uv_fluct[i][j] <<
                    " " << tke[i][j] << " " << ke[i][j] << endl;
        }
        Turbulent_kinetic << "\n";
    }

    // Adding Boundaries
    // BORDO SUPERIORE (i = -1)
    for (int j = 0; j < Nx; j++) {
        double x = j * dx + dx / 2;
        double y = H;
        TemperatureDistribution << x << " " << y << " " << T_avg[0][j] << endl;
        VelocityUDistribution << x << " " << y << " " << "0" << endl;
        VelocityVDistribution << x << " " << y << " " << "0" << endl;
        StreamFunction << x << " " << y << " " << "0" << endl;
        Turbulent_kinetic << x << " " << y << " " << "0" << " " << "0" << " " << "0" << " " <<
                    "0" << " " << "0" << endl;
    }

    // BORDO INFERIORE (i = Ny)
    for (int j = 0; j < Nx; j++) {
        double x = j * dx + dx / 2;
        double y = 0.0;
        TemperatureDistribution << x << " " << y << " " <<  T_avg[Ny-1][j] << endl;
        VelocityUDistribution << x << " " << y << " " << "0" << endl;
        VelocityVDistribution << x << " " << y << " " << "0" << endl;
        StreamFunction << x << " " << y << " " << "0" << endl;
        Turbulent_kinetic << x << " " << y << " " << "0" << " " << "0" << " " << "0" << " " <<
                    "0" << " " << "0" << endl;
    }

    // BORDO SINISTRO (j = -1)
    for (int i = 0; i < Ny; i++) {
        double x = 0.0;
        double y = H - i * dx - dx / 2;
        TemperatureDistribution << x << " " << y << " " << "1" << endl;
        VelocityUDistribution << x << " " << y << " " << "0" << endl;
        VelocityVDistribution << x << " " << y << " " << "0" << endl;
        StreamFunction << x << " " << y << " " << "0" << endl;
        Turbulent_kinetic << x << " " << y << " " << "0" << " " << "0" << " " << "0" << " " <<
                    "0" << " " << "0" << endl;
    }

    // BORDO DESTRO (j = Nx)
    for (int i = 0; i < Ny; i++) {
        double x = 1.0;
        double y = H - i * dx - dx / 2;
        TemperatureDistribution << x << " " << y << " " << "0" << endl;
        VelocityUDistribution << x << " " << y << " " << "0" << endl;
        VelocityVDistribution << x << " " << y << " " << "0" << endl;
        StreamFunction << x << " " << y << " " << "0" << endl;
        Turbulent_kinetic << x << " " << y << " " << "0" << " " << "0" << " " << "0" << " " <<
                    "0" << " " << "0" << endl;
    }

    // ANGOLO ALTO SINISTRO
    TemperatureDistribution << 0.0 << " " << H << " " << "1" << endl;
    VelocityUDistribution << 0.0 << " " << H << " " << "0" << endl;
    VelocityVDistribution << 0.0 << " " << H << " " << "0" << endl;
    StreamFunction << 0.0 << " " << H << " " << "0" << endl;
    Turbulent_kinetic << 0.0 << " " << H << " " << "0" << " " << "0" << " " << "0" << " " <<
                    "0" << " " << "0" << endl;
    // ANGOLO ALTO DESTRO
    TemperatureDistribution << 1.0 << " " << H << " " << "0" << endl;
    VelocityUDistribution << 1.0 << " " << H << " " << "0" << endl;
    VelocityVDistribution << 1.0 << " " << H << " " << "0" << endl;
    StreamFunction << 1.0 << " " << H << " " << "0" << endl;
    Turbulent_kinetic << 1.0 << " " << H << " " << "0" << " " << "0" << " " << "0" << " " <<
                    "0" << " " << "0" << endl;
    // ANGOLO BASSO SINISTRO
    TemperatureDistribution << 0.0 << " " << 0.0 << " " << "1" << endl;
    VelocityUDistribution << 0.0 << " " << 0.0 << " " << "0" << endl;
    VelocityVDistribution << 0.0 << " " << 0.0 << " " << "0" << endl;
    StreamFunction << 0.0 << " " << 0.0 << " " << "0" << endl;
    Turbulent_kinetic << 0.0 << " " << 0.0 << " " << "0" << " " << "0" << " " << "0" << " " <<
                    "0" << " " << "0" << endl;
    // ANGOLO BASSO DESTRO
    TemperatureDistribution << 1.0 << " " << 0.0 << " " << "0" << endl;
    VelocityUDistribution << 1.0 << " " << 0.0 << " " << "0" << endl;
    VelocityVDistribution << 1.0 << " " << 0.0 << " " << "0" << endl;
    StreamFunction << 1.0 << " " << 0.0 << " " << "0" << endl;
    Turbulent_kinetic << 1.0 << " " << 0.0 << " " << "0" << " " << "0" << " " << "0" << " " <<
                    "0" << " " << "0" << endl;


    TemperatureDistribution.close();
    VelocityUDistribution.close();
    VelocityVDistribution.close();
    Uyl2.close();
    Vxl2.close();
    StreamFunction.close();
    Turbulent_kinetic.close();

    /////// Final calculations ///////

    // Average Nusselt number
    calculator->Nusselt(Nx, Ny, Tcold, Thot, L, alpha, dx, Nux_avg, T1star, T1, qx, u1star, u1, dTdx, Nux, dxstar,
                        dystar);
    Infos << "Average Nusselt for Ra = " << Ra << " --> " << Nux_avg << endl;

    // Maximum u* velocity at x = L/2, adimensionalizzata con H
    for (int i = 0; i < Ny; i++) {
        u1star_L2[i] = (u_avg[i][Nx / 2] + u_avg[i][Nx / 2 + 1]) / 2 * H / alpha;
    }
    auto ustar_max = ranges::max_element(u1star_L2);
    int u_index = distance(u1star_L2.begin(), ustar_max);
    double y_u = (Ny - u_index - 0.5) * dx;  // y fisico
    double ystar_u = y_u / H;  // y adimensionale con H
    Infos << "Max u* velocity = " << *ustar_max << " at y/H = " << ystar_u << endl;

    // Maximum v* velocity at y = H/2, adimensionalizzata con H
    for (int j = 0; j < Nx; j++) {
        v1star_L2[j] = (v_avg[Ny / 2][j] + v_avg[Ny / 2 + 1][j]) / 2 * H / alpha;
    }
    auto vstar_max = ranges::max_element(v1star_L2);
    int v_index = distance(v1star_L2.begin(), vstar_max);
    double x_v = (v_index + 0.5) * dx;  // x fisico
    double xstar_v = x_v / H;  // x adimensionale con H (NON più /L!)
    Infos << "Max v* velocity = " << *vstar_max << " at x/H = " << xstar_v << endl;

    // Stop timer and print total duration
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<seconds>(stop - start);
    Infos << "Total Execution Time: " << duration.count() << " seconds" << endl;
    Infos << "Total GS Solver Time: " << total_GS_time << " seconds = " << total_GS_time * 100 / duration.count() <<
            "% of the total execution time" << endl;

    Infos.close();

    if (method == "AUTO") {
        delete calculator_UDS;
        delete calculator_CDS;
    } else {
        delete calculator;
    }

    return 0;
}