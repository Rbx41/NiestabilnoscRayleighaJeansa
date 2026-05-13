#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <random>
#include <iomanip>

using namespace std;

inline int idx(int y, int x, int Columns) {
    return y * Columns + x;
}

// ======================================================
// CIC MASS DEPOSITION
// ======================================================

void aktualizuj_gestosc(
    vector<double>& rho,
    const vector<vector<double>>& gwiazdy,
    int Rows,
    int Columns,
    double masa_SMBH,
    double masa_gwiazdy)
{
    fill(rho.begin(), rho.end(), 0.0);

    // Central black hole
    rho[idx(Rows / 2, Columns / 2, Columns)] += masa_SMBH;

    for (const auto& g : gwiazdy) {

        double x = g[0];
        double y = g[2];

        int x0 = floor(x);
        int y0 = floor(y);

        double dx = x - x0;
        double dy = y - y0;

        // prevent out-of-bounds
        if (x0 < 1 || x0 >= Columns - 2 ||
            y0 < 1 || y0 >= Rows - 2)
            continue;

        // CIC weights
        double w00 = (1.0 - dx) * (1.0 - dy);
        double w10 = dx * (1.0 - dy);
        double w01 = (1.0 - dx) * dy;
        double w11 = dx * dy;

        rho[idx(y0,     x0,     Columns)] += masa_gwiazdy * w00;
        rho[idx(y0,     x0 + 1, Columns)] += masa_gwiazdy * w10;
        rho[idx(y0 + 1, x0,     Columns)] += masa_gwiazdy * w01;
        rho[idx(y0 + 1, x0 + 1, Columns)] += masa_gwiazdy * w11;
    }
}

void wygladz_gestosc(vector<double>& rho, int Rows, int Columns) {
    vector<double> temp = rho;
    
    // 3x3 smoothing kernel
    for (int y = 1; y < Rows - 1; y++) {
        for (int x = 1; x < Columns - 1; x++) {
            rho[idx(y, x, Columns)] = 
                0.25 * temp[idx(y, x, Columns)] + 
                0.125 * (temp[idx(y-1, x, Columns)] + temp[idx(y+1, x, Columns)] + 
                         temp[idx(y, x-1, Columns)] + temp[idx(y, x+1, Columns)]) +
                0.0625 * (temp[idx(y-1, x-1, Columns)] + temp[idx(y-1, x+1, Columns)] + 
                          temp[idx(y+1, x-1, Columns)] + temp[idx(y+1, x+1, Columns)]);
        }
    }
}


// ======================================================
// INITIAL GALAXY
// ======================================================

void inicjalizuj_galaktyke(
    vector<vector<double>>& gwiazdy,
    int N,
    int Rows,
    int Columns,
    double masa_SMBH)
{
    double cx = Columns / 2.0;
    double cy = Rows / 2.0;

    mt19937 gen(42);

    uniform_real_distribution<> dist_u(0.0, 1.0);
    uniform_real_distribution<> dist_theta(0.0, 2.0 * M_PI);

    double Rmax = 80.0;
    double G = 1.0;

    for (int i = 0; i < N; i++) {

        // Proper disk sampling
        double u = dist_u(gen);
        double r = Rmax * sqrt(u);

        double theta = dist_theta(gen);

        double x = cx + r * cos(theta);
        double y = cy + r * sin(theta);

        // Keplerian rotation
        //double v = sqrt(G * masa_SMBH / (r + 0.5));
        double v = 8* sqrt( G * masa_SMBH / (r + 1.0));

        // small velocity dispersion
        //v *= (0.95 + 0.1 * dist_u(gen));

        double vx = -v * sin(theta);
        double vy =  v * cos(theta);

        gwiazdy.push_back({x, vx, y, vy});
    }
}


void aktualizuj_brzegi(
    vector<double>& phi,
    const vector<double>& rho,
    int Rows,
    int Columns,
    double delta)
{
    double G = 1.0;

   
    auto calc_phi_at = [&](int y, int x) {
        double pot = 0.0;
        
       
        for (int r = 1; r < Rows - 1; ++r) {
            for (int c = 1; c < Columns - 1; ++c) {
                double mass = rho[idx(r, c, Columns)];
                
       
                if (mass > 0.0) { 
                    double dx = x - c;
                    double dy = y - r;
                    double dist = sqrt(dx * dx + dy * dy);
                    
                    if (dist > 0) {
                        pot -= mass / dist;
                    }
                }
            }
        }

        return pot * (G / delta);
    };

    // Calculate for Top and Bottom edges
    for (int x = 0; x < Columns; ++x) {
        phi[idx(0, x, Columns)] = calc_phi_at(0, x);
        phi[idx(Rows - 1, x, Columns)] = calc_phi_at(Rows - 1, x);
    }

    // Calculate for Left and Right edges (excluding corners already done)
    for (int y = 1; y < Rows - 1; ++y) {
        phi[idx(y, 0, Columns)] = calc_phi_at(y, 0);
        phi[idx(y, Columns - 1, Columns)] = calc_phi_at(y, Columns - 1);
    }
}



// ======================================================
// SOR POISSON SOLVER
// ======================================================

void nadrelaksacja(
    vector<double>& phi,
    const vector<double>& rho,
    int Rows,
    int Columns,
    double delta)
{
    double lambda = 1.85;

    // FIXED: delta^2
    double factor = 4.0 * M_PI * delta * delta;

    for (int y = 1; y < Rows - 1; y++) {

        for (int x = 1; x < Columns - 1; x++) {

            int id = idx(y, x, Columns);

            double neighbors =
                phi[idx(y + 1, x, Columns)] +
                phi[idx(y - 1, x, Columns)] +
                phi[idx(y, x + 1, Columns)] +
                phi[idx(y, x - 1, Columns)];

            double phi_new =
                0.25 * (neighbors - factor * rho[id]);

            phi[id] =
                (1.0 - lambda) * phi[id]
                + lambda * phi_new;
        }
    }
}

// ======================================================
// FORCE FIELD
// ======================================================

double finite_phi_x_difference(
    const vector<double>& phi,
    int x,
    int y,
    double delta,
    int Columns)
{
    return (
        phi[idx(y, x + 1, Columns)]
        - phi[idx(y, x - 1, Columns)]
    ) / (2.0 * delta);
}

double finite_phi_y_difference(
    const vector<double>& phi,
    int x,
    int y,
    double delta,
    int Rows,
    int Columns)
{
    return (
        phi[idx(y + 1, x, Columns)]
        - phi[idx(y - 1, x, Columns)]
    ) / (2.0 * delta);
}

// ======================================================
// CIC FORCE INTERPOLATION
// ======================================================

double interpolate_grad(
    const vector<double>& phi,
    double x,
    double y,
    double delta,
    int Rows,
    int Columns,
    bool x_grad)
{
    int x0 = floor(x);
    int y0 = floor(y);

    int x1 = x0 + 1;
    int y1 = y0 + 1;

    // safety
    if (x0 < 1) x0 = 1;
    if (y0 < 1) y0 = 1;

    if (x1 > Columns - 2) x1 = Columns - 2;
    if (y1 > Rows - 2) y1 = Rows - 2;

    double dx = x - x0;
    double dy = y - y0;

    double f00, f10, f01, f11;

    if (x_grad) {

        f00 = finite_phi_x_difference(phi, x0, y0, delta, Columns);
        f10 = finite_phi_x_difference(phi, x1, y0, delta, Columns);
        f01 = finite_phi_x_difference(phi, x0, y1, delta, Columns);
        f11 = finite_phi_x_difference(phi, x1, y1, delta, Columns);

    } else {

        f00 = finite_phi_y_difference(phi, x0, y0, delta, Rows, Columns);
        f10 = finite_phi_y_difference(phi, x1, y0, delta, Rows, Columns);
        f01 = finite_phi_y_difference(phi, x0, y1, delta, Rows, Columns);
        f11 = finite_phi_y_difference(phi, x1, y1, delta, Rows, Columns);
    }

    double w00 = (1.0 - dx) * (1.0 - dy);
    double w10 = dx * (1.0 - dy);
    double w01 = (1.0 - dx) * dy;
    double w11 = dx * dy;

    return
        w00 * f00 +
        w10 * f10 +
        w01 * f01 +
        w11 * f11;
}

// ======================================================
// EQUATIONS OF MOTION
// ======================================================

void dydx(
    double y[],
    double dydx_out[],
    const vector<double>& phi,
    int Rows,
    int Columns,
    double delta)
{
    double x = y[0];
    double ypos = y[2];

    // clamp
    if (x < 2.0) x = 2.0;
    if (x > Columns - 3.0) x = Columns - 3.0;

    if (ypos < 2.0) ypos = 2.0;
    if (ypos > Rows - 3.0) ypos = Rows - 3.0;

    dydx_out[0] = y[1];

    dydx_out[1] =
        -interpolate_grad(
            phi,
            x,
            ypos,
            delta,
            Rows,
            Columns,
            true);

    dydx_out[2] = y[3];

    dydx_out[3] =
        -interpolate_grad(
            phi,
            x,
            ypos,
            delta,
            Rows,
            Columns,
            false);
}

// ======================================================
// RK4
// ======================================================

void rungeKutta(
    double h,
    double y[],
    const vector<double>& phi,
    int Rows,
    int Columns,
    double delta)
{
    const int N = 4;

    double k1[N], k2[N], k3[N], k4[N], yt[N];

    dydx(y, k1, phi, Rows, Columns, delta);

    for (int i = 0; i < N; i++)
        yt[i] = y[i] + 0.5 * h * k1[i];

    dydx(yt, k2, phi, Rows, Columns, delta);

    for (int i = 0; i < N; i++)
        yt[i] = y[i] + 0.5 * h * k2[i];

    dydx(yt, k3, phi, Rows, Columns, delta);

    for (int i = 0; i < N; i++)
        yt[i] = y[i] + h * k3[i];

    dydx(yt, k4, phi, Rows, Columns, delta);

    for (int i = 0; i < N; i++) {
        y[i] +=
            (h / 6.0)
            * (k1[i]
            + 2.0 * k2[i]
            + 2.0 * k3[i]
            + k4[i]);
    }
}

// ======================================================
// MAIN
// ======================================================

int main() {

    const int Rows = 301;
    const int Columns = 301;

    double delta = 1.0;

    int LICZBA_GWIAZD = 300000;

    double MASA_SMBH = 40.0;
    double MASA_GWIAZDY = 0.02;

    int steps = 1300;

    double h = 0.0005;

    vector<double> phi(Rows * Columns, 0.0);
    vector<double> rho(Rows * Columns, 0.0);

    vector<vector<double>> gwiazdy;

    cout << "Initializing galaxy...\n";

    inicjalizuj_galaktyke(
        gwiazdy,
        LICZBA_GWIAZD,
        Rows,
        Columns,
        MASA_SMBH);

    //initialize_edges(phi, Rows, Columns);

    ofstream fout("orbity_gwiazd.txt");

    fout << fixed << setprecision(3);

    cout << "Running simulation...\n";

    for (int step = 0; step < steps; step++) {

        // 1. deposit mass
        aktualizuj_gestosc(
            rho,
            gwiazdy,
            Rows,
            Columns,
            MASA_SMBH,
            MASA_GWIAZDY);


        wygladz_gestosc(
            rho,
            Rows,
            Columns
        );

        aktualizuj_brzegi(
                phi,
                rho,
                Rows,
                Columns,
                delta);

        // 2. solve Poisson
        for (int i = 0; i < 100; i++) {

            nadrelaksacja(
                phi,
                rho,
                Rows,
                Columns,
                delta);
        }

        // 3. move particles
        for (auto& g : gwiazdy) {

            rungeKutta(
                h,
                g.data(),
                phi,
                Rows,
                Columns,
                delta);
        }

        // save
        if (step % 5 == 0) {

            for (const auto& g : gwiazdy) {

                fout << g[0] << "\t"
                     << g[2] << "\t";
            }

            fout << "\n";
        }

        if (step % 10 == 0) {

            cout
                << "Step "
                << step
                << " / "
                << steps
                << "\r"
                << flush;
        }
    }

    cout << "\nDone.\n";

    return 0;
}