#include <stdio.h>
#include <stdint.h>
#include <math.h>

#define SIZE 128

#define NEXT(n_, p_n_, pp_n_) \
                n_ = (n_ + 1) % 3; \
                p_n_ = (p_n_ + 1) % 3; \
                pp_n_ = (pp_n_ + 1) % 3 \

const uint32_t N = SIZE - 1;
uint32_t n = 2, p_n = 1, pp_n = 0;


double u[3][SIZE + 1][SIZE + 1][SIZE + 1];

double Lx = 1., Ly = 1., Lz = 1.;
const uint32_t MAX_IT = 220;
const double T = 1.;
const double tau = T / (double) MAX_IT;
double h = 1. / (double) N;
double a_t = 0.;

void compute();
void init();
double u_analytical(double x, double y, double z, double t);
double phi(double x, double y, double z);
void parse(int argc, char **argv);
double delta_h(uint32_t n, uint32_t i, uint32_t j, uint32_t k);
void print_eps(uint32_t it, double t);


int main(int argc, char **argv)
{

    parse(argc, argv);
    init();
    uint32_t it = 2;
    while (it <= MAX_IT)
    {   
        compute();
        print_eps(n, tau * it);
        NEXT(n, p_n, pp_n);
        it++;
    }
    return 0;
}


void compute()
{
    for (uint32_t i = 1; i <= N; i++) 
    {
        for (uint32_t j = 1; j <= N - 1; j++) 
        {
            for (uint32_t k = 1; k <= N - 1; k++)
            {
                u[n][i][j][k] = 2. * u[p_n][i][j][k] - u[pp_n][i][j][k] + tau * tau * delta_h(p_n, i, j, k);
                if (i == 1)
                {
                    u[n][N + 1][j][k] = u[n][i][j][k];
                }
            }
        }
    }

    for (uint32_t j = 1; j <= N - 1; j++)
    {
        for (uint32_t k = 1; k <= N - 1; k++)
        {
            u[n][0][j][k] = u[n][N][j][k];
        }
    }
}


void init()
{
    for (uint32_t i = 0; i <= N + 1; i++) 
    {
        for (uint32_t j = 0; j <= N + 1; j++) 
        {
            for (uint32_t k = 0; k <= N + 1; k++)
            {
                u[pp_n][i][j][k] = phi(h * i, h * j, h * k);
                if (k == 0 || k == N || j == 0 || j == N)
                {
                    u[n][i][j][k] = u[p_n][i][j][k] = u[pp_n][i][j][k] = 0.;
                }
            }
        }
    }

    print_eps(pp_n, 0.);


    for (uint32_t i = 1; i <= N; i++) 
    {
        for (uint32_t j = 1; j <= N - 1; j++) 
        {
            for (uint32_t k = 1; k <= N - 1; k++)
            {
                u[p_n][i][j][k] = u[pp_n][i][j][k] + (tau * tau / 2.) * delta_h(pp_n, i, j, k);
                if (i == 1) 
                {
                    u[p_n][N + 1][j][k] = u[p_n][i][j][k];
                }
            }
        }
    }

    for (uint32_t j = 1; j <= N - 1; j++)
    {
        for (uint32_t k = 1; k <= N - 1; k++)
        {
            u[p_n][0][j][k] = u[p_n][N][j][k];
        }
    }

    print_eps(p_n, tau);
}


double u_analytical(double x, double y, double z, double t)
{
    return sin(x * 2. * M_PI / Lx) * 
           sin(y * M_PI / Ly) *
           sin(z * M_PI / Lz) *
           cos(a_t * t + 2. * M_PI);
}


double phi(double x, double y, double z)
{
    return sin(x * 2. * M_PI / Lx) * 
           sin(y * M_PI / Ly) *
           sin(z * M_PI / Lz);
}


void parse(int argc, char **argv)
{
    if (argc >= 5) 
    {
        sscanf(argv[1], "%lf", &Lx);
        sscanf(argv[2], "%lf", &Ly);
        sscanf(argv[3], "%lf", &Lz);
        // sscanf(argv[4], "%lf", &tau);
    }
    a_t = M_PI * sqrt(4.0 / (Lx * Lx) + 1.0 / (Ly * Ly) + 1.0 / (Lz * Lz));
    h = Lx / (double) N;
}


double delta_h(uint32_t n, uint32_t i, uint32_t j, uint32_t k)
{
    return ((u[n][i - 1][j][k] - 2 * u[n][i][j][k] + u[n][i + 1][j][k]) / (h * h)) +
           ((u[n][i][j - 1][k] - 2 * u[n][i][j][k] + u[n][i][j + 1][k]) / (h * h)) +
           ((u[n][i][j][k - 1] - 2 * u[n][i][j][k] + u[n][i][j][k + 1]) / (h * h));
}


void print_eps(uint32_t it, double t)
{
    double max_eps = 0.;
    for (uint32_t i = 0; i <= N; i++) 
    {
        for (uint32_t j = 0; j <= N; j++) 
        {
            for (uint32_t k = 0; k <= N; k++)
            {
                if (max_eps < fabs(u[it][i][j][k] - u_analytical(h * i, h * j, h * k, t))) 
                {
                    max_eps = fabs(u[it][i][j][k] - u_analytical(h * i, h * j, h * k, t));
                }
            }
        }
    }
    printf("t = %lf, eps = %lf\n", t, max_eps);
}
