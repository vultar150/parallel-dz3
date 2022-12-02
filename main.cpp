#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <cmath>
#include <mpi.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#define TAG_X 1215
#define TAG_Y 1216
#define TAG_Z 1217


#define SIZE 128

uint32_t SIZE_X = SIZE;
uint32_t SIZE_Y = SIZE;
uint32_t SIZE_Z = SIZE;
const double T = 1.0;

/* tau = T / MAX_IT must be <= h / 3^(1/3)
 * where h = L / (SIZE - 1) */
const uint32_t MAX_IT = 1000;

double Lx = 1.0, Ly = 1.0, Lz = 1.0;
double hx, hy, hz;
double a_t;

uint32_t start_x = 0, start_y = 0, start_z = 0;


double u_analytical(double x, double y, double z, double t) {
    return sin(x * 2. * M_PI / Lx) * 
           sin(y * M_PI / Ly) *
           sin(z * M_PI / Lz) *
           cos(a_t * t + 2. * M_PI);
}


double phi(double x, double y, double z) {
    return u_analytical(x, y, z, 0.0);
}


struct Matrix2 {
    uint32_t Nx, Ny;
    double *data;

    Matrix2(uint32_t Nx, uint32_t Ny): Nx(Nx), Ny(Ny) {
        data = new double[Nx * Ny]();
    }

    double& Get(uint32_t i, uint32_t j) {
        return data[i*Ny + j];
    }

    void PrintAll() {
        for (uint32_t i = 0; i < Nx; ++i) {
            for (uint32_t j = 0; j < Ny; ++j) {
                    printf("[%d, %d]: %.10lf\n", i, j, Get(i, j));
            }
        }
    }

    ~Matrix2() {
        delete []data; 
    }
};


struct Matrix3 {
    uint32_t Nx, Ny, Nz;
    uint32_t size;
    double *data;

    Matrix3(uint32_t Nx, uint32_t Ny, uint32_t Nz):
                    Nx(Nx), Ny(Ny), Nz(Nz), size(Nx * Ny * Nz) {
        data = new double[size]();
    }

    double& Get(uint32_t i, uint32_t j, uint32_t k) {
        return data[i*Ny*Nz + j*Nz + k];
    }

    void PrintAll(int num_procs[], int coords[]) {
        printf("Block: [%d, %d, %d]\n", coords[0], coords[1], coords[2]);
        uint32_t LNx = coords[0] == num_procs[0] - 1 ? Nx : Nx - 1;
        uint32_t LNy = coords[1] == num_procs[1] - 1 ? Ny : Ny - 1;
        uint32_t LNz = coords[2] == num_procs[2] - 1 ? Nz : Nz - 1;
        for (uint32_t i = start_x == 0 ? 0 : 1; i < LNx; ++i) {
            for (uint32_t j = start_y == 0 ? 0 : 1; j < LNy; ++j) {
                for (uint32_t k = start_z == 0 ? 0 : 1; k < LNz; ++k) {
                    printf("[%d, %d, %d] = %.10lf\n", 
                            start_x + i, start_y + j, start_z + k, Get(i, j, k));
                }
            }
        }
    }

    ~Matrix3() {
        delete []data; 
    }
};


class Process {
public:
    Process(uint32_t Nx, uint32_t Ny, uint32_t Nz, int rank, int num_procs[], int coords[], 
            double T, uint32_t max_it, MPI_Comm GRID_COMM):
                curr_block(Nx, Ny, Nz), 
                p_block(Nx, Ny, Nz), 
                pp_block(Nx, Ny, Nz),
                left_recv(Nx, Nz), left_send(Nx, Nz),
                right_recv(Nx, Nz), right_send(Nx, Nz),
                up_recv(Nx, Ny), up_send(Nx, Ny),
                down_recv(Nx, Ny), down_send(Nx, Ny),
                rank(rank),
                T(T), max_it(max_it), tau(T / (double) max_it),
                n(2), p_n(1), pp_n(0),
                GRID_COMM(GRID_COMM),
                Nx(Nx), Ny(Ny), Nz(Nz) {
        for (uint32_t i = 0; i < 3; ++i) {
            this->num_procs[i] = num_procs[i];
            this->coords[i] = coords[i];
        }
        blocks[0] = &pp_block;
        blocks[1] = &p_block;
        blocks[2] = &curr_block;
    }

    void Solve() {
        uint32_t it = 0;
        double local_eps = 0.0;
        double global_eps = 0.0;
        double t = 0.0;
        Step0();
        local_eps = GetEps(pp_n, t);
        MPI_Reduce(&local_eps, &global_eps, 1, MPI_DOUBLE, 
                   MPI_MAX, 0, GRID_COMM);
        if (rank == 0) { printf("t = %lf, it = %d, eps = %.10lf\n", t, it, global_eps); }
        it++;
        
        Step1();
        ExchangePlanes_X(p_n);
        t = tau * it;
        local_eps = GetEps(p_n, t);
        MPI_Reduce(&local_eps, &global_eps, 1, MPI_DOUBLE, 
                   MPI_MAX, 0, GRID_COMM);
        if (rank == 0) { printf("t = %lf, it = %d, eps = %.10lf\n", t, it, global_eps); }
        it++;

        while (it <= max_it) {
            PreparePlanes_Y_Z_ToSend(p_n);
            ExchangePlanes_Y_Z();
            CopyPlanes_Y_Z_AfterRecieve(p_n);
            StepN();
            ExchangePlanes_X(n);
            t = tau * it;
            local_eps = GetEps(n, t);
            MPI_Reduce(&local_eps, &global_eps, 1, MPI_DOUBLE, 
                        MPI_MAX, 0, GRID_COMM);
            if (rank == 0) { printf("t = %lf, it = %d, eps = %.10lf\n", t, it, global_eps); }
            it++;
            UpdateIterationNumber();
        }
    }


private:
    int rank;
    uint32_t Nx, Ny, Nz;
    double T, tau;
    double eps;
    uint32_t max_it;
    MPI_Comm GRID_COMM;
    int num_procs[3];
    int coords[3];
    uint32_t n, p_n, pp_n;
    Matrix3 *blocks[3];
    Matrix3 curr_block, p_block, pp_block; // current, prev, pprev
    Matrix2 left_recv, right_recv, up_recv, down_recv;
    Matrix2 left_send, right_send, up_send, down_send;

    void UpdateIterationNumber() {
        n = (n + 1) % 3;
        p_n = (p_n + 1) % 3;
        pp_n = (pp_n + 1) % 3;
    }

    double delta_h(uint32_t ln, uint32_t i, uint32_t j, uint32_t k) {
        return ( ( blocks[ln]->Get(i - 1, j, k) - 
                   2 * blocks[ln]->Get(i, j, k) + 
                   blocks[ln]->Get(i + 1, j, k) ) / 
                 ( hx * hx ) ) +
               ( ( blocks[ln]->Get(i, j - 1, k) - 
                   2 * blocks[ln]->Get(i, j, k) + 
                   blocks[ln]->Get(i, j + 1, k) ) / 
                 ( hy * hy ) ) +
               ( ( blocks[ln]->Get(i, j, k - 1) - 
                   2 * blocks[ln]->Get(i, j, k) + 
                   blocks[ln]->Get(i, j, k + 1) ) / 
                 ( hz * hz ) );
    }

    void Step0() {
        for (uint32_t i = 0; i < Nx; ++i) {
            for (uint32_t j = 0; j < Ny; ++j) {
                for (uint32_t k = 0; k < Nz; ++k) {
                    blocks[pp_n]->Get(i, j, k) = 
                                    phi(hx * (start_x + i), 
                                        hy * (start_y + j), 
                                        hz * (start_z + k));
                }
            }
        }
    }

    void Step1() {
        uint32_t LNy = coords[1] < num_procs[1] - 1 ? Ny : Ny - 1;
        uint32_t LNz = coords[2] < num_procs[2] - 1 ? Nz : Nz - 1;
        for (uint32_t i = 1; i < Nx - 1; ++i) {
            for (uint32_t j = 1; j < LNy - 1; ++j) {
                for (uint32_t k = 1; k < LNz - 1; ++k) {
                    blocks[p_n]->Get(i, j, k) = 
                                    blocks[pp_n]->Get(i, j, k) + 
                                    (tau * tau / 2.) * 
                                    delta_h(pp_n, i, j, k);
                }
            }
        }
    }

    void StepN() {
        uint32_t LNy = coords[1] < num_procs[1] - 1 ? Ny : Ny - 1;
        uint32_t LNz = coords[2] < num_procs[2] - 1 ? Nz : Nz - 1;
        for (uint32_t i = 1; i < Nx - 1; ++i) {
            for (uint32_t j = 1; j < LNy - 1; ++j) {
                for (uint32_t k = 1; k < LNz - 1; ++k) {
                    blocks[n]->Get(i, j, k) = 
                                    2. * blocks[p_n]->Get(i, j, k) - 
                                    blocks[pp_n]->Get(i, j, k) + 
                                    tau * tau * delta_h(p_n, i, j, k);
                }
            }
        }
    }

    void PreparePlanes_Y_Z_ToSend(uint32_t ln) {
        PreparePlanes_Y_ToSend(ln);
        PreparePlanes_Z_ToSend(ln);
    }

    void PreparePlanes_Y_ToSend(uint32_t ln) {
        if (coords[1] > 0) {
            for (uint32_t i = 0; i < Nx; ++i) {
                for (uint32_t k = 0; k < Nz; ++k) {
                    left_send.Get(i, k) = blocks[ln]->Get(i, 1, k);
                }
            }
        }
        if (coords[1] < num_procs[1] - 1) {
            for (uint32_t i = 0; i < Nx; ++i) {
                for (uint32_t k = 0; k < Nz; ++k) {
                    right_send.Get(i, k) = blocks[ln]->Get(i, Ny - 2, k);
                }
            }    
        }
    }

    void PreparePlanes_Z_ToSend(uint32_t ln) {
        if (coords[2] > 0) {
            for (uint32_t i = 0; i < Nx; ++i) {
                for (uint32_t j = 0; j < Ny; ++j) {
                    down_send.Get(i, j) = blocks[ln]->Get(i, j, 1);
                }
            }
        }

        if (coords[2] < num_procs[2] - 1) {
            for (uint32_t i = 0; i < Nx; ++i) {
                for (uint32_t j = 0; j < Ny; ++j) {
                    up_send.Get(i, j) = blocks[ln]->Get(i, j, Nz - 2);
                }
            }
        }
    }

    void CopyPlanes_Y_Z_AfterRecieve(uint32_t ln) { 
        CopyPlanes_Y_AfterRecieve(ln);
        CopyPlanes_Z_AfterRecieve(ln);
    }

    void CopyPlanes_Y_AfterRecieve(uint32_t ln) {
        if (coords[1] > 0) {
            for (uint32_t i = 0; i < Nx; ++i) {
                for (uint32_t k = 0; k < Nz; ++k) {
                    blocks[ln]->Get(i, 0, k) = left_recv.Get(i, k);
                }
            }
        }
        if (coords[1] < num_procs[1] - 1) {
            for (uint32_t i = 0; i < Nx; ++i) {
                for (uint32_t k = 0; k < Nz; ++k) {
                    blocks[ln]->Get(i, Ny - 1, k) = right_recv.Get(i, k);
                }
            }    
        }
    }

    void CopyPlanes_Z_AfterRecieve(uint32_t ln) {
        if (coords[2] > 0) {
            for (uint32_t i = 0; i < Nx; ++i) {
                for (uint32_t j = 0; j < Ny; ++j) {
                    blocks[ln]->Get(i, j, 0) = down_recv.Get(i, j);
                }
            }
        }

        if (coords[2] < num_procs[2] - 1) {
            for (uint32_t i = 0; i < Nx; ++i) {
                for (uint32_t j = 0; j < Ny; ++j) {
                    blocks[ln]->Get(i, j, Nz - 1) = up_recv.Get(i, j);
                }
            }
        }
    }

    void ExchangePlanes_Y_Z() {
        ExchangePlanes_Y();
        ExchangePlanes_Z();
    }

    void ExchangePlanes_X(uint32_t ln) {
        int rank_recv, rank_send;
        MPI_Status status;
        MPI_Cart_shift(GRID_COMM, 0, -1, &rank_recv, &rank_send);
        MPI_Sendrecv(blocks[ln]->data + 1 * Ny * Nz, Ny * Nz, 
                               MPI_DOUBLE, rank_send, TAG_X, 
                     blocks[ln]->data + (Nx - 1) * Ny * Nz, Ny * Nz,
                               MPI_DOUBLE, rank_recv, TAG_X, 
                     GRID_COMM, &status);

        MPI_Cart_shift(GRID_COMM, 0, 1, &rank_recv, &rank_send);
        MPI_Sendrecv(blocks[ln]->data + (Nx - 2) * Ny * Nz, Ny * Nz, 
                               MPI_DOUBLE, rank_send, TAG_X, 
                     blocks[ln]->data, Ny * Nz,
                               MPI_DOUBLE, rank_recv, TAG_X, 
                     GRID_COMM, &status);
    }

    void ExchangePlanes_Y() {
        int rank_recv, rank_send;
        MPI_Status status;
        MPI_Cart_shift(GRID_COMM, 1, -1, &rank_recv, &rank_send);
        if (coords[1] != 0 && coords[1] != num_procs[1] - 1) {
            MPI_Sendrecv(left_send.data, Nx * Nz, MPI_DOUBLE, rank_send, TAG_Y, 
                        right_recv.data, Nx * Nz, MPI_DOUBLE, rank_recv, TAG_Y, 
                        GRID_COMM, &status);
        } else if (coords[1] == 0 && coords[1] != num_procs[1] - 1) {
            MPI_Recv(right_recv.data, Nx * Nz, MPI_DOUBLE, rank_recv, TAG_Y, 
                    GRID_COMM, &status);
        } else if (coords[1] == num_procs[1] - 1 && coords[1] != 0){
            MPI_Send(left_send.data, Nx * Nz, MPI_DOUBLE, rank_send, TAG_Y, 
                    GRID_COMM);
        }

        MPI_Cart_shift(GRID_COMM, 1, 1, &rank_recv, &rank_send);
        if (coords[1] != 0 && coords[1] != num_procs[1] - 1) {
            MPI_Sendrecv(right_send.data, Nx * Nz, MPI_DOUBLE, rank_send, TAG_Y, 
                        left_recv.data, Nx * Nz, MPI_DOUBLE, rank_recv, TAG_Y, 
                        GRID_COMM, &status);
        } else if (coords[1] == 0 && coords[1] != num_procs[1] - 1){
            MPI_Send(right_send.data, Nx * Nz, MPI_DOUBLE, rank_send, TAG_Y, 
                    GRID_COMM);
        } else if (coords[1] == num_procs[1] - 1 && coords[1] != 0){
            MPI_Recv(left_recv.data, Nx * Nz, MPI_DOUBLE, rank_recv, TAG_Y, 
                    GRID_COMM, &status);
        } 
    }

    void ExchangePlanes_Z() {
        int rank_recv, rank_send;
        MPI_Status status;
        MPI_Cart_shift(GRID_COMM, 2, -1, &rank_recv, &rank_send);
        if (coords[2] != 0 && coords[2] != num_procs[2] - 1) {
            MPI_Sendrecv(down_send.data, Nx * Ny, MPI_DOUBLE, rank_send, TAG_Z, 
                        up_recv.data, Nx * Ny, MPI_DOUBLE, rank_recv, TAG_Z, 
                        GRID_COMM, &status);
        } else if (coords[2] == 0 && coords[2] != num_procs[2] - 1){
            MPI_Recv(up_recv.data, Nx * Ny, MPI_DOUBLE, rank_recv, TAG_Z, 
                    GRID_COMM, &status);
        } else if (coords[2] == num_procs[2] - 1 && coords[2] != 0){
            MPI_Send(down_send.data, Nx * Ny, MPI_DOUBLE, rank_send, TAG_Z, 
                    GRID_COMM);
        }

        MPI_Cart_shift(GRID_COMM, 2, 1, &rank_recv, &rank_send);
        if (coords[2] != 0 && coords[2] != num_procs[2] - 1) {
            MPI_Sendrecv(up_send.data, Nx * Ny, MPI_DOUBLE, rank_send, TAG_Z, 
                        down_recv.data, Nx * Ny, MPI_DOUBLE, rank_recv, TAG_Z, 
                        GRID_COMM, &status);
        } else if (coords[2] == 0 && coords[2] != num_procs[2] - 1){
            MPI_Send(up_send.data, Nx * Ny, MPI_DOUBLE, rank_send, TAG_Z, 
                    GRID_COMM);
        } else if (coords[2] == num_procs[2] - 1 && coords[2] != 0){
            MPI_Recv(down_recv.data, Nx * Ny, MPI_DOUBLE, rank_recv, TAG_Z, 
                    GRID_COMM, &status);
        } 
    }

    double GetEps(uint32_t ln, double t) {
        double max_eps = 0.0;
        double tmp = 0.0;

        for (uint32_t i = coords[0] == 0 ? 0 : 1; i < Nx - 1; i++) {
            for (uint32_t j = coords[1] == 0 ? 0 : 1; j < Ny - 1; j++) {
                for (uint32_t k = coords[2] == 0 ? 0 : 1; k < Nz - 1; k++) {
                    tmp = fabs( blocks[ln]->Get(i, j, k) - 
                                u_analytical(hx * (start_x + i), 
                                             hy * (start_y + j), 
                                             hz * (start_z + k), t) );
                    if (max_eps < tmp) {
                        max_eps = tmp;
                    }
                }
            }
        }
        return max_eps;
    }
};


void ParseAndInit(int argc, char **argv, const int num_procs[],
                    const int coords[], uint32_t &Nx, uint32_t &Ny, uint32_t &Nz) {
    if (argc >= 7) {
        sscanf(argv[1], "%lf", &Lx);
        sscanf(argv[2], "%lf", &Ly);
        sscanf(argv[3], "%lf", &Lz);
	sscanf(argv[4], "%d", &SIZE_X);
	sscanf(argv[5], "%d", &SIZE_Y);
	sscanf(argv[6], "%d", &SIZE_Z);
    }
    a_t = M_PI * sqrt(4.0 / (Lx * Lx) + 1.0 / (Ly * Ly) + 1.0 / (Lz * Lz));
    hx = Lx / (double) (SIZE_X - 1);
    hy = Ly / (double) (SIZE_Y - 1);
    hz = Lz / (double) (SIZE_Z - 1);

    uint32_t Nmin_x = (SIZE_X / num_procs[0]);
    uint32_t Nmin_y = (SIZE_Y / num_procs[1]);
    uint32_t Nmin_z = (SIZE_Z / num_procs[2]);

    uint32_t Nextra_x = SIZE_X % num_procs[0];
    uint32_t Nextra_y = SIZE_Y % num_procs[1];
    uint32_t Nextra_z = SIZE_Z % num_procs[2];

    auto getCoordSize = [](const int &coord, const uint32_t &Nextra, const uint32_t &Nmin) {
        if (coord != 0) {
            return coord < Nextra ? Nmin + 3 : Nmin + 2;
        } else {
            return coord < Nextra ? Nmin + 2 : Nmin + 1;
        }
    };

    Nx = getCoordSize(coords[0], Nextra_x, Nmin_x);
    Ny = getCoordSize(coords[1], Nextra_y, Nmin_y);
    Nz = getCoordSize(coords[2], Nextra_z, Nmin_z);

    auto getStartN = [](const int &coord, const uint32_t &Nextra, const uint32_t &Nmin) {
        uint32_t tmp = 0;
        uint32_t start = 0;
        for (uint32_t i = 0; i < coord; ++i) {
            tmp = i < Nextra ? Nmin + 1 : Nmin;
            start += tmp;
        }
        return start > 0 ? start - 1 : start;
    };
    
    start_x = getStartN(coords[0], Nextra_x, Nmin_x);
    start_y = getStartN(coords[1], Nextra_y, Nmin_y);
    start_z = getStartN(coords[2], Nextra_z, Nmin_z);
}


int main(int argc, char **argv) {
    uint32_t Nx, Ny, Nz;
    int num_procs[3] = {0, 0, 0};
    int periods[3] = {1, 0, 0};
    int coords[3];
    int rank_size;
    int rank;

    MPI_Comm GRID_COMM;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &rank_size);
    MPI_Dims_create(rank_size, 3, num_procs);
    MPI_Cart_create(MPI_COMM_WORLD, 3, num_procs, periods, 1, &GRID_COMM);
    MPI_Comm_rank(GRID_COMM, &rank);
    MPI_Cart_coords(GRID_COMM, rank, 3, coords);

    ParseAndInit(argc, argv, num_procs, coords, Nx, Ny, Nz);
    double start_time = MPI_Wtime();
    Process P(Nx, Ny, Nz, rank, num_procs, coords, T, MAX_IT, GRID_COMM);
    P.Solve();
    double time = MPI_Wtime() - start_time, global_time;
    MPI_Reduce(&time, &global_time, 1, MPI_DOUBLE, 
		MPI_MAX, 0, GRID_COMM);
    if (rank == 0) {
	printf("Num processes = [%d, %d, %d]\n", num_procs[0], num_procs[1], num_procs[2]);
	printf("time = %lf\n", global_time);
    }
    

    MPI_Finalize();
    return 0;
}
