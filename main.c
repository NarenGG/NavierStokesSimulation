#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <opencv2/opencv.hpp>

#define DEFINITION 50
#define NUM_TIME_STEPS 1000

typedef struct {
    double u, v, w;
} Velocity;

void initialize_pressure_grid(double pressure_grid[DEFINITION+2][DEFINITION+2][DEFINITION+2]) {
    for (int i = 0; i < DEFINITION + 2; i++) {
        for (int j = 0; j < DEFINITION + 2; j++) {
            for (int k = 0; k < (DEFINITION / 2) + 1; k++) {
                pressure_grid[i][j][k] = 1.0;
            }
        }
    }
}

void initialize_velocity_grid(Velocity velocity_grid[DEFINITION+2][DEFINITION+2][DEFINITION+2]) {
    for (int i = 0; i < DEFINITION + 2; i++) {
        for (int j = 0; j < DEFINITION + 2; j++) {
            for (int k = 0; k < DEFINITION + 2; k++) {
                velocity_grid[i][j][k].u = 0.0;
                velocity_grid[i][j][k].v = 0.0;
                velocity_grid[i][j][k].w = 0.0;
            }
        }
    }
}

void calculate_next_state(double pressure_grid[DEFINITION+2][DEFINITION+2][DEFINITION+2], Velocity velocity_grid[DEFINITION+2][DEFINITION+2][DEFINITION+2], double rho, double mu, double dt, double dx, double dy, double dz, double next_pressure_grid[DEFINITION+2][DEFINITION+2][DEFINITION+2], Velocity next_velocity_grid[DEFINITION+2][DEFINITION+2][DEFINITION+2]) {
    for (int i = 1; i < DEFINITION + 1; i++) {
        for (int j = 1; j < DEFINITION + 1; j++) {
            for (int k = 1; k < DEFINITION + 1; k++) {
                double du_dx = (velocity_grid[i+1][j][k].u - velocity_grid[i-1][j][k].u) / (2 * dx);
                double dv_dy = (velocity_grid[i][j+1][k].v - velocity_grid[i][j-1][k].v) / (2 * dy);
                double dw_dz = (velocity_grid[i][j][k+1].w - velocity_grid[i][j][k-1].w) / (2 * dz);

                double laplacian_p = (
                    (pressure_grid[i+1][j][k] - 2 * pressure_grid[i][j][k] + pressure_grid[i-1][j][k]) / pow(dx, 2) +
                    (pressure_grid[i][j+1][k] - 2 * pressure_grid[i][j][k] + pressure_grid[i][j-1][k]) / pow(dy, 2) +
                    (pressure_grid[i][j][k+1] - 2 * pressure_grid[i][j][k] + pressure_grid[i][j][k-1]) / pow(dz, 2)
                );

                next_pressure_grid[i][j][k] = pressure_grid[i][j][k] + dt * (
                    -rho * (velocity_grid[i][j][k].u * du_dx + velocity_grid[i][j][k].v * dv_dy + velocity_grid[i][j][k].w * dw_dz) +
                    mu * laplacian_p
                );

                next_velocity_grid[i][j][k].u = fmin(fmax(velocity_grid[i][j][k].u + dt * (
                    -velocity_grid[i][j][k].u * du_dx - velocity_grid[i][j][k].v * dv_dy - velocity_grid[i][j][k].w * dw_dz +
                    mu * (laplacian_p / rho)
                ), -DBL_MAX), DBL_MAX);

                next_velocity_grid[i][j][k].v = fmin(fmax(velocity_grid[i][j][k].v + dt * (
                    -velocity_grid[i][j][k].u * du_dx - velocity_grid[i][j][k].v * dv_dy - velocity_grid[i][j][k].w * dw_dz +
                    mu * (laplacian_p / rho)
                ), -DBL_MAX), DBL_MAX);

                next_velocity_grid[i][j][k].w = fmin(fmax(velocity_grid[i][j][k].w + dt * (
                    -velocity_grid[i][j][k].u * du_dx - velocity_grid[i][j][k].v * dv_dy - velocity_grid[i][j][k].w * dw_dz +
                    mu * (laplacian_p / rho)
                ), -DBL_MAX), DBL_MAX);
            }
        }
    }
}

void plot_slices(double grid[DEFINITION+2][DEFINITION+2][DEFINITION+2], const char *filename) {
    cv::Mat slice_xy(DEFINITION, DEFINITION, CV_64F, &grid[DEFINITION/2][1][1]);
    cv::Mat slice_xz(DEFINITION, DEFINITION, CV_64F, &grid[1][DEFINITION/2][1]);
    cv::Mat slice_yz(DEFINITION, DEFINITION, CV_64F, &grid[1][1][DEFINITION/2]);

    cv::imwrite(filename, slice_xy);
}

int main() {
    double pressure_grid[DEFINITION+2][DEFINITION+2][DEFINITION+2] = {0};
    Velocity velocity_grid[DEFINITION+2][DEFINITION+2][DEFINITION+2] = {0};
    double rho = 1.0, mu = 0.01, dt = 0.1, dx = 1.0, dy = 1.0, dz = 1.0;

    initialize_pressure_grid(pressure_grid);
    initialize_velocity_grid(velocity_grid);

    double next_pressure_grid[DEFINITION+2][DEFINITION+2][DEFINITION+2];
    Velocity next_velocity_grid[DEFINITION+2][DEFINITION+2][DEFINITION+2];

    for (int epoch = 1; epoch <= NUM_TIME_STEPS; epoch++) {
        calculate_next_state(pressure_grid, velocity_grid, rho, mu, dt, dx, dy, dz, next_pressure_grid, next_velocity_grid);

        memcpy(pressure_grid, next_pressure_grid, sizeof(double) * (DEFINITION+2) * (DEFINITION+2) * (DEFINITION+2));
        memcpy(velocity_grid, next_velocity_grid, sizeof(Velocity) * (DEFINITION+2) * (DEFINITION+2) * (DEFINITION+2));

        char filename[256];
        sprintf(filename, "data/%d.png", epoch);
        plot_slices(pressure_grid, filename);

        printf("Processed Grid at Epoch %d\n", epoch);
    }

    return 0;
}