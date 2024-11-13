#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct {
    int size;
    float dt, diff, visc, kappa;
    float *s, *density, *Vx, *Vy, *Vx0, *Vy0;
} FluidCube;

static int idx(FluidCube *cube, int x, int y) {
    return x + y * cube->size;
}

FluidCube* FluidCube_create(int size, float diffusion, float viscosity, float dt, float kappa) {
    FluidCube *cube = (FluidCube*)malloc(sizeof(FluidCube));
    cube->size = size;
    cube->dt = dt;
    cube->diff = diffusion;
    cube->visc = viscosity;
    cube->kappa = kappa;

    int n = size * size;
    cube->s = (float*)malloc(n * sizeof(float));
    cube->density = (float*)malloc(n * sizeof(float));
    cube->Vx = (float*)malloc(n * sizeof(float));
    cube->Vy = (float*)malloc(n * sizeof(float));
    cube->Vx0 = (float*)malloc(n * sizeof(float));
    cube->Vy0 = (float*)malloc(n * sizeof(float));

    for (int i = 0; i < n; i++) {
        cube->s[i] = cube->density[i] = cube->Vx[i] = cube->Vy[i] = cube->Vx0[i] = cube->Vy0[i] = 0.0f;
    }

    return cube;
}

void FluidCube_addDensity(FluidCube *cube, int x, int y, float amount) {
    cube->density[idx(cube, x, y)] += amount;
}

void FluidCube_addVelocity(FluidCube *cube, int x, int y, float amountX, float amountY) {
    int idx_ = idx(cube, x, y);
    cube->Vx[idx_] += amountX;
    cube->Vy[idx_] += amountY;
}

static void setBound(FluidCube *cube, int b, float *x) {
    int size = cube->size;
    for (int i = 1; i < size - 1; i++) {
        x[idx(cube, i, 0)] = b == 2 ? -x[idx(cube, i, 1)] : x[idx(cube, i, 1)];
        x[idx(cube, i, size - 1)] = b == 2 ? -x[idx(cube, i, size - 2)] : x[idx(cube, i, size - 2)];
    }
    for (int j = 1; j < size - 1; j++) {
        x[idx(cube, 0, j)] = b == 1 ? -x[idx(cube, 1, j)] : x[idx(cube, 1, j)];
        x[idx(cube, size - 1, j)] = b == 1 ? -x[idx(cube, size - 2, j)] : x[idx(cube, size - 2, j)];
    }

    x[idx(cube, 0, 0)] = 0.5f * (x[idx(cube, 1, 0)] + x[idx(cube, 0, 1)]);
    x[idx(cube, 0, size - 1)] = 0.5f * (x[idx(cube, 1, size - 1)] + x[idx(cube, 0, size - 2)]);
    x[idx(cube, size - 1, 0)] = 0.5f * (x[idx(cube, size - 2, 0)] + x[idx(cube, size - 1, 1)]);
    x[idx(cube, size - 1, size - 1)] = 0.5f * (x[idx(cube, size - 2, size - 1)] + x[idx(cube, size - 1, size - 2)]);
}

static void linSolve(FluidCube *cube, int b, float *x, float *x0, float a, float c, int iter) {
    int size = cube->size;
    float cRecip = 1.0f / c;
    for (int k = 0; k < iter; k++) {
        for (int j = 1; j < size - 1; j++) {
            for (int i = 1; i < size - 1; i++) {
                x[idx(cube, i, j)] = (x0[idx(cube, i, j)] + a * (x[idx(cube, i + 1, j)] + x[idx(cube, i - 1, j)] + x[idx(cube, i, j + 1)] + x[idx(cube, i, j - 1)])) * cRecip;
            }
        }
        setBound(cube, b, x);
    }
}

void FluidCube_diffuse(FluidCube *cube, int b, float *x, float *x0, float diff, float dt, int iter) {
    float a = dt * diff * (cube->size - 2) * (cube->size - 2);
    linSolve(cube, b, x, x0, a, 1 + 4 * a, iter);
}

void FluidCube_project(FluidCube *cube, float *velocX, float *velocY, float *p, float *div, int iter) {
    int size = cube->size;
    for (int j = 1; j < size - 1; j++) {
        for (int i = 1; i < size - 1; i++) {
            div[idx(cube, i, j)] = -0.5f * (velocX[idx(cube, i + 1, j)] - velocX[idx(cube, i - 1, j)] + velocY[idx(cube, i, j + 1)] - velocY[idx(cube, i, j - 1)]) / size;
            p[idx(cube, i, j)] = 0;
        }
    }
    setBound(cube, 0, div);
    setBound(cube, 0, p);
    linSolve(cube, 0, p, div, 1, 4, iter);

    for (int j = 1; j < size - 1; j++) {
        for (int i = 1; i < size - 1; i++) {
            velocX[idx(cube, i, j)] -= 0.5f * (p[idx(cube, i + 1, j)] - p[idx(cube, i - 1, j)]) * size;
            velocY[idx(cube, i, j)] -= 0.5f * (p[idx(cube, i, j + 1)] - p[idx(cube, i, j - 1)]) * size;
        }
    }
    setBound(cube, 1, velocX);
    setBound(cube, 2, velocY);
}

void FluidCube_advect(FluidCube *cube, int b, float *d, float *d0, float *velocX, float *velocY, float dt) {
    int size = cube->size;
    float dtx = dt * (size - 2);
    float dty = dt * (size - 2);

    for (int j = 1; j < size - 1; j++) {
        for (int i = 1; i < size - 1; i++) {
            float tmp1 = dtx * velocX[idx(cube, i, j)];
            float tmp2 = dty * velocY[idx(cube, i, j)];
            float x = i - tmp1;
            float y = j - tmp2;

            if (x < 0.5f) x = 0.5f;
            if (x > size - 0.5f) x = size - 0.5f;
            if (y < 0.5f) y = 0.5f;
            if (y > size - 0.5f) y = size - 0.5f;

            int i0 = (int)x;
            int i1 = i0 + 1;
            int j0 = (int)y;
            int j1 = j0 + 1;

            float s1 = x - i0;
            float s0 = 1.0f - s1;
            float t1 = y - j0;
            float t0 = 1.0f - t1;

            d[idx(cube, i, j)] = s0 * (t0 * d0[idx(cube, i0, j0)] + t1 * d0[idx(cube, i0, j1)]) + s1 * (t0 * d0[idx(cube, i1, j0)] + t1 * d0[idx(cube, i1, j1)]);
        }
    }
    setBound(cube, b, d);
}

void FluidCube_decay(FluidCube *cube, float value) {
    int size = cube->size;
    for (int j = 0; j < size; j++) {
        for (int i = 0; i < size; i++) {
            int idx_ = idx(cube, i, j);
            cube->density[idx_] = fmaxf(0, cube->density[idx_] - value);
        }
    }
}

void FluidCube_step(FluidCube *cube) {
    int size = cube->size;
    float visc = cube->visc;
    float diff = cube->diff;
    float dt = cube->dt;
    float *Vx = cube->Vx;
    float *Vy = cube->Vy;
    float *Vx0 = cube->Vx0;
    float *Vy0 = cube->Vy0;
    float *s = cube->s;
    float *density = cube->density;
    float kappa = cube->kappa;

    FluidCube_diffuse(cube, 1, Vx0, Vx, visc, dt, 20);
    FluidCube_diffuse(cube, 2, Vy0, Vy, visc, dt, 20);
    FluidCube_project(cube, Vx0, Vy0, Vx, Vy, 4);
    FluidCube_advect(cube, 1, Vx, Vx0, Vx0, Vy0, dt);
    FluidCube_advect(cube, 2, Vy, Vy0, Vx0, Vy0, dt);
    FluidCube_project(cube, Vx, Vy, Vx0, Vy0, 4);
    FluidCube_diffuse(cube, 0, s, density, diff, dt, 4);
    FluidCube_advect(cube, 0, density, s, Vx, Vy, dt);
    FluidCube_decay(cube, kappa);
}

int main() {
    int size = 100;
    FluidCube *cube = FluidCube_create(size, 0.0001, 0.0001, 0.1, 0.01);

    FluidCube_addDensity(cube, size / 2, size / 2, 100.0f);
    FluidCube_addVelocity(cube, size / 2, size / 2, 1.0f, 1.0f);

    for (int i = 0; i < 100; i++) {
        FluidCube_step(cube);
    }

    free(cube->s);
    free(cube->density);
    free(cube->Vx);
    free(cube->Vy);
    free(cube->Vx0);
    free(cube->Vy0);
    free(cube);

    return 0;
}