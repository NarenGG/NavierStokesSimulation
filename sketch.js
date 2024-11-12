
class FluidCube {
    constructor(size, diffusion, viscosity, dt, kappa) {
        this.size = size;
        this.dt = dt;
        this.diff = diffusion;
        this.visc = viscosity;
        this.kappa = kappa;

        const n = Math.pow(this.size, 3)

        this.s = new Array(n); for (let i=0; i<n; ++i) this.s[i] = 0;
        this.density = new Array(n); for (let i=0; i<n; ++i) this.density[i] = 0.01;

        this.Vx = new Array(n); for (let i=0; i<n; ++i) this.Vx[i] = 0;
        this.Vy = new Array(n); for (let i=0; i<n; ++i) this.Vy[i] = 0;
        this.Vz = new Array(n); for (let i=0; i<n; ++i) this.Vz[i] = 0;
        this.Vx0 = new Array(n); for (let i=0; i<n; ++i) this.Vx0[i] = 0;
        this.Vy0 = new Array(n); for (let i=0; i<n; ++i) this.Vy0[i] = 0;
        this.Vz0 = new Array(n); for (let i=0; i<n; ++i) this.Vz0[i] = 0;
    }

    idx(x, y, z) {
        return ((x) + (y) * this.size + (z) * Math.pow(this.size, 2));
    }

    addDensity(x, y, z, amount) {
        this.density[this.idx(x, y, z)] += amount;
    }

    addVelocity(x, y, z, amountX, amountY, amountZ) {
        let idx_ = this.idx(x, y, z)
        this.Vx[idx_] += amountX;
        this.Vy[idx_] += amountY;
        this.Vz[idx_] += amountZ;
    }

    setBound(b, x) {
        for(let j = 1; j < this.size - 1; j++) {
            for(let i = 1; i < this.size - 1; i++) {
                x[this.idx(i, j, 0  )] = b == 3 ? -x[this.idx(i, j, 1  )] : x[this.idx(i, j, 1  )];
                x[this.idx(i, j, this.size - 1)] = b == 3 ? -x[this.idx(i, j, this.size - 2)] : x[this.idx(i, j, this.size - 2)];
            }
        }
        for(let k = 1; k < this.size - 1; k++) {
            for(let i = 1; i < this.size - 1; i++) {
                x[this.idx(i, 0  , k)] = b == 2 ? -x[this.idx(i, 1  , k)] : x[this.idx(i, 1  , k)];
                x[this.idx(i, this.size - 1, k)] = b == 2 ? -x[this.idx(i, this.size - 2, k)] : x[this.idx(i, this.size - 2, k)];
            }
        }
        for(let k = 1; k < this.size - 1; k++) {
            for(let j = 1; j < this.size - 1; j++) {
                x[this.idx(0  , j, k)] = b == 1 ? -x[this.idx(1  , j, k)] : x[this.idx(1  , j, k)];
                x[this.idx(this.size - 1, j, k)] = b == 1 ? -x[this.idx(this.size - 2, j, k)] : x[this.idx(this.size - 2, j, k)];
            }
        }
        
        x[this.idx(0, 0, 0)]       = 0.33 * (x[this.idx(1, 0, 0)]
                                    + x[this.idx(0, 1, 0)]
                                    + x[this.idx(0, 0, 1)]);
        x[this.idx(0, this.size - 1, 0)]     = 0.33 * (x[this.idx(1, this.size - 1, 0)]
                                    + x[this.idx(0, this.size - 2, 0)]
                                    + x[this.idx(0, this.size - 1, 1)]);
        x[this.idx(0, 0, this.size - 1)]     = 0.33 * (x[this.idx(1, 0, this.size - 1)]
                                    + x[this.idx(0, 1, this.size - 1)]
                                    + x[this.idx(0, 0, this.size)]);
        x[this.idx(0, this.size - 1, this.size - 1)]   = 0.33 * (x[this.idx(1, this.size - 1, this.size - 1)]
                                    + x[this.idx(0, this.size - 2, this.size - 1)]
                                    + x[this.idx(0, this.size - 1, this.size - 2)]);
        x[this.idx(this.size - 1, 0, 0)]     = 0.33 * (x[this.idx(this.size - 2, 0, 0)]
                                    + x[this.idx(this.size - 1, 1, 0)]
                                    + x[this.idx(this.size - 1, 0, 1)]);
        x[this.idx(this.size - 1, this.size - 1, 0)]   = 0.33 * (x[this.idx(this.size - 2, this.size - 1, 0)]
                                    + x[this.idx(this.size - 1, this.size - 2, 0)]
                                    + x[this.idx(this.size - 1, this.size - 1, 1)]);
        x[this.idx(this.size - 1, 0, this.size - 1)]   = 0.33 * (x[this.idx(this.size - 2, 0, this.size - 1)]
                                    + x[this.idx(this.size - 1, 1, this.size - 1)]
                                    + x[this.idx(this.size - 1, 0, this.size - 2)]);
        x[this.idx(this.size - 1, this.size - 1, this.size - 1)] = 0.33 * (x[this.idx(this.size - 2, this.size - 1, this.size - 1)]
                                    + x[this.idx(this.size - 1, this.size - 2, this.size - 1)]
                                    + x[this.idx(this.size - 1, this.size - 1, this.size - 2)]);
    }

    linSolve(b, x, x0, a, c, iter) {
        const N = this.size;
        const cRecip = 1.0 / c;
        for (let k = 0; k < iter; k++) {
            for (let m = 1; m < N - 1; m++) {
                for (let j = 1; j < N - 1; j++) {
                    for (let i = 1; i < N - 1; i++) {
                        x[this.idx(i, j, m)] =
                            (x0[this.idx(i, j, m)]
                                + a*(    x[this.idx(i+1, j  , m  )]
                                        +x[this.idx(i-1, j  , m  )]
                                        +x[this.idx(i  , j+1, m  )]
                                        +x[this.idx(i  , j-1, m  )]
                                        +x[this.idx(i  , j  , m+1)]
                                        +x[this.idx(i  , j  , m-1)]
                               )) * cRecip;
                    }
                }
            }
            this.setBound(b, x);
        }
    }

    diffuse (b, x, x0, diff, dt, iter) {
        const N = this.size;
        const a = dt * diff * (N - 2) * (N - 2);
        this.linSolve(b, x, x0, a, 1 + 6 * a, iter);
    }

    project(velocX, velocY, velocZ, p, div, iter) {
        const N = this.size;
        for (let k = 1; k < N - 1; k++) {
            for (let j = 1; j < N - 1; j++) {
                for (let i = 1; i < N - 1; i++) {
                    div[this.idx(i, j, k)] = -0.5 * (
                        velocX[this.idx(i + 1, j, k)] - velocX[this.idx(i - 1, j, k)] +
                        velocY[this.idx(i, j + 1, k)] - velocY[this.idx(i, j - 1, k)] +
                        velocZ[this.idx(i, j, k + 1)] - velocZ[this.idx(i, j, k - 1)]
                    ) / N;
                    p[this.idx(i, j, k)] = 0;
                }
            }
        }
        this.setBound(0, div);
        this.setBound(0, p);
        this.linSolve(0, p, div, 1, 6, iter);
    
        for (let k = 1; k < N - 1; k++) {
            for (let j = 1; j < N - 1; j++) {
                for (let i = 1; i < N - 1; i++) {
                    velocX[this.idx(i, j, k)] -= 0.5 * (p[this.idx(i + 1, j, k)] - p[this.idx(i - 1, j, k)]) * N;
                    velocY[this.idx(i, j, k)] -= 0.5 * (p[this.idx(i, j + 1, k)] - p[this.idx(i, j - 1, k)]) * N;
                    velocZ[this.idx(i, j, k)] -= 0.5 * (p[this.idx(i, j, k + 1)] - p[this.idx(i, j, k - 1)]) * N;
                }
            }
        }
        this.setBound(1, velocX);
        this.setBound(2, velocY);
        this.setBound(3, velocZ);
    }

    advect(b, d, d0, velocX, velocY, velocZ, dt) {
        const N = this.size;
        let i0, i1, j0, j1, k0, k1;
        
        const dtx = dt * (N - 2);
        const dty = dt * (N - 2);
        const dtz = dt * (N - 2);
        
        let s0, s1, t0, t1, u0, u1;
        let tmp1, tmp2, tmp3, x, y, z;
        
        const Nfloat = N;
        let ifloat, jfloat, kfloat;
        
        for (let k = 1, kfloat = 1; k < N - 1; k++, kfloat++) {
            for (let j = 1, jfloat = 1; j < N - 1; j++, jfloat++) { 
                for (let i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
                    tmp1 = dtx * velocX[this.idx(i, j, k)];
                    tmp2 = dty * velocY[this.idx(i, j, k)];
                    tmp3 = dtz * velocZ[this.idx(i, j, k)];
                    x = ifloat - tmp1; 
                    y = jfloat - tmp2;
                    z = kfloat - tmp3;
                    
                    if (x < 0.5) x = 0.5; 
                    if (x > Nfloat + 0.5) x = Nfloat + 0.5; 
                    i0 = Math.floor(x); 
                    i1 = i0 + 1;
                    if (y < 0.5) y = 0.5; 
                    if (y > Nfloat + 0.5) y = Nfloat + 0.5; 
                    j0 = Math.floor(y);
                    j1 = j0 + 1; 
                    if (z < 0.5) z = 0.5;
                    if (z > Nfloat + 0.5) z = Nfloat + 0.5;
                    k0 = Math.floor(z);
                    k1 = k0 + 1;
                    
                    s1 = x - i0; 
                    s0 = 1.0 - s1; 
                    t1 = y - j0; 
                    t0 = 1.0 - t1;
                    u1 = z - k0;
                    u0 = 1.0 - u1;
                    
                    const i0i = i0;
                    const i1i = i1;
                    const j0i = j0;
                    const j1i = j1;
                    const k0i = k0;
                    const k1i = k1;
                    
                    d[this.idx(i, j, k)] = 
                        s0 * (t0 * (u0 * d0[this.idx(i0i, j0i, k0i)]
                                  + u1 * d0[this.idx(i0i, j0i, k1i)])
                            + t1 * (u0 * d0[this.idx(i0i, j1i, k0i)]
                                  + u1 * d0[this.idx(i0i, j1i, k1i)]))
                      + s1 * (t0 * (u0 * d0[this.idx(i1i, j0i, k0i)]
                                  + u1 * d0[this.idx(i1i, j0i, k1i)])
                            + t1 * (u0 * d0[this.idx(i1i, j1i, k0i)]
                                  + u1 * d0[this.idx(i1i, j1i, k1i)]));
                }
            }
        }
        this.setBound(b, d);
    }

    decay(value) {
        const N = this.size;
        for (let k = 0; k < N; k++) {
            for (let j = 0; j < N; j++) {
                for (let i = 0; i < N; i++) {
                    const idx = this.idx(i, j, k);
                    this.density[idx] = Math.max(0, this.density[idx] - value);
                }
            }
        }
    }

    step() {
        const N = this.size;
        const visc = this.visc;
        const diff = this.diff;
        const dt = this.dt;
        let Vx = this.Vx;
        let Vy = this.Vy;
        let Vz = this.Vz;
        let Vx0 = this.Vx0;
        let Vy0 = this.Vy0;
        let Vz0 = this.Vz0;
        let s = this.s;
        let density = this.density;
        let kappa = this.kappa;

        this.diffuse(1, Vx0, Vx, visc, dt, 4, N);
        this.diffuse(2, Vy0, Vy, visc, dt, 4, N);
        this.diffuse(3, Vz0, Vz, visc, dt, 4, N);
        
        this.project(Vx0, Vy0, Vz0, Vx, Vy, 4, N);
        
        this.advect(1, Vx, Vx0, Vx0, Vy0, Vz0, dt, N);
        this.advect(2, Vy, Vy0, Vx0, Vy0, Vz0, dt, N);
        this.advect(3, Vz, Vz0, Vx0, Vy0, Vz0, dt, N);
        
        this.project(Vx, Vy, Vz, Vx0, Vy0, 4, N);
        
        this.diffuse(0, s, density, diff, dt, 4, N);
        this.advect(0, density, s, Vx, Vy, Vz, dt, N);
    }

    isOnArcBoundary(x, y, z) {
        const dx = x - this.arcCenter[0];
        const dy = y - this.arcCenter[1];
        const dz = z - this.arcCenter[2];
        const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
        return Math.abs(dist - this.arcRadius) < 1; // Adjust for thickness tolerance
    }

    displayMiddleSlice() {
        const N = this.size;
        const mid = Math.floor(N / 2);

        for (let j = 0; j < N; j++) {
            for (let i = 0; i < N; i++) {
                const densityValue = this.density[this.idx(i, j, mid)];
                const colorValue = color(0, 0, map(densityValue, 0, 1, 0, 255));
                fill(colorValue);
                noStroke();
                rect(i * (width / N), j * (height / N), width / N, height / N);
            }
        }
    }

    addDensityCircle(radius, amount) {
        const N = this.size;
        const cellWidth = width / N;
        const cellHeight = height / N;
        const i = Math.floor(mouseX / cellWidth);
        const j = Math.floor(mouseY / cellHeight);
    
        if (i >= 0 && i < N && j >= 0 && j < N) {
            for (let k = 0; k < N; k++) {
                for (let y = -radius; y <= radius; y++) {
                    for (let x = -radius; x <= radius; x++) {
                        const dist = Math.sqrt(x * x + y * y);
                        if (dist <= radius) {
                            const ni = i + x;
                            const nj = j + y;
                            if (ni >= 0 && ni < N && nj >= 0 && nj < N) {
                                this.addDensity(ni, nj, k, amount);
                            }
                        }
                    }
                }
            }
        }
    }
}

let simCube;
const N = 40;

function setup() {
    createCanvas(400, 400);
    simCube = new FluidCube(N, 0.01, 0.5, 0.01, 0.001);
}

function draw() {
    background(220);
    if (mouseIsPressed) {
        simCube.addDensityCircle(2, 0.25);
    }

    simCube.step();
    simCube.displayMiddleSlice();
}