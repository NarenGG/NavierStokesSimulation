class FluidCube {
    constructor(size, diffusion, viscosity, dt, kappa) {
        this.size = size;
        this.dt = dt;
        this.diff = diffusion;
        this.visc = viscosity;
        this.kappa = kappa;

        const n = Math.pow(this.size, 2);

        this.s = new Array(n).fill(0);
        this.density = new Array(n).fill(0); // Start with no initial density

        this.Vx = new Array(n).fill(0);
        this.Vy = new Array(n).fill(0);
        this.Vx0 = new Array(n).fill(0);
        this.Vy0 = new Array(n).fill(0);
    }

    idx(x, y) {
        return x + y * this.size;
    }

    addDensity(x, y, amount) {
        this.density[this.idx(x, y)] += amount;
    }

    addVelocity(x, y, amountX, amountY) {
        let idx_ = this.idx(x, y);
        this.Vx[idx_] += amountX;
        this.Vy[idx_] += amountY;
    }

    setBound(b, x) {
        for (let i = 1; i < this.size - 1; i++) {
            x[this.idx(i, 0)] = b == 2 ? -x[this.idx(i, 1)] : x[this.idx(i, 1)];
            x[this.idx(i, this.size - 1)] = b == 2 ? -x[this.idx(i, this.size - 2)] : x[this.idx(i, this.size - 2)];
        }
        for (let j = 1; j < this.size - 1; j++) {
            x[this.idx(0, j)] = b == 1 ? -x[this.idx(1, j)] : x[this.idx(1, j)];
            x[this.idx(this.size - 1, j)] = b == 1 ? -x[this.idx(this.size - 2, j)] : x[this.idx(this.size - 2, j)];
        }

        x[this.idx(0, 0)] = 0.5 * (x[this.idx(1, 0)] + x[this.idx(0, 1)]);
        x[this.idx(0, this.size - 1)] = 0.5 * (x[this.idx(1, this.size - 1)] + x[this.idx(0, this.size - 2)]);
        x[this.idx(this.size - 1, 0)] = 0.5 * (x[this.idx(this.size - 2, 0)] + x[this.idx(this.size - 1, 1)]);
        x[this.idx(this.size - 1, this.size - 1)] = 0.5 * (x[this.idx(this.size - 2, this.size - 1)] + x[this.idx(this.size - 1, this.size - 2)]);
    }

    linSolve(b, x, x0, a, c, iter) {
        const N = this.size;
        const cRecip = 1.0 / c;
        for (let k = 0; k < iter; k++) {
            for (let j = 1; j < N - 1; j++) {
                for (let i = 1; i < N - 1; i++) {
                    x[this.idx(i, j)] =
                        (x0[this.idx(i, j)]
                            + a * (x[this.idx(i + 1, j)]
                                + x[this.idx(i - 1, j)]
                                + x[this.idx(i, j + 1)]
                                + x[this.idx(i, j - 1)]
                            )) * cRecip;
                }
            }
            this.setBound(b, x);
        }
    }

    diffuse(b, x, x0, diff, dt, iter) {
        const N = this.size;
        const a = dt * diff * (N - 2) * (N - 2);
        this.linSolve(b, x, x0, a, 1 + 4 * a, iter);
    }

    project(velocX, velocY, p, div, iter) {
        const N = this.size;
        for (let j = 1; j < N - 1; j++) {
            for (let i = 1; i < N - 1; i++) {
                div[this.idx(i, j)] = -0.5 * (
                    velocX[this.idx(i + 1, j)] - velocX[this.idx(i - 1, j)] +
                    velocY[this.idx(i, j + 1)] - velocY[this.idx(i, j - 1)]
                ) / N;
                p[this.idx(i, j)] = 0;
            }
        }
        this.setBound(0, div);
        this.setBound(0, p);
        this.linSolve(0, p, div, 1, 4, iter);

        for (let j = 1; j < N - 1; j++) {
            for (let i = 1; i < N - 1; i++) {
                velocX[this.idx(i, j)] -= 0.5 * (p[this.idx(i + 1, j)] - p[this.idx(i - 1, j)]) * N;
                velocY[this.idx(i, j)] -= 0.5 * (p[this.idx(i, j + 1)] - p[this.idx(i, j - 1)]) * N;
            }
        }
        this.setBound(1, velocX);
        this.setBound(2, velocY);
    }

    advect(b, d, d0, velocX, velocY, dt) {
        const N = this.size;
        let i0, i1, j0, j1;

        const dtx = dt * (N - 2);
        const dty = dt * (N - 2);

        let s0, s1, t0, t1;
        let tmp1, tmp2, x, y;

        const Nfloat = N;
        let ifloat, jfloat;

        for (let j = 1, jfloat = 1; j < N - 1; j++, jfloat++) {
            for (let i = 1, ifloat = 1; i < N - 1; i++, ifloat++) {
                tmp1 = dtx * velocX[this.idx(i, j)];
                tmp2 = dty * velocY[this.idx(i, j)];
                x = ifloat - tmp1;
                y = jfloat - tmp2;

                if (x < 0.5) x = 0.5;
                if (x > Nfloat + 0.5) x = Nfloat + 0.5;
                i0 = Math.floor(x);
                i1 = i0 + 1;
                if (y < 0.5) y = 0.5;
                if (y > Nfloat + 0.5) y = Nfloat + 0.5;
                j0 = Math.floor(y);
                j1 = j0 + 1;

                s1 = x - i0;
                s0 = 1.0 - s1;
                t1 = y - j0;
                t0 = 1.0 - t1;

                const i0i = i0;
                const i1i = i1;
                const j0i = j0;
                const j1i = j1;

                d[this.idx(i, j)] =
                    s0 * (t0 * d0[this.idx(i0i, j0i)] + t1 * d0[this.idx(i0i, j1i)]) +
                    s1 * (t0 * d0[this.idx(i1i, j0i)] + t1 * d0[this.idx(i1i, j1i)]);
            }
        }
        this.setBound(b, d);
    }

    decay(value) {
        const N = this.size;
        for (let j = 0; j < N; j++) {
            for (let i = 0; i < N; i++) {
                const idx = this.idx(i, j);
                this.density[idx] = Math.max(0, this.density[idx] - value);
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
        let Vx0 = this.Vx0;
        let Vy0 = this.Vy0;
        let s = this.s;
        let density = this.density;
        let kappa = this.kappa;

        this.diffuse(1, Vx0, Vx, visc, dt, 20); // Increase iterations for more stability
        this.diffuse(2, Vy0, Vy, visc, dt, 20); // Increase iterations for more stability

        this.project(Vx0, Vy0, Vx, Vy, 4);

        this.advect(1, Vx, Vx0, Vx0, Vy0, dt);
        this.advect(2, Vy, Vy0, Vx0, Vy0, dt);

        this.project(Vx, Vy, Vx0, Vy0, 4);

        this.diffuse(0, s, density, diff, dt, 4);
        this.advect(0, density, s, Vx, Vy, dt);
        this.decay(kappa);
    }

    display() {
        const N = this.size;

        for (let j = 0; j < N; j++) {
            for (let i = 0; i < N; i++) {
                const densityValue = this.density[this.idx(i, j)];
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
            for (let y = -radius; y <= radius; y++) {
                for (let x = -radius; x <= radius; x++) {
                    const dist = Math.sqrt(x * x + y * y);
                    if (dist <= radius) {
                        const ni = i + x;
                        const nj = j + y;
                        if (ni >= 0 && ni < N && nj >= 0 && nj < N) {
                            this.addDensity(ni, nj, amount);
                        }
                    }
                }
            }
        }
    }
}

class DyeGenerator {
    constructor(fluidCube, x, y, dyeAmount, velocity) {
        this.fluidCube = fluidCube;
        this.x = x;
        this.y = y;
        this.dyeAmount = dyeAmount;
        this.velocity = velocity;
    }

    generate(angle) {
        const velocityX = this.velocity * Math.cos(angle);
        const velocityY = this.velocity * Math.sin(angle);

        this.fluidCube.addDensity(this.x, this.y, this.dyeAmount);
        this.fluidCube.addVelocity(this.x, this.y, velocityX, velocityY); // Adjust velocity
    }
}

class ArcBarrier {
    constructor(centerX, centerY, radius, thickness, startAngle, endAngle) {
        this.centerX = centerX;
        this.centerY = centerY;
        this.radius = radius;
        this.thickness = thickness;
        this.startAngle = startAngle;
        this.endAngle = endAngle;
    }

    isInside(x, y) {
        const dx = x - this.centerX;
        const dy = y - this.centerY;
        const dist = Math.sqrt(dx * dx + dy * dy);
        const angle = Math.atan2(dy, dx);
        return dist >= this.radius - this.thickness && dist <= this.radius + this.thickness && angle >= this.startAngle && angle <= this.endAngle;
    }

    reflectVelocity(fluidCube) {
        const N = fluidCube.size;
        for (let j = 0; j < N; j++) {
            for (let i = 0; i < N; i++) {
                if (this.isInside(i, j)) {
                    const idx = fluidCube.idx(i, j);
                    const dx = i - this.centerX;
                    const dy = j - this.centerY;
                    const dist = Math.sqrt(dx * dx + dy * dy);
                    const nx = dx / dist;
                    const ny = dy / dist;

                    const dotProduct = fluidCube.Vx[idx] * nx + fluidCube.Vy[idx] * ny;
                    fluidCube.Vx[idx] -= 2 * dotProduct * nx;
                    fluidCube.Vy[idx] -= 2 * dotProduct * ny;
                }
            }
        }
    }

    display(fluidCube) {
        const N = fluidCube.size;
        noStroke();
        fill(255); // Red color with some transparency
        for (let j = 0; j < N; j++) {
            for (let i = 0; i < N; i++) {
                if (this.isInside(i, j)) {
                    rect(i * (width / N), j * (height / N), width / N, height / N);
                }
            }
        }
    }
}

/*
let simCube;
let dyeGenerator;
let thetaSlider;
const thetaInitial = -25 * Math.PI / 180;
const N = 100;
let arcBarrier;

function setup() {
    createCanvas(400, 400);
    thetaSlider = createSlider(-Math.PI * 2, Math.PI * 2, thetaInitial, Math.PI / 180);
    thetaSlider.position(10, 10);
    thetaSlider.size(80);

    simCube = new FluidCube(N, 0, 0.0000000000000000000000000001, 0.1, 0.0); // Increase viscosity and reduce diffusion to reduce dispersion
    dyeGenerator = new DyeGenerator(simCube, Math.round(N / 3), N / 2, 2, 0.25); // High velocity
    arcBarrier = new ArcBarrier(N / 2, N / 2, 10, 1, -Math.PI / 2, Math.PI / 2); // Example semi-circle arc barrier
}

function draw() {
    background(220);

    theta = thetaSlider.value();
    console.log(theta);

    dyeGenerator.generate(theta);

    arcBarrier.reflectVelocity(simCube);

    if (mouseIsPressed) {
        //simCube.addDensityCircle(2, 0.25);
    }

    simCube.step();
    simCube.display();
    arcBarrier.display(simCube);
}
*/