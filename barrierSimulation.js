
let simCube;
let dyeGenerator;
let thetaSlider;
let dtSlider;
const thetaInitial = -25 * Math.PI / 180;
const N = 100;
let arcBarrier;

function setup() {
    createCanvas(400, 400);
    thetaSlider = createSlider(-Math.PI * 2, Math.PI * 2, thetaInitial, Math.PI / 180);
    thetaSlider.position(10, 10);
    thetaSlider.size(80);
    dtSlider = createSlider(0.001, 1, 0.1, 0.01);
    dtSlider.position(10, 30);
    dtSlider.size(80);

    simCube = new FluidCube(N, 0, 0.0000000000000000000000000001, dtSlider.value(), 0.0); // Increase viscosity and reduce diffusion to reduce dispersion
    dyeGenerator = new DyeGenerator(simCube, Math.round(N / 3), N / 2, 2, 0.25); // High velocity
    arcBarrier = new ArcBarrier(N / 2, N / 2, 10, 1, -Math.PI / 2, Math.PI / 2); // Example semi-circle arc barrier
}

function draw() {
    background(220);

    theta = thetaSlider.value();
    simCube.dt = dtSlider.value();

    dyeGenerator.generate(theta);

    arcBarrier.reflectVelocity(simCube);

    if (mouseIsPressed) {
        //simCube.addDensityCircle(2, 0.25);
    }

    simCube.step();
    simCube.display();
    arcBarrier.display(simCube);
}