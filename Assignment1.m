m = 0.26*9.11e-31; % Effective mass of electron
height = 100e-9;
length = 200e-9;
num_particles = 100; % Initalize the number of particles within the region
T =300; % In kelvin
k = 1.3806e-23; % Boltzmann constant

% thermal velocity vth
vth = sqrt(2*k*T/m);

% mean free path
t = 0.2e-12 * vth;

% Initalize the vectors to zero in the x-y plane
state = zeros(num_particles,4);

% Assign each particle with fixed velocity given by vth