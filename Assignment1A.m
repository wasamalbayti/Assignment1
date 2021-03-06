m = 0.26*9.11e-31; % Effective mass of electron
length = 200e-9;
width = 100e-9;
num_particles = 20; % Initalize the number of particles within the region
T =300; % Temp in kelvin
k = 1.3806e-23; % Boltzmann constant

% thermal velocity vth
vth = sqrt(2*k*T/m);

% mean free path
t = 0.2e-12 * vth;

% Time step
delta_T = 1e-15;

% Assigning each particle a random location in the x ? y plane within the region
inital_xPosition = length * rand(num_particles,1);
inital_yPosition = width * rand(num_particles,1);

% Assign each particle with fixed velocity given by vth
angle = 2*pi*rand(num_particles,1); 
inital_xVelocity = ones(num_particles,1);
inital_yVelocity = ones(num_particles,1);
for i=1:num_particles
    inital_xVelocity(i) = vth*cos(angle(i));
    inital_yVelocity(i) = vth*sin(angle(i));
end

%simulate
for i=1:1000
    
    % Updating particle location using Newton�s laws of motion
    inital_xPosition = inital_xPosition + inital_xVelocity.*delta_T;
    inital_yPosition = inital_yPosition + inital_yVelocity.*delta_T;
    
    % boundary condition where the particle reflects at the same angle in 
    inital_xPosition(inital_xPosition<=0) = inital_xPosition(inital_xPosition<=0)+1;
    inital_xPosition(inital_xPosition>=1) = inital_xPosition(inital_xPosition>=1)-1;
    
    % Periodic boundary condition where the particle jumps to the opposite edge
    inital_yVelocity(inital_yPosition>=width)= -inital_yVelocity(inital_yPosition>=width);
    inital_yVelocity(inital_yPosition<=0) = -inital_yVelocity(inital_yPosition<=0);
    
    % Calculating the average temperature 
    temp = m*mean(sqrt((sum(abs(inital_xVelocity))/num_particles)^2+(sum(abs(inital_yVelocity))/num_particles)^2).^2)/k;
    
    
    %Plot the movement of the particles
    subplot(2,1,1)
    plot(inital_xPosition,inital_yPosition,'.')
    title('Particle trajectories')
    xlabel('X position')
    ylabel('Y position')
    xlim([0 length])
    ylim([0 width])
    hold on 
    pause(.01)
        
    % Temperature plot
    subplot(2,1,2)
    plot(i*delta_T,temp,'.')
    title('Temperature plot')
    ylabel('Temperature (k)')
    xlabel('Time (s)')
    hold on
end




