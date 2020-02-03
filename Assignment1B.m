m = 0.26*9.11e-31; % Effective mass of electron
length = 200e-9;
width = 100e-9;
num_particles = 20; % Initalize the number of particles within the region
T =300; % Temp in kelvin
k = 1.3806e-23; % Boltzmann constant

% thermal velocity vth
vth = sqrt(2*k*T/m)

% mean free path
t = 0.2e-12 * vth

% Time step
delta_T = 1e-15;

% Assigning each particle a random location in the x ? y plane within the region
inital_xPosition = length * rand(num_particles,1);
inital_yPosition = width * rand(num_particles,1);


%Part 2
%Assigning a random velocity to each of the particles at the start
inital_xVelocity=randn(num_particles,1).*vth/sqrt(2);
inital_yVelocity=randn(num_particles,1).*vth/sqrt(2);
velAvg = sqrt(inital_xVelocity.^2 + inital_yVelocity.^2);

total=0;
collisionsNum=0;
meanVelocity=mean(velAvg);


for i=1:1000
    
    % Updating particle location using Newton’s laws of motion
    inital_xPosition = inital_xPosition + inital_xVelocity.*delta_T;
    inital_yPosition = inital_yPosition + inital_yVelocity.*delta_T;
    
    % boundary condition where the particle reflects at the same angle in
    inital_xPosition(inital_xPosition<=0) = inital_xPosition(inital_xPosition<=0)+1;
    inital_xPosition(inital_xPosition>=1) = inital_xPosition(inital_xPosition>=1)-1;
    
    % Periodic boundary condition where the particle jumps to the opposite edge
    inital_yVelocity(inital_yPosition>=width)= -inital_yVelocity(inital_yPosition>=width);
    inital_yVelocity(inital_yPosition<=0) = -inital_yVelocity(inital_yPosition<=0);
    
    % exponential scattering probability for electrons
    prob = 1-exp(-delta_T/(0.2e-12));
    
    %Scatter electrons
    if (prob > rand())
        inital_xVelocity=randn(num_particles,1).*vth/sqrt(2);
        inital_yVelocity=randn(num_particles,1).*vth/sqrt(2);
        velAvg = sqrt(inital_xVelocity.^2 + inital_yVelocity.^2);
        
        collisionsNum=collisionsNum+1;
        diff=i-count;
        count=i;
        total=total+(diff*delta_T);
        meanVelocity=meanVelocity+(mean(velAvg));
        
        %Running Mean Free Path
        averageMFP=(total/collisionsNum)*(meanVelocity/collisionsNum)
        %mean time between collisons
        meanCollisons = total/collisionsNum
        
        
    end
    figure(1)
    subplot(3,1,1)
    plot(inital_xPosition,inital_yPosition,'.','MarkerSize', 0.1)
    title('Movement of Electrons')
    xlim([0 length])
    ylim([0 width])
    hold on
    pause(.01)
    
    Temp=m*mean(velAvg.^2)/k;
    subplot(3,1,2)
    plot(i*delta_T,Temp,'x')
    title('Temperature ');
    ylabel('Temperature')
    xlabel('Time')
    hold on
    
    subplot(3,1,3)
    histogram(velAvg,5)
    title('Electron Speed')
    xlabel('Speed (m/s)')
    ylabel('count')
    hold on
    
end



