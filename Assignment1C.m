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

%diffusive is assigned to be 0,
diffusive = 0;
%Assigning a random velocity to each of the particles at the start
inital_xVelocity=randn(num_particles,1).*vth/sqrt(2);
inital_yVelocity=randn(num_particles,1).*vth/sqrt(2);
velAvg = sqrt(inital_xVelocity.^2 + inital_yVelocity.^2);
% Create an arry for previous x and y positons
prev_xPosition = inital_xPosition();
prev_yPosition = inital_yPosition();


% Simulation
for i = 1:1000
    % Updating particle location using Newton’s laws of motion
    inital_xPosition = inital_xPosition + inital_xVelocity.*delta_T;
    inital_yPosition = inital_yPosition + inital_yVelocity.*delta_T;
    
    % boundary condition where the particle reflects at the same angle in
    inital_xPosition(inital_xPosition<=0) = inital_xPosition(inital_xPosition<=0)+1;
    inital_xPosition(inital_xPosition>=1) = inital_xPosition(inital_xPosition>=1)-1;
    
    % Periodic boundary condition where the particle jumps to the opposite edge
    inital_yVelocity(inital_yPosition>=width)= -inital_yVelocity(inital_yPosition>=width);
    inital_yVelocity(inital_yPosition<=0) = -inital_yVelocity(inital_yPosition<=0);
    
    xposition = (inital_xPosition<1.3*length/2 & inital_xPosition>0.7*length/2);
    yposition =(inital_yPosition<width/3 | inital_yPosition>2*width/3);
    
    xPosyPos=(xposition & yposition);
    
    
    if(diffusive==1 & sum(xPosyPos)~=0)
        if((prev_xPosition(xPosyPos)<1.3*length/2 & prev_xPosition(xPosyPos)>0.7*length/2)) % is electrons position in between the bottles?
            inital_yVelocity(xPosyPos)=-inital_yVelocity(xPosyPos);
        else
            inital_xVelocity(xPosyPos)=-inital_xVelocity(xPosyPos);
        end
        
    elseif(diffusive==0 && sum(xPosyPos)~=0)
        if((prev_xPosition(xPosyPos)<1.3*length/2 & prev_xPosition(xPosyPos)>0.7*length/2) & inital_yVelocity(xPosyPos)>0)
            inital_yPosition(xPosyPos)=inital_yPosition(xPosyPos)-2*(inital_yPosition(xPosyPos)+2*width/3);
            inital_xVelocity(xPosyPos)=randn().*vth/sqrt(2);
            inital_yVelocity(xPosyPos)=-abs(randn().*vth/sqrt(2));
            
        elseif(inital_xVelocity(xPosyPos)>0)
            inital_xPosition(xPosyPos)=inital_xPosition(xPosyPos)-2*(inital_xPosition(xPosyPos)-0.7*length/2);
            inital_xVelocity(xPosyPos)=-abs(randn().*vth/sqrt(2));
            inital_yVelocity(xPosyPos)=randn().*vth/sqrt(2);
            
        elseif((prev_xPosition(xPosyPos)<1.3*length/2 & prev_xPosition(xPosyPos)>0.7*length/2) & inital_yVelocity(xPosyPos)<0)
            inital_yPosition(xPosyPos)=inital_yPosition(xPosyPos)+2*(width/3-inital_yPosition(xPosyPos));
            inital_xVelocity(xPosyPos)=randn().*vth/sqrt(2);
            inital_yVelocity(xPosyPos)=abs(randn().*vth/sqrt(2));
            
        else
            inital_xPosition(xPosyPos)=inital_xPosition(xPosyPos)+2*(1.3*length/2-inital_xPosition(xPosyPos));
            inital_xVelocity(xPosyPos)=abs(randn().*vth/sqrt(2));
            inital_yVelocity(xPosyPos)=abs(randn().*vth/sqrt(2));
            
        end
    end
    
    % exponential scattering probability for electrons
    prob = 1-exp(-delta_T/(0.2e-12));
    
    %Scatter electrons
    if (prob > rand())
        inital_xVelocity=randn(num_particles,1).*vth/sqrt(2);
        inital_yVelocity=randn(num_particles,1).*vth/sqrt(2);
        velAvg = sqrt(inital_xVelocity.^2 + inital_yVelocity.^2);
        collisionsNum=0;
        meanVelocity=0;
        collisionsNum=collisionsNum+1;
        diff=i-count;
        count=i;
        total=total+(diff*delta_T);
        meanVelocity=meanVelocity+(mean(velAvg));
        
        %Running Mean Free Path for part b
        averageMFP=(total/collisionsNum)*(meanVelocity/collisionsNum)
        %mean time between collisons
        meanCollisons = total/collisionsNum
        
    end
    
    figure(1)
    plot(inital_xPosition, inital_yPosition,'.','MarkerSize', 0.1)
    title('Movement of Electrons')
    xlim([0 length])
    ylim([0 width])
    
    % Setting up bottneck
    line([0.7*length/2 0.7*length/2],[w 2*width/3])
    line([1.3*length/2 1.3*length/2],[w 2*width/3])
    line([0.7*length/2 1.3*length/2],[width width])
    line([0.7*length/2 1.3*length/2],[2*width/3 2*width/3])
    line([0.7*length/2 0.7*length/2],[0 width/3])
    line([1.3*length/2 1.3*length/2],[0 width/3])
    line([0.7*length/2 1.3*length/2],[0 0])
    line([0.7*length/2 1.3*length/2],[width/3 width/3])
    hold on
    pause(0.01)
    
    
end


[xgradient,ygradient] = meshgrid(0:(length/10):length, 0:(width/10):width);

%This array will count the number of particles in each division of the frame
electronV=zeros(10,10);

temperatureV=zeros(10,10);
electrons_num=0;
velocityTot=0;

for i=1:10
    min_x=xgradient(1,i);
    max_x=xgradient(1,i+1);
    for k =1:10
        min_y=ygradient(k,1);
        max_y=ygradient(k+1,1);
        
        for j=1:num_particles
            %Check to see if particle is within this division of the frame
            if( (inital_yPosition(j)>min_y & inital_yPosition(j)<max_y) & (inital_xPosition(j)<max_x & inital_xPosition(j)>min_x))
                velocityTot=velocityTot+sqrt(inital_xVelocity(j)^2+inital_yVelocity(j)^2);
                electrons_num=electrons_num+1;
                electronV(i,k)=electronV(i,k)+1;
                
            end
        end
    end
end

%Plot the electron density map
figure(2)
surf(electronV)
title('Final Position Electron Density Map')
zlabel('Number of Electrons in Section')
hold on
pause(0.01)

%Plot the temperature map
figure(3)
surf(temperatureV)
title('Final Temperature Density Map')
zlabel('Temperature in Section')
hold on
pause(0.01)
