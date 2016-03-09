%% Overland Flow
% JSB 3/6/16

clear all
figure(1)
clf

%% Initialize

% Constants
%R = 100/(3600*24*365) % Constant everywhere Rainfall rate (m/s)
I =  50/(3600*24*365); % Infiltration rate (m/s)
n = .006; %roughness of surface

% Arrays
% time array setup
dt = 1; % time step (second)
tmax = 60*60*4; % max time (hours)
t = 0:dt:tmax; % time array

% distance array setup
dx = 50; % distance step (m)
xmax = 10000; % max dist (m) 
x = 0:dx:xmax; % dist array

% non-constant rainfall on hillslope
R = x; 
R(R<=500) = 100/(3600*24*365); % Everywhere less than 500m x gets this much rainfall (m/s)
R(R>500) = 60/(3600*24*365); % Everywhere more than 500m x gets this much rainfall (m/s)

% Initial overflow
Hnaught = 0; % H at time zero 
H = Hnaught * ones(size(x)); %initial H

% Initial topography 
zbmax = 200; % max bedrock
S = zbmax/xmax % Slope of initial bedrock linear
zb = zbmax-(S*x); % initial bedrock topo linear
%zb = -(.0003*(x.^2))+200; % initial bedrock topo
%S = diff(zb)/dx; % Slope of initial bedrock
z = zb+H; % Initial hill slope topo

% plot animation
n1=100; %number of plots
tplot = tmax/n1;

%% Analytical solution

%% Run


imax = length(t);
% In the time loop
for i = 1:imax;

Hedge=H(1:end-1) + 0.5*diff(H); % H at edge     
umean = (1/n)*(Hedge.^(2/3)).*(abs(S).^(1/2)); % mean water velocity

% calculate Q and diff Q
Q = umean.*Hedge; % water flux
Q = [0 Q Q(end)]; % fill in array
dQdx = diff(Q)/dx; % rate of change of flux

% rate of change water
dHdt = -dQdx+R-I; % rate of change of water height
H = H + (dHdt*dt); % water height
H = max(H,0); % cannot have negative water

% recalculate water + bedrock topo
z = zb+H; % total topography

%% PLOT
if (rem(t(i),tplot)==0) % decrease plotting frequency - speed up animation
figure(1)
plot(x,H,'-.','linewidth',3) % plot water height over time
hold off

   xlabel('Distance (m)','fontname','arial','fontsize',21) % x label
   ylabel('Height (m)','fontname','arial','fontsize',21) % y label
   set(gca,'fontsize',18,'fontname','arial') % axes number labels
   title(['Over flow evolution after ',num2str(t(i)),'seconds']) % title - accumulates model time
   axis([0 1000 0 .02]) % hold axes constant
   %set(gca,'YDIR','reverse')
   pause(.05)

end
end



