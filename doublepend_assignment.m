
%define mass, gravity, and lengths.
g = 9.81; % in m/s^2
L = 0.4325; % in m (L1 = L2 = L)
m = 0.339; % in kg (m1 = m2 = m)

% make time array
tArr = linspace(0,20,601); % create time arr for a 20 sec movie at 30 fps

% Solve the system of equations "doublepend" using ode45.
options = odeset('RelTol',1e-14,'AbsTol',[1e-14 1e-14 1e-14 1e-14]);
[t,q] = ode45(@doublepend,tArr,[pi/2 pi/2 0 0],options,L,g);

% export data for visualization
T = table(t,q);
writetable(T,'output.dat');