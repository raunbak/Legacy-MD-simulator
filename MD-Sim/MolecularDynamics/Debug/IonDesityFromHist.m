clear all;
radius = 99;

NOP = 1000%round(2*pi*radius);       % Antal pixels som kanten udgør
theta = linspace(0, 2*pi, NOP);
rho = ones(1, NOP) * radius;
[X1 Y1] = pol2cart(theta, rho); % Koordinater på kanten

X1_temp = X1 + 100;   % Ryk cirklen til midten af histogrammet.
Y1_temp = Y1 + 100;

a = roipoly(200,200,X1_temp,Y1_temp);

radius = 100;
NOP = 1000%round(2*pi*radius);       % Antal pixels som kanten udgør
theta = linspace(0, 2*pi, NOP);
rho = ones(1, NOP) * radius;
[X1 Y1] = pol2cart(theta, rho); % Koordinater på kanten

X1_temp = X1 + 100;   % Ryk cirklen til midten af histogrammet.
Y1_temp = Y1 + 100;

a = roipoly(200,200,X1_temp,Y1_temp)- a;