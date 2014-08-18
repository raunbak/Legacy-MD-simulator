clear all;clc;
disp('Running Simulation... not really just loading data')
data = importdata('HistogramData.txt',',',1);
%% Loading the data in.

puredata = data.data;

disp('Creating histogram from data')
Hist = zeros(250,250,250);
size(Hist)
n = 1;
for i = 1:size(Hist,1);
    for j = 1:size(Hist,2);
       for k = 1:size(Hist,3);
           
            Hist(i,j,k) = puredata(n);
            n = n + 1;
       end
    end
end


z = 2;
ThinknessOfSlice = (250/2)-z:(250/2)+z;

Slice = Hist(:,:,ThinknessOfSlice);

%% Finding the number of Ion in a radius. For high temp data will be lost.

%Adding the slices in the z-dicektion.
SliceSum = sum(Slice,3);

Counts = zeros( size(Hist,1)/2,1); % List to hold the number of hits in a radius.
Density = zeros( size(Hist,1)/2,1);


radius = 0; % Start radius
NotDone = true; %Variable for the while loop.
NOP = 1000; % Number of data points for the circle data. (high means good circles)
theta = linspace(0, 2*pi, NOP); % angles for the ROI

i = 1; % Index for counts

while (NotDone) 
rho = ones(1, NOP) * radius;
[X1 Y1] = pol2cart(theta, rho); % Koordinater på kanten
X1_temp = X1 + (size(Hist,1)/2);  % Ryk cirklen til midten af histogrammet.
Y1_temp = Y1 + (size(Hist,1)/2);

a = roipoly( size(Hist,1),size(Hist,1),X1_temp,Y1_temp);

radius = radius +1;


rho = ones(1, NOP) * radius;
[X1 Y1] = pol2cart(theta, rho); % Koordinater på kanten
X1_temp = X1 + (size(Hist,1)/2);   % Ryk cirklen til midten af histogrammet.
Y1_temp = Y1 + (size(Hist,1)/2);

b = roipoly(size(Hist,1),size(Hist,1),X1_temp,Y1_temp);

ROI = b - a;

Counts(i) = sum(sum(ROI.*SliceSum));
Density(i) = Counts(i) / (pi*( radius^2 -(radius-1)^2)*2*z);

i = i +1;

if (radius ==  (size(Hist,1)/2) )
    NotDone = false;
end 

end

disp('Found em')

%%
f_1 = figure;

bar(Counts)
xlabel('Radius (in blocks)');
ylabel('Sum of counts from historogram in ring');

export_fig(f_1,'Slice','-pdf','-nocrop','-transparent')

%%
f_2 = figure;

plot(Density(1:end),'xk')
xlabel('Radius');
ylabel('Time aveage density [arb]');
export_fig(f_2,'Density','-pdf','-nocrop','-transparent')


