clear all;clc;
disp('Running Simulation... not really just loading data')
Pdata = importdata('HistogramData.txt',',',1);
Vdata = importdata('VelocityHistogramData.txt',',',1);

%% Loading the data in.

pureP = Pdata.data;
pureV = Vdata.data;

disp('Creating histogram from data')
Hist = zeros( max(pureP(:,2))+1, max(pureP(:,3))+1, max(pureP(:,4))+1);
VHist = zeros( max(pureP(:,2))+1, max(pureP(:,3))+1, max(pureP(:,4))+1);

size(Hist)
n = 1;
for i = 1:size(Hist,1);
    for j = 1:size(Hist,2);
       for k = 1:size(Hist,3);
           
            Hist(i,j,k) = pureP(n);
            VHist(i,j,k) = pureV(n);
            n = n + 1;
       end
    end
end

z = 2;
ThinknessOfSlice = (size(Hist,3)/2)-z:(size(Hist,3)/2)+z;

PSlice = Hist(:,:,ThinknessOfSlice);
VSlice = VHist(:,:,ThinknessOfSlice);

%%
%Adding the slices in the z-dicektion.
PSliceSum = sum(PSlice,3);
VSliceSum = sum(VSlice,3);

Counts = zeros((size(Hist,1)/2),1); % List to hold the number of hits in a radius.
AvgV = zeros((size(Hist,1)/2),1);

radius = 0; % Start radius
NotDone = true; %Variable for the while loop.
NOP = 1000; % Number of data points for the circle data. (high means good circles)
theta = linspace(0, 2*pi, NOP); % angles for the ROI

i = 1; % Index for counts

while (NotDone) 
rho = ones(1, NOP) * radius;
[X1 Y1] = pol2cart(theta, rho); % Koordinater på kanten
X1_temp = X1 +(size(Hist,1)/2);   % Ryk cirklen til midten af histogrammet.
Y1_temp = Y1 + (size(Hist,1)/2);

a = roipoly((size(Hist,1)),(size(Hist,1)),X1_temp,Y1_temp);

radius = radius +1;


rho = ones(1, NOP) * radius;
[X1 Y1] = pol2cart(theta, rho); % Koordinater på kanten
X1_temp = X1 + (size(Hist,1)/2);   % Ryk cirklen til midten af histogrammet.
Y1_temp = Y1 + (size(Hist,1)/2);

b = roipoly((size(Hist,1)),(size(Hist,1)),X1_temp,Y1_temp);

ROI = b - a;

Counts(i) = sum(sum(ROI.*PSliceSum));
if (Counts(i) == 0) 
    AvgV(i) = 0;

else
    
    AvgV(i) = (sum(sum(ROI.*VSliceSum)) / Counts(i));
end


i = i +1;

if (radius == 100)
    NotDone = false;
end 

end

disp('Found em')

%%
f_1 = figure;
plot(AvgV,'xk')
xlabel('Radius (in blocks)');
ylabel('AvgV ');
%bar(Density)
export_fig(f_1,'VAvg','-pdf','-nocrop','-transparent')

