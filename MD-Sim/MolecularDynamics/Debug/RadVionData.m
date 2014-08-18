clear all;clc;
disp('Running Simulation... not really just loading data')
data = importdata('IonData.txt',',',2);
%%

puredata = data.data;


Vrad = sqrt(puredata(:,5).^2+puredata(:,6).^2);
Vtot = sqrt(puredata(:,5).^2+puredata(:,6).^2+puredata(:,7).^2);

plot(Vrad)
figure
plot(Vtot)