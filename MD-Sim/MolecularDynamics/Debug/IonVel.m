clear all;clc;
disp('Running Simulation... not really just loading data')
data = importdata('IonData.txt',',',2);

%%
puredata = data.data;

Ions = [];
length =  max((puredata(:,5))) + abs(min(puredata(:,5)));
procentLength = length/100;

z_width = procenLength*2; % 2% procent of max length is the slice.

Total_radius = 8.96312e-005;

IonsInSlice = 0;

for i = 1:size(puredata,4)
    
    if (puredata(i,5) <= z_width )
        IonsInSlice = IonsInSlice +1;
        
        Ions  = puredata(i,:);
    end
end