clear all;clc;
disp('Running Simulation... not really just loading data')
Pdata = importdata('HistogramData.txt',',',1);
Vdata = importdata('VelocityHistogramData.txt',',',1);

%% Loading the data in.

pureP = Pdata.data;
pureV = Vdata.data;

disp('Creating histogram from data')
Hist =zeros( max(pureP(:,2))+1, max(pureP(:,3))+1, max(pureP(:,4))+1);
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

z = 1;
ThinknessOfSlice = (size(Hist,1)/2)-z:(size(Hist,1)/2)+z;

%%
PSlice = Hist(125:126,125:126,:);
VSlice = VHist(125:126,125:126,:);

%%
TPSlice = sum(sum(PSlice,1));
TVSlice = sum(sum(VSlice,1));
disp('Did this');
%%
AvgV = zeros(0,size(PSlice,3));
for i = 1:size(PSlice,3)
   
    if(TPSlice(i) == 0)
        AvgV(i) = 0;
    else
        
        AvgV(i) = TVSlice(i) / TPSlice(i); 
        
    end
    
end

AvgV = (AvgV(size(Hist,3)/2 + 1:size(Hist,3)) + AvgV(size(Hist,3)/2:-1:1))./2;
%%
f_1 = figure;
plot(AvgV,'xk')
xlabel('Length of Crystal');
ylabel('AvgV ');
%bar(Density)
export_fig(f_1,'Time avg Vavg','-pdf','-nocrop','-transparent')

