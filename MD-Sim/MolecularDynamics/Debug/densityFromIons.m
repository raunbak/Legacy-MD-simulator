clear all;clc;
disp('Running Simulation... not really just loading data')
data = importdata('IonData.txt',',',2);
%%

puredata = data.data;

%%
radiusOfCrystal = 8.96312e-5; 
zdistance =  (2*radiusOfCrystal / 100) * 5; % zdistance will 

Slice = [];
n = 1;
for i = 1:size(puredata,1);
    if ( abs(puredata(i,5)) <= zdistance )
        Slice(n,:) = puredata(i,:);
        n = n+1;
    end
end


dist = sqrt(Slice(:,3).^2 + Slice(:,4).^2);
maxdistance = max(dist);


radii = linspace(0,10*10^-5,100);
NumberOfIonsInR = zeros(size(radii,2),1);


for i = 2:size(radii,2)
    
    for n = 1:size(Slice,1)
        d = sqrt(Slice(n,3)^2 + Slice(n,4)^2);
        if ( d <= radii(i) && d > radii(i-1) )
            
            NumberOfIonsInR(i) = NumberOfIonsInR(i) + 1;
            
        end    
    end
end

%%

V = (pi*( radii(2:end).^2 - radii(1:end-1).^2 )* zdistance)';
%%
Mass = NumberOfIonsInR * 40 * 1.66053878283e-27; % I Kg.
Density = Mass(2:end) ./ V;

%%

plot(radii(2:end),Density)

