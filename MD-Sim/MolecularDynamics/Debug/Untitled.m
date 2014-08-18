clear all;clc;
disp('Running Simulation... not really just loading data')
data = importdata('IonData.txt',',',2);
%%

puredata = data.data;

zdistance =  max((puredata(:,5))) / 4;

Slice = [];
n = 1;
for i = 1:size(puredata,1);
    if ( abs(puredata(i,5)) <= zdistance )
        Slice(n,:) = puredata(i,:);
        n = n+1;
    end
end

maxradius = 5.3244e-5; 
%sqrt(max(abs(puredata(:,3)))^2 + max(abs(puredata(:,4)))^2);
blocks = 2*100;

NumberofIonsInBlock = zeros(blocks,1);
density = zeros(blocks,1);
for k = 1:blocks;
    for i = 1:size(Slice,1);
        d = sqrt(puredata(i,3)^2 + puredata(i,4)^2);
        if ( d < k * maxradius/(blocks/2) && d > (k-1)*maxradius/(blocks/2))
            NumberofIonsInBlock(k) = NumberofIonsInBlock(k) + 1;
        end
        
    end
    Rmax = k * maxradius/(blocks/2);
    Rmin = (k-1)*maxradius/(blocks/2);
    density(k) = NumberofIonsInBlock(k) / (pi*( Rmax^2 -Rmin^2)*2*zdistance);
end

radiuslist = linspace(1/blocks,2,blocks);

%%
f_1 = figure;
plot(radiuslist,density)
%bar(radiuslist,NumberofIonsInBlock)
export_fig(f_1,'CutImage','-pdf','-nocrop','-transparent')
%plot(NumberofIonsInBlock)
%((int) ((ions[i].Position(0))/PixelToDistance+((double) HistNx)/2))