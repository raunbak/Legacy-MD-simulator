clear all;clc;
disp('Running Simulation... not really just loading data')
data = importdata('HistogramData.txt',',',1);
%%

puredata = data.data;
%A = zeros(20, 10, 3);
%size(A)

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
%%
Counts = zeros(size(Hist,1)/2,1);
Howmany = zeros(size(Hist,1)/2,1);
r = 1;
n = 1;
while ( r < 125)
    
    
    for i = 125-r:r+125;
        for j = 125-r:r+125;
            for k = 125-r:r+125;
           
                d = sqrt( (i-125)^2 + (j-125)^2 + (k-125)^2 );
                if( d <= r  && d >= r - 1)
                    
                    Counts(n) = Counts(n) +  Hist(i,j,k) ;
                    
                    Howmany(n) = Howmany(n) + 1;
                    
                    
                end
                
            end
        end
    end

    %if ( n > 1)
    %Counts(n) = Counts(n) - Counts(n-1);
    %end
    n = n + 1;
    r = r + 1
    
end

%%

R = 1:size(Hist,1)/2;
r = 0:((size(Hist,1)/2)-1);
V = ((4/3) * pi * (R.^3 - r.^3))';

Density = Counts ./ V ;

%%
f_1 = figure;
hold on
plot(Density,'k')
plot(Density,'xk')
xlabel('Radius [bins]');
ylabel('Time aveage density [arb]');

plot([100.6742 100.6742],[0 40],'r')

%plot([94 94]',[0 40],'r')
axis([0 130 0 40])

hold off
export_fig(f_1,'Sphere','-pdf','-nocrop','-transparent')
