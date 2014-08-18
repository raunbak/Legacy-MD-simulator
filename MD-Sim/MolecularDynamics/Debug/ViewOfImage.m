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

z = 10;

disp('Creating image')
Image = HistogramToImageSlice(Hist,120,130);

disp('Adjusting gray scale')
Image = mat2gray(Image);
%%

Slice = Image(:,250/2-2:250/2+2);
Slice = sum(Slice,2);

plot(Slice)
