clear all;clc; close all;
disp('Data load')
data = importdata('TemperatureData.txt',',',1);
%%
Steps = 45;%105;
T = 0.093;
puredata = data.data;
dim = 2;
f_1 = figure;
hold on
plot(puredata(1:1000,dim),'.')


Total = [];
count = 0;
for i = 0:7500
    plot([round(Steps/2)+Steps*i round(Steps/2)+Steps*i]',[0 100],'g')
    if (i >= 150)
    Total = [Total  puredata(round(Steps/2)+Steps*i,dim)];
    count = count + 1;
    end
end
plot(puredata(:,1),T*ones(length(puredata(:,1)),1),'r');

%
%LowestPartOfPeriode = floor((Steps/2) + 0.5);
%
%for t = 1:24830
%    if (mod(Steps,t+ LowestPartOfPeriode-1) == 0 && mod(t,Steps) ~= 0)
%    plot([t t],[0 100])
%    end
%end
%axis([0 5*10^4 0 200])
AvgTemp = sum(Total) / count;
AvgTempstd = std(Total);


AvgTemp
AvgTempstd
procent_afvigelse = (T - AvgTemp) / T * 100
text(2.442e4,78,['Procent afvigelse ',num2str(procent_afvigelse),'%'])

%axis([2.43e4 2.455e4 -0.5 80])
xlabel('Tid [Skridt]');
ylabel('T [K]');
hold off

export_fig(f_1,'Temperatur','-pdf','-nocrop','-transparent')


%%
