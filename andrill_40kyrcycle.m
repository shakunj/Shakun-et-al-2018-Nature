% This script calculates 10Be and 26Al concentrations in a bedrock profile
% at sea level and 2000 m asl in Antarctica for various exposure and 
% erosion scenarios over 40 kyr glacial cycles (as shown in Shakun et al., 
% 2018, Nature, Extended Data Figure 5). Nuclide production only occurs 
% when exposed, and erosion only occurs when ice covered. 

clear all
clc

%% at sea level
P10_0 = xlsread('production rate profile.xlsx','sea level','F2:F74076'); % sea-level high-latitude 10Be prod rate in 1 cm depth steps from 0 to 740 m
P26_0 = xlsread('production rate profile.xlsx','sea level','G2:G74076'); % sea-level high-latitude 26Al prod rate in 1 cm depth steps from 0 to 740 m
onoff = xlsread('andrill_onoff_40kyr.xlsx','onoff','C2:AB2501'); % matrix defining exposure (1) and burial (0) scenarios through 40 kyr cycles over 8 Myr
z = [0:70000]'; % depth range from 0 to 70000 cm (i.e., 0 to 700 m)

for er = 0:5:100 % erosion rate in mm/ky = m/My
    for perc = 0:5:100 % percent of 40 kyr glacial cycle that is ice free
        bedrock_10Be_0 = [zeros(70020,1); 0]; % initialize a bedrock profile with no 10Be
        bedrock_26Al_0 = [zeros(70020,1); 0]; % initialize a bedrock profile with no 26Al
            for t = 2:2:5000 % time steps in 2 ky intervals 
                bedrock_10Be_0(1:length(z),t/2+1) = onoff(t/2,perc/5+1)*P10_0(1+er/5:length(z)+er/5)*2000 + (bedrock_10Be_0(1+abs(onoff(t/2,perc/5+1)-1)*er/5:length(z)+abs(onoff(t/2,perc/5+1)-1)*er/5,t/2)*exp(-1*4.99e-7*2000));
                bedrock_26Al_0(1:length(z),t/2+1) = onoff(t/2,perc/5+1)*P26_0(1+er/5:length(z)+er/5)*2000 + (bedrock_26Al_0(1+abs(onoff(t/2,perc/5+1)-1)*er/5:length(z)+abs(onoff(t/2,perc/5+1)-1)*er/5,t/2)*exp(-1*9.83e-7*2000));
                %                                                production (if ice free)                                                          erosion (if ice covered)                                            decay 
            end
        summary_10Be_0(perc/5+1,er/5+1) = mean(bedrock_10Be_0(1,1500:2500)); % save average 10Be from 3 to 5 Myr into simulation and fill in a matrix by percent exposure during glacial cycles (y) versus erosion rate (x)
        summary_26Al_0(perc/5+1,er/5+1) = mean(bedrock_26Al_0(1,1500:2500)); % save average 26Al from 3 to 5 Myr into simulation and fill in a matrix by percent exposure during glacial cycles (y) versus erosion rate (x)
    end
end

%% at 2000 m asl
P10_2000 = xlsread('production rate profile.xlsx','2000 m asl','F2:F74076'); % 2000 m asl high-latitude 10Be prod rate in 1 cm depth steps from 0 to 740 m
P26_2000 = xlsread('production rate profile.xlsx','2000 m asl','G2:G74076'); % 2000 m asl high-latitude 26Al prod rate in 1 cm depth steps from 0 to 740 m

for er = 0:5:100 % erosion rate in mm/ky = m/My
    for perc = 0:5:100 % percent of 40 kyr glacial cycle that is ice free
        bedrock_10Be_2000 = [zeros(70020,1); 0]; % initialize a bedrock profile with no 10Be
        bedrock_26Al_2000 = [zeros(70020,1); 0]; % initialize a bedrock profile with no 26Al
            for t = 2:2:5000 % time steps in 2 ky intervals 
                bedrock_10Be_2000(1:length(z),t/2+1) = onoff(t/2,perc/5+1)*P10_2000(1+er/5:length(z)+er/5)*2000 + (bedrock_10Be_2000(1+abs(onoff(t/2,perc/5+1)-1)*er/5:length(z)+abs(onoff(t/2,perc/5+1)-1)*er/5,t/2)*exp(-1*4.99e-7*2000));
                bedrock_26Al_2000(1:length(z),t/2+1) = onoff(t/2,perc/5+1)*P26_2000(1+er/5:length(z)+er/5)*2000 + (bedrock_26Al_2000(1+abs(onoff(t/2,perc/5+1)-1)*er/5:length(z)+abs(onoff(t/2,perc/5+1)-1)*er/5,t/2)*exp(-1*9.83e-7*2000));
                %                                                production (if ice free)                                                          erosion (if ice covered)                                            decay 
            end
        summary_10Be_2000(perc/5+1,er/5+1) = mean(bedrock_10Be_2000(1,1500:2500)); % save average 10Be from 3 to 5 Myr into simulation and fill in a matrix by percent exposure during glacial cycles (y) versus erosion rate (x)
        summary_26Al_2000(perc/5+1,er/5+1) = mean(bedrock_26Al_2000(1,1500:2500)); % save average 26Al from 3 to 5 Myr into simulation and fill in a matrix by percent exposure during glacial cycles (y) versus erosion rate (x)
    end
end

% Plot results
x = [0:5:100]'; % erosion rates from 0 to 100 m/Myr
y = [0:5:100]'; % percent of exposure during glacial cycle
colormap(hot)
subplot(2,2,1)
[C,h] = contourf(x,y,log10(summary_10Be_0),[3:0.2:6]);
set(h,'LineColor','none')
xlabel('Erosion rate (m/Myr)'), ylabel('% of glacial cycle exposed'), title('10Be, 0 m asl')
subplot(2,2,2)
[C,h] = contourf(x,y,log10(summary_10Be_2000),[3:0.2:6]);
set(h,'LineColor','none')
xlabel('Erosion rate (m/Myr)'), ylabel('% of glacial cycle exposed'), title('10Be, 2000 m asl')
subplot(2,2,3)
[C,h] = contourf(x,y,log10(summary_26Al_0),[3:0.2:6]);
set(h,'LineColor','none')
xlabel('Erosion rate (m/Myr)'), ylabel('% of glacial cycle exposed'), title('26Al, 0 m asl')
subplot(2,2,4)
[C,h] = contourf(x,y,log10(summary_26Al_2000),[3:0.2:6]);
set(h,'LineColor','none')
xlabel('Erosion rate (m/Myr)'), ylabel('% of glacial cycle exposed'), title('26Al, 2000 m asl')