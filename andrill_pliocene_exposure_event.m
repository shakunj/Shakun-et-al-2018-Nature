% This script calculates 10Be and 26Al concentrations in a bedrock profile
% at sea level and 2000 m asl for the past 5.3 Myr assuming a single
% mid-Pliocene exposure event (as shown in Shakun et al., 2018, Nature, 
% Extended Data Figure 6). The user sets the duration of the exposure
% event and the glacial erosion rate.

clear all
clc
er = 100; % erosion rate in mm/ky = m/My
dose = 200; % length of exposure event in ky
P10_0 = xlsread('production rate profile.xlsx','sea level','F2:F74076'); % sea level high-latitude 10Be prod rate in 1 cm depth steps from 0 to 740 m
P10_2000 = xlsread('production rate profile.xlsx','2000 m asl','F2:F74076'); % 2000 m asl high-latitude 10Be prod rate in 1 cm depth steps from 0 to 740 m
z = [0:70000]'; % depth range from 0 to 70000 cm (i.e., 0 to 700 m)
age = [5300:-2:0]';

bedrock_10Be_0 = [zeros(70020,1); 0]; % initialize a bedrock profile with no 10Be
bedrock_10Be_2000 = [zeros(70020,1); 0]; % initialize a bedrock profile with no 10Be

for t = 2:2:1698 % time steps in 2 ky, starting at 5300 ka, covered
    bedrock_10Be_0(1:length(z),t/2+1) = 0*P10_0(1+er/5:length(z)+er/5)*2000 + (bedrock_10Be_0(1+er/5:length(z)+er/5,t/2)*exp(-1*4.99e-7*2000));
    bedrock_10Be_2000(1:length(z),t/2+1) = 0*P10_2000(1+er/5:length(z)+er/5)*2000 + (bedrock_10Be_2000(1+er/5:length(z)+er/5,t/2)*exp(-1*9.83e-7*2000));
                %                                             no production                                erosion                    decay 
end
for t = 1700:2:1700+dose-2 % time steps in 2 ky, exposure starting at 3600 ka (AND-1B Pliocene diatomite inter(x)val)
    bedrock_10Be_0(1:length(z),t/2+1) = 1*P10_0(1+er/5:length(z)+er/5)*2000 + (bedrock_10Be_0(1:length(z),t/2)*exp(-1*4.99e-7*2000));
    bedrock_10Be_2000(1:length(z),t/2+1) = 1*P10_2000(1+er/5:length(z)+er/5)*2000 + (bedrock_10Be_2000(1:length(z),t/2)*exp(-1*9.83e-7*2000));
                %                                                production                             no erosion                    decay 
end
for t = 1700+dose:2:5300 % time steps in 2 ky, cover(x)ed until 0 ka
    bedrock_10Be_0(1:length(z),t/2+1) = 0*P10_0(1+er/5:length(z)+er/5)*2000 + (bedrock_10Be_0(1+er/5:length(z)+er/5,t/2)*exp(-1*4.99e-7*2000));
    bedrock_10Be_2000(1:length(z),t/2+1) = 0*P10_2000(1+er/5:length(z)+er/5)*2000 + (bedrock_10Be_2000(1+er/5:length(z)+er/5,t/2)*exp(-1*9.83e-7*2000));
                %                                             no production                                erosion                    decay 
end
surface_10Be_0 = bedrock_10Be_0(1,:)'; % time series of surface 10Be evolution
surface_10Be_2000 = bedrock_10Be_2000(1,:)'; % time series of surface 10Be evolution

% Average the simulated time series to the resolution of the measured 
% AND-1B record
sample_ages = [320	384 % upper and lower ages of AND-1B samples in ka
420     422
530     648
680     1410
1650	2930
2970	3120
3190	5300];

index = 2651-sample_ages/2;
index = flipud(index);
index = fliplr(index);

surface_10Be_0_binned = NaN(length(age),1);
surface_10Be_2000_binned = NaN(length(age),1);
for z = 1:length(index) % averages of simulated 10Be concentrations over the time intervals spanned by each sample
    surface_10Be_0_binned(index(z,1):index(z,2)) = mean(surface_10Be_0(index(z,1):index(z,2)));
    surface_10Be_2000_binned(index(z,1):index(z,2)) = mean(surface_10Be_2000(index(z,1):index(z,2)));
end

% Plot results 
semilogy(age,surface_10Be_0_binned)
hold on
semilogy(age,surface_10Be_2000_binned)
xlabel('Age (ka)'), ylabel('10Be (atoms/g)'), title('10Be conc with a Pliocene exposure event')