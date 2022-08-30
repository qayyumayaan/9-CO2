%% Section1 : Import data.mat
close all
clear
clc

%{

Written by Ayaan Qayyum 8/23/22. Equations modelled from https://pubmed.ncbi.nlm.nih.gov/14608002/

%}

%% Importing Data
%load("/Users/qayyuma/Documents/9 CO2/VVA001_v3.mat")
load("G:\Study\9 CO2\VVA001_v3.mat");

%{

filterspec = '*.mat';
Title = 'Pick file generated from the Writing Script:';
[infile,pathname] = uigetfile(filterspec,Title,"C:\Users\amazi\Downloads\Mentorship\");
FILE = append(pathname,infile);
load(FILE);

%}

%% Variables

[CPks, locs] = BreathDef(CO2);
[HR] = heartrate(EKG);
[Q] = COest(BP,HR); % Cardiac output from blood pressure and heart rate

Vv = 5; % venous blood volume, L/min

V_D = 150; % ml, Anatomical dead space
RQ = .85; % Respiratory quotient, VCO2/VO2
Va = 1.3; % ml, Arterial blood volume
c = .1316; %/mmHg

%{
%x = 0;
%f_x = .53*(1.266-exp(-.0257*x));
%f_x_1 = -log(1.266 - (x/.53))/.0257;
%}

CO2vn1 = ones([9 1]); % %, Venous partial CO2 content
CO2an1 = ones([9 1]); % %, Arterial CO2 content
PkCO2n = 42*ones([9 1])/9; % mmHg, Lung compartment k partial CO2 pressure
PkCO2n1 = ones([9 1]); % mmHg, Lung compartment k partial CO2 pressure
CO2v = ones([9 1]);

[~, T] = size(locs);
NCO2 = zeros([9 T]);
PETCO2 = ones([9 T]);
PkCO2 = zeros([9 T]);

Trespn = 0; % s, Respiratory interval. See line ~124. 
%PETCO2n = CPks(1,L); % mmHg, end-tidal partial CO2 pressure
%PETCO2n1 = PETCO2n; % End-tidal partial CO2 pressure
VO2n = 62; % L/min, Pulmonary O2 Uptake

SVn = 70*[.58; 3.21; 5.84; 8.47; 11.10; 13.73; 16.36; 18.99; 21.62]/100; % mL, Stroke volume per breath. was just 70
V_Tn = 500*[4.58; 6.63; 8.48; 10.13; 11.62; 12.96; 14.17; 15.25; 16.22]/100; % ml, Tidal volume. was just .5
FRC = 3*[6.58; 8.64; 10.11; 11.16; 11.90; 12.43; 12.81; 13.08; 13.27]/100; % mL Functional residual capacity. was just 3
V_Cap = 75*[6.58; 8.64; 10.11; 11.16; 11.90; 12.43; 12.81; 13.08; 13.27]/100; % ml, Lung capillary blood volume. was just 75

%SVng_k = ones([9,1]) * 550; % mL, Stroke volume per breath, per segment
g_k = ones([9 1]);
h_k = ones([9 1]);
w_k = ones([9 1]);


%% Computation

for k = 1:9
    g_k(k,1) = -.0205 + .0263*k;
end

for k = 1:9
    h_k(k,1) = .266*(1.102*exp(-.1063*k));
end

for k = 1:9
    w_k(k,1) = .1055*(1.36708 - exp(-.3393*k));
end

% numbered equations

%for L = 1:3
for L = 1:T

PETCO2n1 = CPks(1,L); % mmHg, end-tidal partial CO2 pressure

C = CO2vn1.*SVn; % 2

A = CO2an1.*SVn; % 3

B = VO2n.*RQ.*(Trespn/60); % 4, est. CO2 produced per breath

CO2vn = CO2vn1 + (A + B - C) / Vv ; % 1

D = sum(.53*(1.266-exp(-.0257*PkCO2n1))) .* times(g_k,SVn); % 6

E = CO2an1.*SVn; % 7

CO2an = CO2an1 + (D-E)/Va; % 5

F = .53*(1.266-exp(-.0257*PkCO2n1)).* w_k .* V_Cap + c .* PkCO2n1 .* w_k .* FRC + c .* PETCO2n1.* w_k.* V_D; % 8

G = CO2vn1 .* SVn .* g_k; % 9

a = .53*(1.266-exp(-.0257*PETCO2n1)) ./ (c * PETCO2n1); % 10

b = (w_k .* FRC + h_k .* V_Tn) ./ (a.*((w_k .* V_Cap + g_k .* SVn) + w_k .* FRC + h_k .* V_Tn)); % 11

CO2kn = b.*(F + G) ./ (w_k .* FRC + h_k .* V_Tn); %12

PkCO2n = CO2kn./c;

PETCO2n = PkCO2n * sum(h_k); % 13

NCO2(:,L) = CO2kn(:,1);
PkCO2(:,L) = PkCO2n(:,1);
CO2v(:,L) = CO2vn(:,1);
CO2a(:,L) = CO2an(:,1);
PETCO2(:,L) = PETCO2n(:,1);

fprintf("Breath %d computed. ",L)

PkCO2n1 = PkCO2n;
CO2vn1 = CO2vn;
CO2an1 = CO2an;
PETCO2n1 = PETCO2n;
Trespn = (locs(1,L)-Trespn)/1000; % s, Respiratory interval

end
disp("Finished.")
%PkCO2 = PkCO2
%CO2v = CO2v
%CO2a = CO2a
%PETCO2 = PETCO2


% End Tidal CO2 records the CO2 output. One breath is from one defined peak to another.


%% SBP and DBP Calculator (Liljestrand & Zander formula)
function [CO_EST] = COest(BP,HR) 

[SBPpeaks, SBPlocation] = findpeaks(BP,"MinPeakProminence",20,"MinPeakHeight",110,"MinPeakDistance",10); 
    % findpeaks(BP,"MinPeakProminence",20,"MinPeakHeight",110,"MinPeakDistance",10); 
SBP = sum(SBPpeaks)/size(SBPlocation,2);

BP_DBP = -1 * BP;
[DBPpeaks, DBPlocation] = findpeaks(BP_DBP,"MinPeakProminence",15,"MinPeakHeight",-80,"MinPeakDistance",100);
    % findpeaks(BP_DBP,"MinPeakProminence",15,"MinPeakHeight",-80,"MinPeakDistance",100)
DBP = -sum(DBPpeaks)/size(DBPlocation,2);
PP = SBP - DBP;

[CO_EST] = HR * PP / (SBP + DBP);

%disp("CO_EST")
end

%% Heart Rate Calculator

function [HR] = heartrate(EKG)
Sekg = smoothdata(EKG(1,:),"movmedian",30);
[~, Locs] = findpeaks(Sekg,"MinPeakProminence",20,"MinPeakWidth",10,"MinPeakHeight",100);
%findpeaks(Sekg,"MinPeakProminence",20,"MinPeakWidth",10,"MinPeakHeight",100)
o1 = size(Locs,2);
v = zeros([1 o1-1]);

for o = 1:(o1-1) % makes a set v that contains heart bpm calculated between peak distances. 
    v(1,Locs(1,o):Locs(1,o+1)) = 1000 * 60 / (Locs(o+1)-Locs(o)); % turns peak distances into bpm from seconds
end

%[Heart_Rate] = Prominent_Peaks(:,1:(size(Prominent_Peaks,2)-1)); % line is needed for plotting against time bc it makes sure the sets are the same length
[HR] = v;
end


%% finding the beginning of a breath
function [CPks, locs] = BreathDef(CO2)

smoothed_data = smoothdata(CO2(1,:),"movmedian",10);
%findpeaks(smoothed_data,"MinPeakHeight",38,"MinPeakProminence",2,"MinPeakDistance",2500); % for testing the smoothing only
diff_smoothed = diff(smoothed_data);
smooooothed = smoothdata(diff_smoothed,"movmedian",200);
[~, N] = size(smooooothed);

for P = 1:N
    if smooooothed(1,P) > 0
        smooooothed(1,P) = 0;
    end
end

smooooothed = -1*smooooothed;
twoderiv = diff(smooooothed);
[~, N1] = size(twoderiv);

for P1 = 1:N1
    if twoderiv(1,P1) < 0
        twoderiv(1,P1) = 0;
    end
end

twoderiv = smoothdata(twoderiv,"movmedian",5);
[~, locs] = findpeaks(twoderiv,"MinPeakHeight",.0001,"MinPeakProminence",.002,"MinPeakDistance",90);

CPks = CO2(1,locs);

%{
plot(twoderiv)
findpeaks(twoderiv,"MinPeakHeight",.0001,"MinPeakProminence",.002,"MinPeakDistance",90);
%}

end

