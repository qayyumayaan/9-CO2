%% Section1 : Import data.mat
close all
clear
clc

%{

Written by Ayaan Qayyum 8/23/22 - 9/3/22. Equations modelled from https://pubmed.ncbi.nlm.nih.gov/14608002/

%}

%% Importing Data
%load("/Users/qayyuma/Documents/9 CO2/VVA001_v3.mat")
%load("C:\Users\amazi\Documents\GitHub\9-CO2\VVA001_v3.mat");
%load("C:\Users\amazi\Downloads\Mentorship\VVA014  Template 2022 Labchart 8 V8 RollTilt TRANS0.18VS_clean.mat");
%load("/Users/qayyuma/Downloads/VVA014  Template 2022 Labchart 8 V8 RollTilt TRANS0.18VS_clean.mat");

if ispc
    Platform = "PC";
elseif ismac
    Platform = "MAC";
elseif isunix
    Platform = "UNIX";
end

initials = char(inputdlg("Please enter your initials."));
if initials == ""
    return
end
filterspec = '*.mat';
Title = 'Pick file generated from VVA Cleandata:';
[infile,pathname] = uigetfile(filterspec,Title);
if pathname == 0
    return
end
FILE = append(pathname,infile);
load(FILE);

close all
delete(findall(0));

%% Variables

CO2Pklocs = ietCO2;

Vv = 4000; % ml, venous blood volume, from cardiac output
V_D = 150; % ml, Anatomical dead space
VCO2 = 248; % ml min-1, exhaled volume of CO2
VO2n = 263; % ml min-1, Pulmonary O2 Uptake. Q * (ConcO2art - ConcO2ven)
RQ = .9; % Respiratory quotient, taken from paper. VCO2/VO2
Va = 1300; % ml, Arterial blood volume
c = .1316; %/mmHg

CO2vn1 = ones([9 1]); % %, Venous partial CO2 content
CO2an1 = ones([9 1]); % %, Arterial CO2 content
PkCO2n1 = ones([9 1]); % mmHg, Lung compartment k partial CO2 pressure. Necessary output!
CO2v = ones([9 1]);
CO2a = ones([9 1]);

BN = length(CO2Pklocs); % # of breaths
NCO2 = zeros([9 BN]);
PETCO2 = ones([1 BN]);
PkCO2 = zeros([9 BN]); % mmHg, Lung compartment k partial CO2 pressure. Necessary output!

Trespn = 0; % s, Respiratory interval. See line ~129. 

V_Tn = 700*[4.58; 6.63; 8.48; 10.13; 11.62; 12.96; 14.17; 15.25; 16.22]/100; % ml, Tidal volume. was just .5
FRC = 3000*[6.58; 8.64; 10.11; 11.16; 11.90; 12.43; 12.81; 13.08; 13.27]/100; % mL Functional residual capacity. was just 3
V_Cap = 200*[6.58; 8.64; 10.11; 11.16; 11.90; 12.43; 12.81; 13.08; 13.27]/100; % ml, Lung capillary blood volume. was just 75

g_k = ones([9 1]);
h_k = ones([9 1]);
w_k = ones([9 1]);


%% Computation

f = @(x) .53*(1.266-exp(-.0257*x));
%f_inv = @(x) -log(1.266 - (x/.53))/.0257;

for k = 1:9
    g_k(k,1) = -.0205 + .0263*k;
end

for k = 1:9
    h_k(k,1) = .226*(1.102*exp(-.1063*k));
end

for k = 1:9
    w_k(k,1) = .10055*(1.36708 - exp(-.3393*k));
end

% numbered equations

%[CO_EST_ADJ] = COest(SBP,DBP,HR);
%CO_ADJ_ml = CO_EST_ADJ * 1000; 
%SVn = (mean(CO_EST_ml)./mean(HR)).*[.58; 3.21; 5.84; 8.47; 11.10; 13.73; 16.36; 18.99; 21.62]/100; % mL, Stroke volume per breath.

[SV_adj] = SVnEst(SBP,DBP); % mL, Stroke volume per breath.
SV_adj_ml = SV_adj * 1000;
SVn = SV_adj_ml.*[.58; 3.21; 5.84; 8.47; 11.10; 13.73; 16.36; 18.99; 21.62]/100;

for L = 1:BN

    PETCO2n1 = etCO2(L); % mmHg, end-tidal partial CO2 pressure

    C = CO2vn1.*SVn; % 2, ml * %

    A = CO2an1.*SVn; % 3, ml * %

    B = VO2n.*RQ.*(Trespn/60); % 4, ml, est. CO2 produced per breath

    CO2vn = CO2vn1 + (A + B - C) / Vv ; % 1, %

    D = f(PkCO2n1) .* g_k .* SVn; % 6

    E = CO2an1.*SVn; % 7, ml * %

    CO2an = CO2an1 + (D-E)/Va; % 5, %

    F = f(PkCO2n1).* w_k .* V_Cap + c .* PkCO2n1 .* w_k .* FRC + c .* PETCO2n1.* w_k.* V_D; % 8

    G = CO2vn1 .* SVn .* g_k; % 9

    a = f(PETCO2n1) ./ (c * PETCO2n1); % 10

    b = (w_k .* FRC + h_k .* V_Tn) ./ (a.*((w_k .* V_Cap + g_k .* SVn) + w_k .* FRC + h_k .* V_Tn)); % 11

    CO2kn = b.*(F + G) ./ (w_k .* FRC + h_k .* V_Tn); % 12, %

    PkCO2n = CO2kn.*c; % mmHg

    PETCO2n = sum(PkCO2n)*sum(h_k); % 13, mmHg

    NCO2(:,L) = CO2kn(:,1);
    PkCO2(:,L) = PkCO2n(:,1);
    CO2v(:,L) = CO2vn(:,1);
    CO2a(:,L) = CO2an(:,1);
    PETCO2(1,L) = PETCO2n(:,1);

    fprintf("Breath %d computed. ",L)

    %PkCO2n./PkCO2n1 % value details if system is stable. 
    PkCO2n1 = PkCO2n;
    CO2vn1 = CO2vn;
    CO2an1 = CO2an;


    Trespn = (CO2Pklocs(1,L)-Trespn)/1000; % s, Respiratory interval

end
disp("Finished.")
disp(PkCO2);


subplot(1,2,1)
plot(etCO2)
hold on
plot(PETCO2)
ylabel('End Tidal CO_2 (mmHg)')
xlabel('breath')
legend('Real ETCO2','Estimated ETCO2')

subplot(1,2,2)
plot(PkCO2')
title('Pressure in 9 compartments')
ylabel('P_{CO_2} (mmHg)')
xlabel('breath')

%CO2v 
%CO2a 
%PETCO2

%% Saving
saveResponse = questdlg("Would you like to SAVE?");
switch saveResponse
    case "Yes"
        saveDir = uigetdir(pathname,"Please select a directory to SAVE.");
        saveFile = erase(infile,".mat");
        saveName = append(saveFile,"_9CO2_",initials,".mat");
        disp("Saving in Progress. Please wait a moment.");
        if Platform == "PC"
            nameSaveDir = append(saveDir,"\",saveName);
        end
        if Platform == "MAC"
            nameSaveDir = append(saveDir,"/",saveName);
        end
        save(nameSaveDir);
    otherwise
        return
end

msgscriptend = append("File saved as ", saveName," at ", saveDir, ". Thank you for using this script! ");
End = msgbox(msgscriptend);

% End Tidal CO2 records the CO2 output. One breath is from one defined peak to another.
%% SV Approximation via the Liljestrand & Zander formula, irrespective of Heart Rate. 
function [SV_adj] = SVnEst(SBP,DBP) 

    SBPavg = sum(SBP)/length(SBP);

    DBPavg = sum(DBP)/length(DBP);

    PPavg = SBPavg - DBPavg; % Pulse Pressure

    SV = (PPavg / (SBPavg + DBPavg)); % mL, Stroke volume per breath, as defined in paper: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5317099/
    SV_adj = SV / 3.548;

end

%% OLD SBP and DBP Calculator (Initial Implementation of the Liljestrand & Zander formula)
function [CO_EST_ADJ] = COest(SBP,DBP,HR) 

    SBPavg = sum(SBP)/length(SBP);

    DBPavg = sum(DBP)/length(DBP);

    PPavg = SBPavg - DBPavg; % Pulse Pressure

    CO_EST = HR * PPavg / (SBPavg + DBPavg);
    CO_EST_ADJ = CO_EST / 3.548;

end