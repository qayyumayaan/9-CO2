%% Section1 : Import data.mat
close all
clear
clc


% Written by Ayaan Qayyum 8/23/22 - 9/3/22, 10/23/22. Equations modelled from https://pubmed.ncbi.nlm.nih.gov/14608002/


%% Importing Data
if ispc
    platform = "PC";
elseif ismac
    platform = "MAC";
elseif isunix
    platform = "UNIX";
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
file = append(pathname,infile);
load(file);

close all
% delete(findall(0));

%% Variables

CO2Pklocs = ietCO2;

VenousBloodVol = 4000; % ml
anatomicDeadSpaceVol = 150; % ml
VCO2 = 248; % ml min-1, exhaled volume of CO2
pulmonaryO2Uptake = 263; % ml min-1 
respiratoryQuotient = .9; % unitless, from paper. VCO2/VO2
arterialBloodVol = 1300; % ml
c = .1316; % mmHg

venousPartialCO2 = ones([9 1]);
arterialCO2 = ones([9 1]);
partialCO2PressurePerCompartment = ones([9 1]); % mmHg, Lung compartment k partial CO2 pressure. Necessary output!
CO2v = ones([9 1]);
CO2a = ones([9 1]);

numBreaths = length(CO2Pklocs);
NCO2 = zeros([9 numBreaths]);
estETCO2 = ones([1 numBreaths]);
lungCompartmentPartialCO2 = zeros([9 numBreaths]); % mmHg, per compartment k. Necessary output!

alveolarVolFRC = [6.58, 8.64, 10.11, 11.16, 11.90, 12.43, 12.81, 13.08, 13.27]'; % Alveolar vol. (% FRCk) while standing

tidalVol = 700*[4.58, 6.63, 8.48, 10.13, 11.62, 12.96, 14.17, 15.25, 16.22]'/100; % ml
functionalResidualCapacity = 3000*alveolarVolFRC/100; % ml 
lungCapillaryBloodVol = 200*alveolarVolFRC/100; % ml

g_k = ones([9 1]);
h_k = ones([9 1]);
weightFunction = ones([9 1]);


%% Computation

f = @(x) .53*(1.266-exp(-.0257*x));
% f_inv = @(x) -log(1.266 - (x/.53))/.0257;

k = 1:9;
g_k(k) = (.0263*k -.0205); % in the upright position

h_k(k) = .226*(1.102*exp(-.1063*k)); % in the upright position

weightFunction(k) = .10055*(1.36708 - exp(-.3393*k));

respiratoryInterval = 0; % used in CO2ProducedPerBreath

% numbered equations

% [CO_EST_ADJ] = COest(SBP,DBP,HR);
% CO_ADJ_ml = CO_EST_ADJ * 1000; 
% SVn = (mean(CO_EST_ml)./mean(HR)).*[.58; 3.21; 5.84; 8.47; 11.10; 13.73; 16.36; 18.99; 21.62]/100; % mL, Stroke volume per breath.

[StrokeVolPerBreath_adj] = SVnEst(SBP,DBP); % mL, Stroke volume per breath.
StrokeVolPerBreath_adj_ml = StrokeVolPerBreath_adj * 1000;
StrokeVolPerBreath = StrokeVolPerBreath_adj_ml.*[.58; 3.21; 5.84; 8.47; 11.10; 13.73; 16.36; 18.99; 21.62]/100;

for L = 1:numBreaths

    estETCO2inBreath = etCO2(L); % mmHg, end-tidal partial CO2 pressure

    C = venousPartialCO2.*StrokeVolPerBreath; % 2, ml * %

    A = arterialCO2.*StrokeVolPerBreath; % 3, ml * %

    CO2ProducedPerBreath = pulmonaryO2Uptake.*respiratoryQuotient.*(respiratoryInterval/60); % 4, ml, est. CO2 produced per breath

    CO2vn = venousPartialCO2 + (A + CO2ProducedPerBreath - C) / VenousBloodVol ; % 1, %

    D = f(partialCO2PressurePerCompartment) .* g_k .* StrokeVolPerBreath; % 6

    E = arterialCO2.*StrokeVolPerBreath; % 7, ml * %

    CO2an = arterialCO2 + (D-E)/arterialBloodVol; % 5, %

    F = f(partialCO2PressurePerCompartment).* weightFunction .* lungCapillaryBloodVol + c .* partialCO2PressurePerCompartment .* weightFunction .* functionalResidualCapacity + c .* estETCO2inBreath.* weightFunction.* anatomicDeadSpaceVol; % 8

    G = venousPartialCO2 .* StrokeVolPerBreath .* g_k; % 9

    a = f(estETCO2inBreath) ./ (c * estETCO2inBreath); % 10

    b = (weightFunction .* functionalResidualCapacity + h_k .* tidalVol) ./ (a.*((weightFunction .* lungCapillaryBloodVol + g_k .* StrokeVolPerBreath) + weightFunction .* functionalResidualCapacity + h_k .* tidalVol)); % 11

    CO2kn = b.*(F + G) ./ (weightFunction .* functionalResidualCapacity + h_k .* tidalVol); % 12, %

    PkCO2n = CO2kn.*c; % mmHg

    estETCO2n = sum(PkCO2n)*sum(h_k); % 13, mmHg

    % variables indexed to save per loop
    NCO2(:,L) = CO2kn(:);
    lungCompartmentPartialCO2(:,L) = PkCO2n(:);
    CO2v(:,L) = CO2vn(:);
    CO2a(:,L) = CO2an(:);
    estETCO2(:,L) = estETCO2n(:);

    fprintf("Breath %d computed. ",L)

    %PkCO2n./PkCO2n1 % value details if system is stable. 
    partialCO2PressurePerCompartment = PkCO2n;
    venousPartialCO2 = CO2vn;
    arterialCO2 = CO2an;

    respiratoryInterval = (CO2Pklocs(1,L)-respiratoryInterval)/1000; % s

end
disp("Finished.")
disp(lungCompartmentPartialCO2);


subplot(1,2,1)
plot(etCO2)
hold on
plot(estETCO2)
ylabel('End Tidal CO_2 (mmHg)')
xlabel('breath')
legend('Real ETCO2','Estimated ETCO2')

subplot(1,2,2)
plot(lungCompartmentPartialCO2')
title('Pressure in 9 compartments')
ylabel('P_{CO_2} (mmHg)')
xlabel('breath')

% CO2v 
% CO2a 
% estETCO2

%% Saving
saveResponse = questdlg("Would you like to SAVE?");
switch saveResponse
    case "Yes"
        saveDir = uigetdir(pathname,"Please select a directory to SAVE.");
        saveFile = erase(infile,".mat");
        saveName = append(saveFile,"_9CO2_",initials,".mat");
        disp("Saving in Progress. Please wait a moment.");
        if platform == "PC"
            nameSaveDir = append(saveDir,"\",saveName);
        end
        if platform == "MAC" || platform == "UNIX"
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