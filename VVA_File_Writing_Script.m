%close all
clc
clear

%load("C:\Users\amazi\Downloads\Mentorship\VVA_Template_2022_Labchart_8_V2_VVA001_VR_Data");

filterspec = '*.mat';
Title = 'Pick data.mat file with the VVA data:';
[infile,pathname] = uigetfile(filterspec,Title,"C:\Users\amazi\Downloads\Mentorship\VVA_Template_2022_Labchart_8_V2_VVA001_VR_Data");
FILE = append(pathname,infile);
load(FILE);

a = 179399;
b = 211049;
ab = b-a;

data = data_block1(:,a:b);
tick = ticktimes_block1(1,a:b);

[S, ~] = size(data);

for K = 1:S
    subplot(S,1,K);
    plot(tick,data(K,:));
end

context = ["BP","ECG","CO2"];
changelog = ["v1: VVA001 from 179399:211049. 8/19" "v2: Added EKG, BP, and MCA. 8/25" "v3: Removed MCA, added changelog and context block. Made this permanent script. Split data and bug fixes. Changed script save location. 8/25"];
%SAVEdata = data([2; 3; 5],:);

BP = data(2,:);
EKG = data(3,:);
CO2 = data(5,:);

%{
[S, ~] = size(SAVEdata);

for K = 1:S
    subplot(S,1,K);
    plot(tick,SAVEdata(K,:));
end
%}

%ticktimes = ticktimes_block1(1,179399:211049);
%data = data_block1([2,3,5],179399:211049);
save("G:\Study\9 CO2\VVA001_v3","tick","BP","EKG","CO2","context","changelog","-v7.3")
fprintf("Saved to %s!",pathname)