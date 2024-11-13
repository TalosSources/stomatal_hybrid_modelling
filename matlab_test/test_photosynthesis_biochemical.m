% Test script for photosynthesis_biochemical function

addpath('../matlab_code/t-and-c/T&C_file_model/T&C_CODE/');

% Define input values (example values provided)
Cc = 200; % µmol mol⁻¹
IPAR = 1500; % µmol m⁻² s⁻¹
Csl = 350; % µmol mol⁻¹
ra = 100; % s m⁻¹
rb = 50; % s m⁻¹
Ts = 25; % °C
Pre = 100; % Pa
Ds = 1.5; % kPa
Psi_L = -1; % MPa
Psi_sto_50 = -1.5; % MPa
Psi_sto_00 = -3; % MPa
CT = 3; % Unitless
Vmax = 100; % µmol m⁻² s⁻¹
DS = 650; % kJ mol⁻¹
Ha = 60; % kJ mol⁻¹
FI = 0.15; % Unitless
Oa = 210; % µmol mol⁻¹
Do = 2.5; % kPa
a1 = 0.7; % Unitless
go = 0.05; % mol m⁻² s⁻¹
gmes = 1; % Unitless
rjv = 0.3; % Unitless

% Call the function
[CcF, An, rs, Rdark, F755nm, GAM, gsCO2] = photosynthesis_biochemical(...
    Cc, IPAR, Csl, ra, rb, Ts, Pre, Ds, ...
    Psi_L, Psi_sto_50, Psi_sto_00, ...
    CT, Vmax, DS, Ha, FI, Oa, Do, a1, go, gmes, rjv);

% Display the results
fprintf('Results:\n');
fprintf('CcF = %.2f\n', CcF);
fprintf('An = %.2f\n', An);
fprintf('rs = %.2f\n', rs);
fprintf('Rdark = %.2f\n', Rdark);
fprintf('F755nm = %.2f\n', F755nm);
fprintf('GAM = %.2f\n', GAM);
fprintf('gsCO2 = %.2f\n', gsCO2);
