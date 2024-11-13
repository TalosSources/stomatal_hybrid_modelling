% Test script for photosynthesis_biochemical function

addpath('../matlab_code/t-and-c/T&C_file_model/T&C_CODE/');

% Define input values (example values provided)
Cc = 277.0525; % µmol mol⁻¹
IPAR = 282.8877; % µmol m⁻² s⁻¹
Csl = 384.4200; % µmol mol⁻¹
ra = 10.4263; % s m⁻¹
rb = 10.1525; % s m⁻¹
Ts = 18.5000; % °C
Pre = 843.3400; % Pa
Ds = 1360.6628; % kPa
Psi_L = -0.0255; % MPa
Psi_sto_50 = -2.5000; % MPa
Psi_sto_00 = -0.5000; % MPa
CT = 3; % Unitless
Vmax = 42; % µmol m⁻² s⁻¹
DS = 0.5; % kJ mol⁻¹
Ha = 72; % kJ mol⁻¹
FI = 0.0810; % Unitless
Oa = 210000; % µmol mol⁻¹
Do = 800; % kPa
a1 = 5; % Unitless
go = 0.01; % mol m⁻² s⁻¹
gmes = inf; % Unitless
rjv = 2.1; % Unitless

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
