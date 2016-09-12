function Tm = calculate_Tm_like_PrimerSelect(DNA)
%
% CALCULATE_TM_LIKE_PRIMERSELECT calculates the melting temperature of a 
% DNA oligos. 
%

%parameters
R = 1.987; %cal/K/mol
primer = 1E-6; % M --> use 1 uM primer concentration as default
Na = 0.05; %M
c = 0.1518; %empirical correction factor from fit

%DNA = {randseq(round(20 + (49-20).*rand))} %oligos of length (20-49)

prop = oligoprop(DNA, 'Salt', 0.05, 'Temp', 25, 'PrimerConc',primer);

%thermodynamic parameter
DH = prop.Thermo(1,1); %enthalpy H is in kcal/mol
DS = prop.Thermo(1,2); %entropy S is in cal/(K mol)
DG = prop.Thermo(1,3); %free energy G is in kcal/mol

%Formula is taken from vector NTI PDF file
Tm = DH/((DS -10.8)/1000 + R/1000*log(primer/4)) - 273.15 + 16.6*log10(Na); 

%Tm with correction; c was empirally determined
Tm = Tm + c;

