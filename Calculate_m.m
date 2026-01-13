function [ m ] = Calculate_m(Temperature, lambda, substance)
% contributors over the years:
% Daniel Jakubczyk, Anastasiya Derkachova, Tho Do Duc, Gennadiy Derkachov, Kwasi Nyandey, Sima Alikhanzadeh-Arani 

% substance:
% pentaethylene glycol - PEG,
% hexaethylene glycol - HEG,
% dipropylene glycol (mixture of isomers) - DPG,
% tripropylene glycol - TPG,
% ethylene glycol - EG,
% diethylene glycol - DEG,
% triethylene glycol - TEG,
% tetraethylene glycol - TetEG,
% propylene glycol - PG,
% Glycerol,
% fused silica - fSiO2,
% BK7 glass,
% EG-Rhodamine 6G mixture -  EG+6G,
% DMSO,
% Benzonitryl,
% MetOH,
% Canola oil,
% polystyrene microspheres - uPS,
% H2O,
% titania nanoparticles - nanoTiO2,
% Au (in bulk),
% Ag (in bulk)

% lambda in nm
% Temperature in deg Celsius

T_ref = 22;

switch substance
    
case 'PEG'   
    % Derkachova et al., https://doi.org/10.1016/j.optmat.2026.117874
    x=lambda/1e3;
    A=1.59771;
    B_IR=0.00957;
    C_IR=4.12834;
    B_UV=0.51063;
    C_UV=0.02021;
    A_T = -3.631e-4;
    B_T = -4.78e-6;
    C_T = -0.25; 
    m=sqrt(A+B_IR*x^2/(x^2-C_IR)+B_UV*x^2/(x^2-C_UV)) +(Temperature-T_ref)*(A_T+B_T/(x+C_T));

case 'HEG'
   % Derkachova et al., https://doi.org/10.1016/j.optmat.2026.117874
   x=lambda/1e3;
   A=1.64934; 
   B_IR=0.00528;
   C_IR=3.00459;
   B_UV=0.46337;
   C_UV=0.02221;
   A_T = -3.486e-4; 
   B_T = -5.24e-6; 
   C_T = -0.243; 
   m=sqrt(A+B_IR*x^2/(x^2-C_IR)+B_UV*x^2/(x^2-C_UV)) +(Temperature-T_ref)*(A_T+B_T/(x+C_T));

case 'DPG'
    % Derkachova et al., https://doi.org/10.1016/j.optmat.2026.117874
    x=lambda/1e3;
    A=1.64629;
    B_IR=0.02138;
    C_IR=6.07253;
    B_UV=0.40305;
    C_UV=0.02319;
    A_T = -3.367e-4;
    B_T = -3.84e-6; 
    C_T = -0.296; 
    m=sqrt(A+B_IR*x^2/(x^2-C_IR)+B_UV*x^2/(x^2-C_UV)) +(Temperature-T_ref)*(A_T+B_T/(x+C_T));

 case 'TPG'
    % Derkachova et al., https://doi.org/10.1016/j.optmat.2026.117874
    x=lambda/1e3;
    A=1.60003; 
    B_IR=0.00401;
    C_IR=2.61436;
    B_UV=0.45704;
    C_UV=0.0213;
    A_T = - 3.512e-4;   
    B_T = -4.54e-6; 
    C_T = -0.267; 
    m=sqrt(A+B_IR*x^2/(x^2-C_IR)+B_UV*x^2/(x^2-C_UV)) +(Temperature-T_ref)*(A_T+B_T/(x+C_T));
    
case 'EG' 
    % Jakubczyk et al., https://doi.org/10.1038/s41597-023-02819-3, https://doi.org/10.1038/s41597-024-03547-y
    x=lambda/1e3;
    A=1.34238;
    B_IR=0.0137;
    C_IR=3.12;
    B_UV=0.68263;
    C_UV=0.01332;
    A_T = -2.643e-4;
    B_T = -7.5e-6;
    C_T = -0.17; 
    m=sqrt(A+B_IR*x^2/(x^2-C_IR)+B_UV*x^2/(x^2-C_UV)) +(Temperature-T_ref)*(A_T+B_T/(x+C_T));
    
case 'DEG'
    % Jakubczyk et al., https://doi.org/10.1038/s41597-023-02819-3, https://doi.org/10.1038/s41597-024-03547-y
    x=lambda/1e3;
    A=1.60169;
    B_IR=0.02833;
    C_IR=6.07658;
    B_UV=0.46338;
    C_UV=0.02;
    A_T = -3.132e-4; 
    B_T = -5.61e-6; 
    C_T = -0.245;
    m=sqrt(A+B_IR*x^2/(x^2-C_IR)+B_UV*x^2/(x^2-C_UV)) +(Temperature-T_ref)*(A_T+B_T/(x+C_T));

case 'TEG'
    % Jakubczyk et al., https://doi.org/10.1038/s41597-023-02819-3, https://doi.org/10.1038/s41597-024-03547-y
    x=lambda/1e3;
    A=1.12259;
    B_IR=0.00475;
    C_IR=2.38854;
    B_UV=0.96913;
    C_UV=0.0108;
   A_T = -3.122e-4; 
    B_T = -6.3e-6; 
    C_T = -0.22;
    m=sqrt(A+B_IR*x^2/(x^2-C_IR)+B_UV*x^2/(x^2-C_UV)) +(Temperature-T_ref)*(A_T+B_T/(x+C_T));

case 'TetEG'
    % Jakubczyk et al., https://doi.org/10.1038/s41597-023-02819-3, https://doi.org/10.1038/s41597-024-03547-y
    x=lambda/1e3;
    A=1.48424;
    B_IR=0.00354;
    C_IR=2.25089;
    B_UV=0.61562;
    C_UV=0.01689;   
    A_T = -3.570e-4; 
    B_T = -6.51e-6; 
    C_T = -0.235; 
    m=sqrt(A+B_IR*x^2/(x^2-C_IR)+B_UV*x^2/(x^2-C_UV)) +(Temperature-T_ref)*(A_T+B_T/(x+C_T));

case 'PG'
    % Jakubczyk et al., https://doi.org/10.1038/s41597-023-02819-3, https://doi.org/10.1038/s41597-024-03547-y
    x=lambda/1e3;
    A=1.15131;
    B_IR=0.00824;
    C_IR=2.7529;
    B_UV=0.87613;
    C_UV=0.01092;
    A_T = -3.027e-4;  
    B_T = -1.2e-6; 
    C_T = -0.72; 
    m=sqrt(A+B_IR*x^2/(x^2-C_IR)+B_UV*x^2/(x^2-C_UV)) +(Temperature-T_ref)*(A_T+B_T/(x+C_T));
      
case 'Glycerol'
    % Jakubczyk et al., https://doi.org/10.1038/s41597-023-02819-3, https://doi.org/10.1038/s41597-024-03547-y
    x=lambda/1e3;
    A=1.6062;
    B_IR=0.06222;
    C_IR=7.98826;
    B_UV=0.53803;
    C_UV=0.0181;
    A_T = -2.395e-4;
    B_T = -6.2e-6;
    C_T = -0.18; 
    m=sqrt(A+B_IR*x^2/(x^2-C_IR)+B_UV*x^2/(x^2-C_UV)) +(Temperature-T_ref)*(A_T+B_T/(x+C_T));
    
case 'fSiO2'
    % fused silica
    % three term temperature-dependent effective Sellmeier model from
    % Leviton et al.
    TKelvin = Temperature+273.15;
    lambda = lambda/1000;

    S1 = 1.10127 - TKelvin * 4.94251e-5 + (TKelvin^2) * 5.27414e-7 - (TKelvin^3) * 1.597e-9 + (TKelvin^4) * 1.75949e-12;
    S2 = 1.78752e-5 + TKelvin * 4.76391e-5 - (TKelvin^2) * 4.49019e-7 + (TKelvin^3) * 1.44546e-9 - (TKelvin^4) * 1.57223e-12;
    S3 = 7.93552e-1 - TKelvin * 1.27815e-3 + (TKelvin^2) * 1.84595e-5 - (TKelvin^3) * 9.20275e-8 + (TKelvin^4) * 1.48829e-10;
    
    lambda1 = -8.906e-2 + TKelvin * 9.0873e-6 - (TKelvin^2) * 6.53638e-8 + (TKelvin^3) * 7.77072e-11 + (TKelvin^4) * 6.84605e-14;
    lambda2 = 2.97562e-1 - TKelvin * 8.59578e-4 + (TKelvin^2) * 6.59069e-6 - (TKelvin^3) * 1.09482e-8 + (TKelvin^4) * 7.85145e-13;
    lambda3 = 9.34454 - TKelvin * 7.09788e-3 + (TKelvin^2) * 1.01968e-4 - (TKelvin^3) * 5.07660e-7 + (TKelvin^4) * 8.21348*10e-10;
    
    m = (1 + (S1*lambda^2)/(lambda^2 - lambda1^2) + (S2*lambda^2)/(lambda^2 - lambda2^2) + (S3*lambda^2)/(lambda^2 - lambda3^2))^0.5;

case 'BK7'
    % dispersion parameters
    B1=1.03961212;
    B2=0.231792344;
    B3=1.01046945;
    C1=6.00069867e-3;
    C2=2.00179144e-2;
    C3=103.560653;
    lambda = lambda/1000;
    
    m = sqrt((B1*lambda^2)/(lambda^2-C1)+(B2*lambda^2)/(lambda^2-C2)+(B3*lambda^2)/(lambda^2-C3)+1);
   
case 'EG+6G'
        % Sellmeier dispersion parameters EG+temperature dependence of solution EG+6G:      
        x=lambda/1e3;
        A=1.34238;
        B_IR=0.0137;
        C_IR=3.12;
        B_UV=0.68263;
        C_UV=0.01332;
        DeltaT = 3.435e-4;
        
        m=sqrt(A+B_IR*x^2/(x^2-C_IR)+B_UV*x^2/(x^2-C_UV)) - DeltaT*(Temperature-T_ref);              
    
% case 'DMI' % very approximate
%     A0 = 1.4664;
%     % Cauchy dispersion parameters:
%     A1 = 1984.54821;
%     A2 = 104174632;
%     m = -0.00026*(Temperature-T_ref)+A0+A1/lambda^2+A2/lambda^4;

case 'DMSO'
        % our fit to:
        % Landolt-Börnstein III/47: Optical Constants
        % & Strecker, W., Spitaler, R.: Chem.Ber. 59 (1926) 1754
        x=lambda/1e3;
        A=1;
        B_IR=0.04419;
        C_IR=46390.67309;
        B_UV=1.09101;
        C_UV=12215.43949;
        
        m=sqrt(A+B_IR*x^2/(x^2-C_IR)+B_UV*x^2/(x^2-C_UV));
    
    
%     A0 = 1.4585;
%     % Cauchy dispersion parameters:
%     A1 = 8584;
%     A2 = -485600000;
%     m = -0.00044*(Temperature-T_ref)+A0+A1/lambda^2+A2/lambda^4;

case 'Benzonitryl' %C6H5CN
    A0 = 1.51158;
    % Cauchy dispersion parameters:
    A1 = 2600;
    A2 = 1.0919e9;
    m =-4.20518e-4*(Temperature-T_ref)+A0+A1/lambda^2+A2/lambda^4;   
    
case 'MetOH' % added by Tho Do Duc
        % Cauchy dispersion parameters:
        A0=1.3195;
        A1=3053.64419;
        A2=-34163639.3011;
        A3=2.622128e+12;
        m=-0.00040*(Temperature-T_ref)+A0+A1/(lambda^2)+A2/(lambda^4)+A3/(lambda^6);    

case 'Canola'
    % Rapeseed oil @22+/-6.5 deg.C
    % Cauchy dispersion parameters:
    A0=1.439298;
    A1=9339.331;
    A2=566636000;
    m = A0+A1/lambda^2+A2/lambda^4;

case 'uPS'
    % PS microspheres @24 deg. C
    % Cauchy dispersion parameters from Ma et al.:
    A0 = 1.5725;
    A1 = 3108;
    A2 = 347790000;
    m = A0+A1/lambda^2+A2/lambda^4;
    
case 'nanoTiO2'
    % probably amorphous
    % 24 deg. C, 35 nm diameter nanoparticles dispersed in water
    % M. N. Polyanskiy. Refractiveindex.info database of optical constants.
    % Sci. Data 11, 94 (2024) https://doi.org/10.1038/s41597-023-02898-2
    x=lambda/1e3;
    A = 1;
    B_IR = 4.6796;
    C_IR = 0.0401;
    m=sqrt(A+B_IR*x^2/(x^2-C_IR));
    
case 'H2O'
    % dispersion parameters for water from Harvey at al. , DOI: 10.1063/1.556029:
    a0 = 0.244257733;
    a1 = 9.74634476e-3;
    a2 = -3.73234996e-3;
    a3 = 2.68678472e-4;
    a4 = 1.5892057e-3;
    a5 = 2.45934259e-3;
    a6 = 0.90070492;
    a7 = -1.66626219e-2;
    lambda_uv = 0.229202;
    lambda_IR = 5.432937;
    laser = lambda/589;
    
    T_K = 273.15 + Temperature;
    
        if ((0<=Temperature)&&(Temperature<=40)) % this "if" concerns density of water
        % Thiesen formula for 0 - 40 degC from Tanaka
        b1 = -3.983035;
        b2 = 301.797;
        b3 = 522528.9;
        b4 = 69.34881;
        %     b5 = 999.974950; % Standard Mean Ocean Water (isotopic composition)
        b5 = 999.972; % tap water (Chappuis)
        Rho = ...
            b5 * ( 1 - ((Temperature + b1)^2)*(Temperature + b2)/(b3 * (Temperature + b4)));
        Rho = Rho/1000;
    else
        % from McCutcheon book
        Rho = 1 - (Temperature + 288.9414)/(508929.2*((Temperature + 68.12963))) ...
            *((Temperature - 3.9863)^2);
        %     Rho*1000
    end
    F = Rho * (a0 + a1*Rho + a2*T_K/273.15 + a3*laser^2*T_K/273.15 + a4/laser^2 ...
        + a5/(laser^2 - lambda_uv^2) + a6/(laser^2 - lambda_IR^2) + a7*Rho^2);
    m = sqrt((2*F+1)/(1-F));
    
    case 'Au'
            c = 299792458;
            h = 6.6260755e-34;
            e = 1.60217733e-19;
            omega = ((c*h)./(lambda.*e))*10^9;
            
            omega_p = 8.6;
            Gamma = 0.072;
            eps_zero = 1;
                        
            eps_up = 5.6;
            om_c = 2.4;
            S = 5.87;
            
            DeltaIm = eps_up./(1 + exp(S.*(om_c - omega)));
            
            AA = 7.78;
            Gamma_L = 0.977;
            omega_c = 2.4;
            
            DeltaRe = AA./( (omega - omega_c).^2 + Gamma_L);
                                                   
            eps = eps_zero - omega_p^2./(omega.^2 + 1i.*Gamma.*omega) + DeltaRe + 1i.*DeltaIm;
            m = sqrt(eps);
 
    case 'Ag'
        
        c = 299792458;
        h = 6.6260755e-34;
        e = 1.60217733e-19;
        omega = ((c*h)./(lambda.*e))*10^9;
            
        omega_p = 9;
        Gamma = 0.02;
        eps_zero = 1;
              
        om_c = 4.05;
        S = 7.9;
        A1 = 0.1;
        A2 = 6.27;
        A3 = 0.54;
        
        DeltaIm = A1 + (A2-A3.*omega)./(1 + exp(S.*(om_c - omega)));
        
        AA = 4.38;
        Gamma_L = 1;
        omega_c = 3.9;
        
        DeltaRe = AA./( (omega - omega_c).^2 + Gamma_L);
        
        eps = eps_zero - omega_p^2./(omega.^2 + 1i.*Gamma.*omega) + DeltaRe + 1i.*DeltaIm;
        
        m = sqrt(eps);
        
    otherwise
        m = 0;
           s = sprintf('Liquid not found \nm = %1.f !!!',m);
       errordlg(s,'m Error');     
end




