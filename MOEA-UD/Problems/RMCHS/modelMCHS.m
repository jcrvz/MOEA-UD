function [soo,moo] = modelMCHS(design_var, material, fluid, Q)
% soo = S_gen*T_a/Q;
% moo = (S_gen_$k*T_a/Q)_{$k in {ht,ff}}

    % Set default values
    if nargin < 4
        Q = 15; % Default value for Q
    end
    if nargin < 3
        fluid = 'H2O'; % Default fluid
    end
    if nargin < 2
        material = 'Al'; % Default material
    end
    
    %% Variables de diseño
    % #-
    alpha_c = design_var(1,1);
    beta 	= design_var(1,2);
    G_d 	= design_var(1,3);
    
    %% Parámetros variables
    
    % Altura del canal (m)
    H_c         = 1.7e-3;
    
    % Mitad del ancho de canal (m)
    w_c         = H_c*alpha_c/2;
    
    % Mitad del ancho de la aleta (m)
    w_p         = w_c/beta;
    
    %% Parámetros definidos
    
    % Conductividad térmica del sólido (W/m K)
    switch material
        case 'Al'  % Aluminio
            
            % Conductividad térmica del material (W/m K)
            k_m     = 237;
            
            % Densidad del material (kg/m3)
            rho_m   = 2707;
            
        case 'Cu'  % Cobre
            
            % Conductividad térmica del material (W/m K)
            k_m     = 401;
            
            % Densidad del material (kg/m3)
            rho_m   = 8954;
            
        case 'SiC' % Carburo de silicio
            
            % Conductividad térmica del material (W/m K)
            k_m     = 270;
            
            % Densidad del material (kg/m3)
            rho_m   = 3300;
            
        case 'AlN' % Nitruro de Aluminio
            
            % Conductividad térmica del material (W/m K)
            k_m     = 320;
            
            % Densidad del material (kg/m3)
            rho_m   = 3300;
            
        case 'HTCG'% HTCG
            
            % Conductividad térmica del material (W/m K)
            k_m     = 1900;
            
            % Densidad del material (kg/m3)
            rho_m   = 1000;
            
        case 'kb1k'% HTCG
            
            % Conductividad térmica del material (W/m K)
            k_m     = 1000;
            
            % Densidad del material (kg/m3)
            rho_m   = 1000;
            
        case 'kb2k'% HTCG
            
            % Conductividad térmica del material (W/m K)
            k_m     = 2000;
            
            % Densidad del material (kg/m3)
            rho_m   = 1000;
            
        case 'Si' % Silicio
            
            % Conductividad térmica del material (W/m K)
            k_m     = 148;
            
            % Densidad del material (kg/m3)
            rho_m   = 2330;
        otherwise
            error('Material no definido');
            
    end
    
    %% Propiedades del fluido
    switch fluid
        case 'Air'  % Aire
            
            % Conductividad térmica del fluido (W/m K)
            k_f     = 0.0261;
            
            % Densidad del fluido (aire) (kg/m3)
            rho_f   = 1.1614;
            
            % Viscosidad cinemática del fluido (aire)
            nu      = 1.58e-5;
            
            % Calor específico del aire a presión constante (J/kg K)
            c_p     = 1007;
            
        case 'Air+10HR'  % Aire con 10% de humedad
            
            % Conductividad térmica del fluido  (W/m K)
            k_f     = 0.02671;
            
            % Densidad del fluido  (kg/m3)
            rho_f   = 1.17434;
            
            % Viscosidad cinemática del fluido (aire)
            nu      = 1.80045e-5/rho_f;
            
            % Calor específico  a presión constante (J/kg K)
            c_p     = 1008.45;
            
        case 'Air+50HR'   % Aire con 50% de humedad
            
            % Conductividad térmica del fluido  (W/m K)
            k_f     = 0.027655;
            
            % Densidad del fluido  (kg/m3)
            rho_f   = 1.16808;
            
            % Viscosidad cinemática del fluido (aire)
            nu      = 1.78858e-5/rho_f;
            
            % Calor específico  a presión constante (J/kg K)
            c_p     = 1025.1;
            
        case 'Air+90HR'  % Aire con 90% de humedad
            
            % Conductividad térmica del fluido  (W/m K)
            k_f     = .028583;
            
            % Densidad del fluido  (kg/m3)
            rho_f   = 1.16182;
            
            % Viscosidad cinemática del fluido (aire)
            nu      = 1.77672e-5/rho_f;
            
            % Calor específico  a presión constante (J/kg K)
            c_p     = 1042.24;
            
        case 'NH3'         % Amoniaco
            
            % Conductividad térmica del fluido  (W/m K)
            k_f     = 0.027;
            
            % Densidad del fluido  (kg/m3)
            rho_f   = 0.7;
            
            % Viscosidad cinemática del fluido (aire)
            nu      = 1.4654e-5;
            
            % Calor específico  a presión constante (J/kg K)
            c_p     = 2158;
            
        case 'H2O'     % Agua
            
            % Conductividad térmica del fluido  (W/m K)
            k_f     = 0.625;
            
            % Densidad del fluido  (kg/m3)
            rho_f   = 994.2;
            
            % Viscosidad dinámica del fluido (aire)
            nu      = 7.25e-4/rho_f;
            
            % Calor específico  a presión constante (J/kg K)
            c_p     = 4178;
            
        case 'H2O+TiO2_01wt'    % Agua con nanopartículas de TiO2:
            
            % fb - fluido base, np - nanopartículas, nf - nanofluido
            % Conductividad térmica del fluido  (W/m K)
            k_fb    = 0.625;
            
            % Densidad del fluido  (kg/m3)
            rho_fb  = 994.2;
            
            % Viscosidad dinámica del fluido (aire)
            mu_fb   = 7.25e-4;
            
            % Calor específico  a presión constante (J/kg K)
            c_p_fb  = 4178;
            
            % Conductividad térmica de las nanopartículas  (W/m K)
            k_np    = 8.4;
            
            % Densidad de las nanopartículas (kg/m3)
            rho_np  = 4157;
            
            % Calor específico  a presión constante (J/kg K)
            c_p_np  = 710;
            
            % Fracción de volumen
            phi_nf  = 0.001; % 0.1%
            n       = 3;
            
            % Conductividad térmica del nanofluido  (W/m K)
            k_f = (k_np + (n - 1)*k_fb - (n - 1)*phi_nf*(k_fb - k_np))*...
                k_fb/(k_np + (n - 1)*k_fb + phi_nf*(k_fb - k_np));
            
            % Densidad del nanofluido  (kg/m3)
            rho_f   = (1 + phi_nf)*rho_fb + phi_nf*rho_np;
            
            % Viscosidad dinámica del nanofluido
            mu_nf   = mu_fb/(1 - phi_nf)^2.5;
            
            % Viscosidad cinemática del nanofluido
            nu      = mu_nf/rho_f;
            
            % Calor específico  a presión constante (J/kg K)
            c_p     = ((1 - phi_nf)*rho_fb*c_p_fb + ...
                phi_nf*rho_np*c_p_np)/rho_f;
            
        case 'H2O+TiO2_05wt'    % Agua con nanopartículas de TiO2:
            
            % fb - fluido base, np - nanopartículas, nf - nanofluido
            % Conductividad térmica del fluido  (W/m K)
            k_fb    = 0.625;
            
            % Densidad del fluido  (kg/m3)
            rho_fb  = 994.2;
            
            % Viscosidad dinámica del fluido (aire)
            mu_fb   = 7.25e-4;
            
            % Calor específico  a presión constante (J/kg K)
            c_p_fb  = 4178;
            
            % Conductividad térmica de las nanopartículas  (W/m K)
            k_np    = 8.4;
            
            % Densidad de las nanopartículas (kg/m3)
            rho_np  = 4157;
            
            % Calor específico  a presión constante (J/kg K)
            c_p_np  = 710;
            
            % Fracción de volumen
            phi_nf  = 0.005; % 0.5%
            n       = 3;
            
            % Conductividad térmica del nanofluido  (W/m K)
            k_f = (k_np + (n - 1)*k_fb - (n - 1)*phi_nf*(k_fb - k_np))*...
                k_fb/(k_np + (n - 1)*k_fb + phi_nf*(k_fb - k_np));
            
            % Densidad del nanofluido  (kg/m3)
            rho_f   = (1 + phi_nf)*rho_fb + phi_nf*rho_np;
            
            % Viscosidad dinámica del nanofluido
            mu_nf   = mu_fb/(1 - phi_nf)^2.5;
            
            % Viscosidad cinemática del nanofluido
            nu      = mu_nf/rho_f;
            
            % Calor específico  a presión constante (J/kg K)
            c_p     = ((1 - phi_nf)*rho_fb*c_p_fb + ...
                phi_nf*rho_np*c_p_np)/rho_f;
            
        case 'H2O+TiO2_09wt'    % Agua con nanopartículas de TiO2:
            
            % fb - fluido base, np - nanopartículas, nf - nanofluido
            % Conductividad térmica del fluido  (W/m K)
            k_fb    = 0.625;
            
            % Densidad del fluido  (kg/m3)
            rho_fb  = 994.2;
            
            % Viscosidad dinámica del fluido (aire)
            mu_fb   = 7.25e-4;
            
            % Calor específico  a presión constante (J/kg K)
            c_p_fb  = 4178;
            
            % Conductividad térmica de las nanopartículas  (W/m K)
            k_np    = 8.4;
            
            % Densidad de las nanopartículas (kg/m3)
            rho_np  = 4157;
            
            % Calor específico  a presión constante (J/kg K)
            c_p_np  = 710;
            
            % Fracción de volumen
            phi_nf  = 0.009; % 0.9%
            n       = 3;
            
            % Conductividad térmica del nanofluido  (W/m K)
            k_f = (k_np + (n - 1)*k_fb - (n - 1)*phi_nf*(k_fb - k_np))*...
                k_fb/(k_np + (n - 1)*k_fb + phi_nf*(k_fb - k_np));
            
            % Densidad del nanofluido  (kg/m3)
            rho_f   = (1 + phi_nf)*rho_fb + phi_nf*rho_np;
            
            % Viscosidad dinámica del nanofluido
            mu_nf   = mu_fb/(1 - phi_nf)^2.5;
            
            % Viscosidad cinemática del nanofluido
            nu      = mu_nf/rho_f;
            
            % Calor específico  a presión constante (J/kg K)
            c_p     = ((1 - phi_nf)*rho_fb*c_p_fb + ...
                phi_nf*rho_np*c_p_np)/rho_f;
            
        case 'H2O+TiO2_10wt'    % Agua con nanopartículas de TiO2:
            
            % fb - fluido base, np - nanopartículas, nf - nanofluido
            % Conductividad térmica del fluido  (W/m K)
            k_fb    = 0.625;
            
            % Densidad del fluido  (kg/m3)
            rho_fb  = 994.2;
            
            % Viscosidad dinámica del fluido (aire)
            mu_fb   = 7.25e-4;
            
            % Calor específico  a presión constante (J/kg K)
            c_p_fb  = 4178;
            
            % Conductividad térmica de las nanopartículas  (W/m K)
            k_np    = 8.4;
            
            % Densidad de las nanopartículas (kg/m3)
            rho_np  = 4157;
            
            % Calor específico  a presión constante (J/kg K)
            c_p_np  = 710;
            
            % Fracción de volumen
            phi_nf  = 0.01; % 1.0%
            n       = 3;
            
            % Conductividad térmica del nanofluido  (W/m K)
            k_f = (k_np + (n - 1)*k_fb - (n - 1)*phi_nf*(k_fb - k_np))*...
                k_fb/(k_np + (n - 1)*k_fb + phi_nf*(k_fb - k_np));
            
            % Densidad del nanofluido  (kg/m3)
            rho_f   = (1 + phi_nf)*rho_fb + phi_nf*rho_np;
            
            % Viscosidad dinámica del nanofluido
            mu_nf   = mu_fb/(1 - phi_nf)^2.5;
            
            % Viscosidad cinemática del nanofluido
            nu      = mu_nf/rho_f;
            
            % Calor específico  a presión constante (J/kg K)
            c_p     = ((1 - phi_nf)*rho_fb*c_p_fb + ...
                phi_nf*rho_np*c_p_np)/rho_f;
            
        case 'H2O+TiO2_50wt'    % Agua con nanopartículas de TiO2:
            
            % fb - fluido base, np - nanopartículas, nf - nanofluido
            % Conductividad térmica del fluido  (W/m K)
            k_fb    = 0.625;
            
            % Densidad del fluido  (kg/m3)
            rho_fb  = 994.2;
            
            % Viscosidad dinámica del fluido (aire)
            mu_fb   = 7.25e-4;
            
            % Calor específico  a presión constante (J/kg K)
            c_p_fb  = 4178;
            
            % Conductividad térmica de las nanopartículas  (W/m K)
            k_np    = 8.4;
            
            % Densidad de las nanopartículas (kg/m3)
            rho_np  = 4157;
            
            % Calor específico  a presión constante (J/kg K)
            c_p_np  = 710;
            
            % Fracción de volumen
            phi_nf  = 0.05; % 5%
            n       = 3;
            
            % Conductividad térmica del nanofluido  (W/m K)
            k_f = (k_np + (n - 1)*k_fb - (n - 1)*phi_nf*(k_fb - k_np))*...
                k_fb/(k_np + (n - 1)*k_fb + phi_nf*(k_fb - k_np));
            
            % Densidad del nanofluido  (kg/m3)
            rho_f   = (1 + phi_nf)*rho_fb + phi_nf*rho_np;
            
            % Viscosidad dinámica del nanofluido
            mu_nf   = mu_fb/(1 - phi_nf)^2.5;
            
            % Viscosidad cinemática del nanofluido
            nu      = mu_nf/rho_f;
            
            % Calor específico  a presión constante (J/kg K)
            c_p     = ((1 - phi_nf)*rho_fb*c_p_fb + ...
                phi_nf*rho_np*c_p_np)/rho_f;
            
        case 'H2O+TiO2_90wt'    % Agua con nanopartículas de TiO2:
            
            % fb - fluido base, np - nanopartículas, nf - nanofluido
            % Conductividad térmica del fluido  (W/m K)
            k_fb    = 0.625;
            
            % Densidad del fluido  (kg/m3)
            rho_fb  = 994.2;
            
            % Viscosidad dinámica del fluido (aire)
            mu_fb   = 7.25e-4;
            
            % Calor específico  a presión constante (J/kg K)
            c_p_fb  = 4178;
            
            % Conductividad térmica de las nanopartículas  (W/m K)
            k_np    = 8.4;
            
            % Densidad de las nanopartículas (kg/m3)
            rho_np  = 4157;
            
            % Calor específico  a presión constante (J/kg K)
            c_p_np  = 710;
            
            % Fracción de volumen
            phi_nf  = 0.09; % 9%
            n       = 3;
            
            % Conductividad térmica del nanofluido  (W/m K)
            k_f = (k_np + (n - 1)*k_fb - (n - 1)*phi_nf*(k_fb - k_np))*...
                k_fb/(k_np + (n - 1)*k_fb + phi_nf*(k_fb - k_np));
            
            % Densidad del nanofluido  (kg/m3)
            rho_f   = (1 + phi_nf)*rho_fb + phi_nf*rho_np;
            
            % Viscosidad dinámica del nanofluido
            mu_nf   = mu_fb/(1 - phi_nf)^2.5;
            
            % Viscosidad cinemática del nanofluido
            nu      = mu_nf/rho_f;
            
            % Calor específico  a presión constante (J/kg K)
            c_p     = ((1 - phi_nf)*rho_fb*c_p_fb + ...
                phi_nf*rho_np*c_p_np)/rho_f;
            
        case 'H2O+Al_01wt'    % Agua con nanopartículas de TiO2:
            
            % fb - fluido base, np - nanopartículas, nf - nanofluido
            % Conductividad térmica del fluido  (W/m K)
            k_fb    = 0.625;
            
            % Densidad del fluido  (kg/m3)
            rho_fb  = 994.2;
            
            % Viscosidad dinámica del fluido (aire)
            mu_fb   = 7.25e-4;
            
            % Calor específico  a presión constante (J/kg K)
            c_p_fb  = 4178;
            
            % Conductividad térmica de las nanopartículas  (W/m K)
            k_np    = 36;
            
            % Densidad de las nanopartículas (kg/m3)
            rho_np  = 3380;
            
            % Calor específico  a presión constante (J/kg K)
            c_p_np  = 773;
            
            % Fracción de volumen
            phi_nf  = 0.001; % 0.1%
            n       = 3;
            
            % Conductividad térmica del nanofluido  (W/m K)
            k_f = (k_np + (n - 1)*k_fb - (n - 1)*phi_nf*(k_fb - k_np))*...
                k_fb/(k_np + (n - 1)*k_fb + phi_nf*(k_fb - k_np));
            
            % Densidad del nanofluido  (kg/m3)
            rho_f   = (1 + phi_nf)*rho_fb + phi_nf*rho_np;
            
            % Viscosidad dinámica del nanofluido
            mu_nf   = mu_fb/(1 - phi_nf)^2.5;
            
            % Viscosidad cinemática del nanofluido
            nu      = mu_nf/rho_f;
            
            % Calor específico  a presión constante (J/kg K)
            c_p     = ((1 - phi_nf)*rho_fb*c_p_fb + ...
                phi_nf*rho_np*c_p_np)/rho_f;
            
        case 'H2O+Al_05wt'    % Agua con nanopartículas de TiO2:
            
            % fb - fluido base, np - nanopartículas, nf - nanofluido
            % Conductividad térmica del fluido  (W/m K)
            k_fb    = 0.625;
            
            % Densidad del fluido  (kg/m3)
            rho_fb  = 994.2;
            
            % Viscosidad dinámica del fluido (aire)
            mu_fb   = 7.25e-4;
            
            % Calor específico  a presión constante (J/kg K)
            c_p_fb  = 4178;
            
            % Conductividad térmica de las nanopartículas  (W/m K)
            k_np    = 36;
            
            % Densidad de las nanopartículas (kg/m3)
            rho_np  = 3380;
            
            % Calor específico  a presión constante (J/kg K)
            c_p_np  = 773;
            
            % Fracción de volumen
            phi_nf  = 0.005; % 0.5%
            n       = 3;
            
            % Conductividad térmica del nanofluido  (W/m K)
            k_f = (k_np + (n - 1)*k_fb - (n - 1)*phi_nf*(k_fb - k_np))*...
                k_fb/(k_np + (n - 1)*k_fb + phi_nf*(k_fb - k_np));
            
            % Densidad del nanofluido  (kg/m3)
            rho_f   = (1 + phi_nf)*rho_fb + phi_nf*rho_np;
            
            % Viscosidad dinámica del nanofluido
            mu_nf   = mu_fb/(1 - phi_nf)^2.5;
            
            % Viscosidad cinemática del nanofluido
            nu      = mu_nf/rho_f;
            
            % Calor específico  a presión constante (J/kg K)
            c_p     = ((1 - phi_nf)*rho_fb*c_p_fb + ...
                phi_nf*rho_np*c_p_np)/rho_f;
            
        case 'H2O+Al_09wt'    % Agua con nanopartículas de TiO2:
            
            % fb - fluido base, np - nanopartículas, nf - nanofluido
            % Conductividad térmica del fluido  (W/m K)
            k_fb    = 0.625;
            
            % Densidad del fluido  (kg/m3)
            rho_fb  = 994.2;
            
            % Viscosidad dinámica del fluido (aire)
            mu_fb   = 7.25e-4;
            
            % Calor específico  a presión constante (J/kg K)
            c_p_fb  = 4178;
            
            % Conductividad térmica de las nanopartículas  (W/m K)
            k_np    = 36;
            
            % Densidad de las nanopartículas (kg/m3)
            rho_np  = 3380;
            
            % Calor específico  a presión constante (J/kg K)
            c_p_np  = 773;
            
            % Fracción de volumen
            phi_nf  = 0.009; % 0.9%
            n       = 3;
            
            % Conductividad térmica del nanofluido  (W/m K)
            k_f = (k_np + (n - 1)*k_fb - (n - 1)*phi_nf*(k_fb - k_np))*...
                k_fb/(k_np + (n - 1)*k_fb + phi_nf*(k_fb - k_np));
            
            % Densidad del nanofluido  (kg/m3)
            rho_f   = (1 + phi_nf)*rho_fb + phi_nf*rho_np;
            
            % Viscosidad dinámica del nanofluido
            mu_nf   = mu_fb/(1 - phi_nf)^2.5;
            
            % Viscosidad cinemática del nanofluido
            nu      = mu_nf/rho_f;
            
            % Calor específico  a presión constante (J/kg K)
            c_p     = ((1 - phi_nf)*rho_fb*c_p_fb + ...
                phi_nf*rho_np*c_p_np)/rho_f;
            
        case 'H2O+Al_10wt'    % Agua con nanopartículas de TiO2:
            
            % fb - fluido base, np - nanopartículas, nf - nanofluido
            % Conductividad térmica del fluido  (W/m K)
            k_fb    = 0.625;
            
            % Densidad del fluido  (kg/m3)
            rho_fb  = 994.2;
            
            % Viscosidad dinámica del fluido (aire)
            mu_fb   = 7.25e-4;
            
            % Calor específico  a presión constante (J/kg K)
            c_p_fb  = 4178;
            
            % Conductividad térmica de las nanopartículas  (W/m K)
            k_np    = 36;
            
            % Densidad de las nanopartículas (kg/m3)
            rho_np  = 3380;
            
            % Calor específico  a presión constante (J/kg K)
            c_p_np  = 773;
            
            % Fracción de volumen
            phi_nf  = 0.01; % 1%
            n       = 3;
            
            % Conductividad térmica del nanofluido  (W/m K)
            k_f = (k_np + (n - 1)*k_fb - (n - 1)*phi_nf*(k_fb - k_np))*...
                k_fb/(k_np + (n - 1)*k_fb + phi_nf*(k_fb - k_np));
            
            % Densidad del nanofluido  (kg/m3)
            rho_f   = (1 + phi_nf)*rho_fb + phi_nf*rho_np;
            
            % Viscosidad dinámica del nanofluido
            mu_nf   = mu_fb/(1 - phi_nf)^2.5;
            
            % Viscosidad cinemática del nanofluido
            nu      = mu_nf/rho_f;
            
            % Calor específico  a presión constante (J/kg K)
            c_p     = ((1 - phi_nf)*rho_fb*c_p_fb + ...
                phi_nf*rho_np*c_p_np)/rho_f;
            
        case 'H2O+Al_50wt'    % Agua con nanopartículas de TiO2:
            
            % fb - fluido base, np - nanopartículas, nf - nanofluido
            % Conductividad térmica del fluido  (W/m K)
            k_fb    = 0.625;
            
            % Densidad del fluido  (kg/m3)
            rho_fb  = 994.2;
            
            % Viscosidad dinámica del fluido (aire)
            mu_fb   = 7.25e-4;
            
            % Calor específico  a presión constante (J/kg K)
            c_p_fb  = 4178;
            
            % Conductividad térmica de las nanopartículas  (W/m K)
            k_np    = 36;
            
            % Densidad de las nanopartículas (kg/m3)
            rho_np  = 3380;
            
            % Calor específico  a presión constante (J/kg K)
            c_p_np  = 773;
            
            % Fracción de volumen
            phi_nf  = 0.05; % 5%
            n       = 3;
            
            % Conductividad térmica del nanofluido  (W/m K)
            k_f = (k_np + (n - 1)*k_fb - (n - 1)*phi_nf*(k_fb - k_np))*...
                k_fb/(k_np + (n - 1)*k_fb + phi_nf*(k_fb - k_np));
            
            % Densidad del nanofluido  (kg/m3)
            rho_f   = (1 + phi_nf)*rho_fb + phi_nf*rho_np;
            
            % Viscosidad dinámica del nanofluido
            mu_nf   = mu_fb/(1 - phi_nf)^2.5;
            
            % Viscosidad cinemática del nanofluido
            nu      = mu_nf/rho_f;
            
            % Calor específico  a presión constante (J/kg K)
            c_p     = ((1 - phi_nf)*rho_fb*c_p_fb + ...
                phi_nf*rho_np*c_p_np)/rho_f;
            
        case 'H2O+Al_90wt'    % Agua con nanopartículas de TiO2:
            
            % fb - fluido base, np - nanopartículas, nf - nanofluido
            % Conductividad térmica del fluido  (W/m K)
            k_fb    = 0.625;
            
            % Densidad del fluido  (kg/m3)
            rho_fb  = 994.2;
            
            % Viscosidad dinámica del fluido (aire)
            mu_fb   = 7.25e-4;
            
            % Calor específico  a presión constante (J/kg K)
            c_p_fb  = 4178;
            
            % Conductividad térmica de las nanopartículas  (W/m K)
            k_np    = 36;
            
            % Densidad de las nanopartículas (kg/m3)
            rho_np  = 3380;
            
            % Calor específico  a presión constante (J/kg K)
            c_p_np  = 773;
            
            % Fracción de volumen
            phi_nf  = 0.09; % 9%
            n       = 3;
            
            % Conductividad térmica del nanofluido  (W/m K)
            k_f = (k_np + (n - 1)*k_fb - (n - 1)*phi_nf*(k_fb - k_np))*...
                k_fb/(k_np + (n - 1)*k_fb + phi_nf*(k_fb - k_np));
            
            % Densidad del nanofluido  (kg/m3)
            rho_f   = (1 + phi_nf)*rho_fb + phi_nf*rho_np;
            
            % Viscosidad dinámica del nanofluido
            mu_nf   = mu_fb/(1 - phi_nf)^2.5;
            
            % Viscosidad cinemática del nanofluido
            nu      = mu_nf/rho_f;
            
            % Calor específico  a presión constante (J/kg K)
            c_p     = ((1 - phi_nf)*rho_fb*c_p_fb + ...
                phi_nf*rho_np*c_p_np)/rho_f;
        otherwise
            error('Fluido no definido');
            
    end
    
    % Número de Prandtl
    Pr          = nu*rho_f*c_p/k_f;
    
    % Espesor de la base
    H_b         = 0.1e-3;
    
    % Ancho del disipador de calor (m)
    W_d         = 51e-3;
    
    % Longitud del canal en dirección del flujo (m)
    L_d         = 51e-3;
    
    % Ancho del chip (m)
    W_i         = W_d;
    
    % Longitud del chip (m)
    L_i         = L_d;
    
    % Densidad de flujo de calor en W/m2
    %q           = 5e4;
    
    % Temperatura absoluta del ambiente (K)
    T_a         = 300;
    
    % Por conducción en el contacto interfaz por área
    R_iA        = 2.75e-4;
    
    %% Parámetros de geometría
    
    % Relación de aspecto del canal
    % alpha_c     = 2*w_c/H_c;
    
    % Relación de espacio de la aleta
    % beta        = w_c/w_p;
    
    % Número de microcanales
    N           = floor((W_d/2 - w_p)/(w_c + w_p));
    
    % Diámetro hidráulico
    A_c_        = 2*w_c*H_c;
    P_          = 2*(H_c + 2*w_c);
    D_h         = 4*A_c_/P_;
    
    %% Parámetros de flujo de calor y de fluido
    
    % Flujo másico del fluido por canal
    mdot        = rho_f*G_d/(2*N);
    
    % Velocidad promedio del flujo (m/s)
    U_avg       = mdot/(rho_f*w_c*H_c);
    
    % Número de Reynolds basado en el diámetro hidráulico
    Re_Dh       = D_h*U_avg/nu;
    
    % Área de la superficie inferior del disipador (m2)
    A_b         = W_d*L_d;
    
    % Área del chip
    A_i         = W_i*L_i;
    
    % Flujo de transferencia de calor total (W)
    q           = Q/A_i;
    
    %% Número de Nusselt y factor de fricción según el régimen
    
    % if Re_Dh < 2300
        % Flujo laminar completamente desarrollado
        
        % Número de Nusselt para flujo laminar
        Nu_Dh       = 2.253 + 8.164*(1/(alpha_c + 1))^(1.5);
        
        % Factor de fricción para flujo laminar
        f           = sqrt((3.2*(Re_Dh*D_h/L_d)^(0.57))^2 + ...
            (4.7 + 19.64*(alpha_c^2 + 1)/(alpha_c + 1)^2)^2)/Re_Dh;
        
    % else
    %     % Flujo turbulento completamente desarrollado
    %     
    %     % Factor de fricción para flujo turbulento
    %     f               = 1/(0.79*log(Re_Dh) - 1.64)^2;
    %     
    %     % Número de Nusselt para flujo turbulento
    %     Nu_Dh           = (f/2)*(Re_Dh - 1000)*Pr/...
    %         (1 + 12.7*(f/2)^.5*(Pr^(2/3) - 1));
    % end
    
    %% Resistencia térmica del disipador de calor
    
    % Resistencia de la pasta
    R_i          = R_iA/A_i;
    
    %       Coeficiente de película (W/K m2)
    h_avg       = k_f*Nu_Dh/D_h;
    
    %       Parámetro de la aleta con altura del canal
    mHc         = sqrt(2*(2*w_p + L_d)*h_avg/(k_m*2*w_p*L_d))*H_c;
    
    %       Eficiencia de la aleta
    eta_p       = tanh(mHc)/mHc;
    
    %       Área efectiva para la transferencia de calor
    A_eff       = 2*N*(eta_p*H_c + w_c)*L_d;
    
    % Resistencia por convección
    R_conv      = 1/(h_avg*A_eff);
    
    % Resistencia por conducción en el fluido
    R_f         = 1/(rho_f*G_d*c_p);
    
    % Resistencia por constricción
    R_const     = ((1 + beta)/(pi*k_m*A_b*beta))*...
        log(1/sin(0.5*pi/(1 + beta)))*H_c*alpha_c;
    
    %       Parámetros de geometría
    a           = sqrt(A_i/pi);
    b           = sqrt(A_b/pi);
    
    %       Relaciones de geometría
    tau         = H_b/b;
    epsilon     = a/b;
    
    % Resistencia equivalente entre la base y el fluido
    R_o         = R_conv + R_f;  %+ R_const;
    
    %       Número de Biot
    Bi          = 1/(pi*k_m*b*R_o);
    
    %       Parámetros adimensionales
    lambda_c    = pi + 1/(sqrt(pi)*epsilon);
    Phi_c       = (tanh(lambda_c*tau) + lambda_c/Bi)/...
        (1 + lambda_c*tanh(lambda_c*tau)/Bi);
    Psi_avg     = epsilon*tau/sqrt(pi) + 0.5*Phi_c*(1 - epsilon)^(3/2);
    
    % Resistencia por dispersion en la base
    R_b         = Psi_avg/(sqrt(pi)*k_m*a);
    
    % Resistencia térmica equivalente del disipador de calor
    R_eq        = R_o; % + R_b + R_i;
    
    %% Temperatura en la interfaz
    T_i         = T_a + R_eq*Q;
    
    %% Caída de presión
    
    % Suma de pérdidas en la entrada y en la salida
    k_ce        = 1.79 - 2.32*(beta/(1 + beta)) + 0.53*(beta/(1 + beta))^2;
    
    % Caída de presión en el disipador
    DeltaP      = 0.5*rho_f*U_avg^2*(f*(L_d/D_h) + k_ce);
    
    % Longitud del tubo
    L_tu        = 0.5;
    
    % Diámetro del tubo
    D_tu        = 1.9e-2;
    
    % Velocidad promedio en los tubos de alimentación
    U_avgtub    = 4*G_d/(pi*D_tu^2);
    
    % Número de Reynolds en el tubo
    Re_Dtu      = U_avgtub*D_tu/nu;
    
    % Factor de fricción del tubo
    f_tu        = 4*(0.09290 + 1.01612/(L_tu/D_tu))*...
        Re_Dtu^(-0.268-0.3193/(L_tu/D_tu));
    
    % Relación de áreas de la sección transversal del tubo y del disipador
    A_tuhs      = 0.25*pi*D_tu^2/(W_d*H_c);
    
    % Caída de presión en los tubos
    DeltaP_tu   = 0.5*rho_f*U_avgtub^2*(0.42 + (1 - A_tuhs^2)^2 + ...
        0.42*(1 - A_tuhs^2) + 1 + 2*f_tu*L_tu/D_tu);
    
    % Caída de presión total
    DeltaP_total = DeltaP; %+ DeltaP_tu;
    
    % Potencia de bombeo
    Phi         = G_d*DeltaP_total;
    
    %% Tasas de generación de entropía
    
    % Tasa de generación de entropía por transferencia de calor
    S_genht     = Q^2*R_eq/(T_a*T_i);
    
    % Tasa de generación de entropía por fricción del fluido
    S_genff     = G_d*DeltaP_total/T_a;
    
    % Tasa de generación de entropía total
    S_gen       = S_genht + S_genff;
    
    %% Masa del disipador
    Md          = rho_m*L_d*(2*W_d*H_b + 2*w_p*H_c*(N + 1));
    
    %% Variables de salida
    soo         = S_gen*T_a/Q;
    moo         = [S_genht*T_a/Q, S_genff*T_a/Q];
    
    % Los parámetros adicionales de salida
    outpar = [S_gen,S_genht,S_genff,U_avg,Phi,Md,h_avg,Re_Dh,R_i,R_b,R_const,R_conv,R_f,N,T_i];
end
