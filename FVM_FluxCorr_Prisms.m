%% Metainformation:
%
% Author:           Philipp Selzer
% Version:          2.0
% Date:             20.11.2019
% Last update:      24.03.2022
% Purpose:          This code is the FVM-flux-correction computing conforming
%                   velocity fields for consistent particle tracking.
%                   It has to be run once per model.

% Correspondance:   philipp.selzer@gmx.net
% Copyright(c) 2019: Philipp Selzer
% License: An extended version of the LGPLv3 (http://www.gnu.org/licenses/lgpl-3.0.en.html), see
% attached license-file

% Acknowledgement: The author thanks Jonas Allgeier and Daniel Erdal for
% their valuable comments and suggestions, as well as for testing my codes
% and their help in the debugging process. The author thanks his PhD
% supervisors Olaf A. Cirpka and René Therrien for their valuable 
% contributions and their advise.


%%
clear all
close all
clc

tic
format long

prefix = 'floodplain_model'; % TO ADAPT
filename_mprops = 'floodplain_calibration.mprops';

%% ===================================================================
% Read all relevant HGS output files
%% ===================================================================

% read the binary node ids of the elements
n_dim = 3; % number of spatial dimensions, should stay 3
no_nodes_per_element = 6; % This code only considers triangular prisms, this
% value should stay 6

no_elem =138565;% TO ADAPT
no_layers = 35;% TO ADAPT
no_nodes = 74412;% TO ADAPT

no_top_nodes = no_nodes/(no_layers + 1);

dim = 1; % for reading scalar values

plotting = true;

%% TO DO: Loop over all relevant time_steps

transient =  true; % if transient == false FVM_corr_target_time_IDs should be FVM_corr_target_time_IDs=1;
% and target_times can be any value, e.g. "Inf"
surface_domain = false; % This flag should stay, at best, false. The consideration of a surface
% domain within the FVM flux reconstruction still needs some improvement
% for mimicking the primal solution accurately. Results in the current
% state which look suboptimally may occur.

%%  read grid and material properties

elem_per_layer = no_elem/no_layers;
elem_top_layer = elem_per_layer;

% read the binary elements
filename1 = strcat(prefix,'o.elements_pm');
[node_ids] = read_HGS_binary_elements(no_nodes_per_element,no_elem,filename1);

% read the binary coordinates
filename2 = strcat(prefix,'o.coordinates_pm');
[coordinates] = read_HGS_binary_coordinates(n_dim,no_nodes,filename2);


K_from_file = true;

if K_from_file
    
    filename4 = strcat(prefix,'o.ElemK_pm.0001');
    
    xyz_cr = 4; % 3 dimensions + carriage return characters due to an individual
    % write statement per element in Fortran 95
    
    [K_sat] = read_HGS_binary_Ksat(xyz_cr,no_elem,filename4);
    
    k_xx = K_sat(:,1);
    k_yy = K_sat(:,2);
    k_zz = K_sat(:,3);
    
end

%% solver settings (the solver is the bicgstab with an ilu-preconditioner)

tol = 1e-14;
maxit = 1000;

%% FLAGS ALTERING THE FLOW SOLUTION IN THE UNSATURATED ZONE

unsaturated_flow = true; % TO ADAPT LINE, True: if a part of the domain is unsaturated, false: fully saturated flow everywhere

% Weighting scheme of relative permeability on the faces, only important if unsaturated_flow == true
lincentro_upface = true; %Should always stay true, alternative schemes were deleted for this version, but are, of course, possible

% % % enforce_unsatflow_correctness = false; % This flag is very efficient, but numerically brutal...
% If you consider this: THE PARAMETERS OF THE HGS SOLUTION ARE NOT FULLY CONSISTENT.
% YOUR HGS SOLUTION MAY BE FLAWED. CHECK THE NUMERICAL PARAMETERS AND THE
% GRID (REFINEMENT IN THE HORIZONTAL AND SHAPE OF TRIANGLES).
% ESPECIALLY HAVE A LOOK AT YOUR VAN GENUCHTEN/BROOKS-COREY PARAMETERS. IF
% POSSIBLE YOU MAY CONSIDER A SMOOTH TRANSITION AT HARSH MATERIAL
% INTERFACES.

% % % if enforce_unsatflow_correctness
% % %     disp('You use the flag "enforce_unsatflow_correctness": ')
% % %     disp('Have you considered refining your grid in the horizontal?')
% % %     disp('Are the problematic regions discretized by nearly equilateral triangles?')
% % %     disp('Are your convergence criteria of the underlying FEM solution strict enough?')
% % %     disp('Note: If you apply this flag, the most upper layer will not give physically realistic results for the flow over the top of the domain!')
% % % end

%% Zones
zones_from_file = true; % IF possible, set this to true as well
% if true - zones are matched from hydraulic conductivity, otherwise zones are read from HGS output-file

if zones_from_file
    zones = importdata('zones.dat');
    zone_for_element = zones(:,2);
    unique_zones = unique(zone_for_element);
    no_of_zones = max(unique_zones);
else
    no_of_zones = 12;% IMPORTANT: PLUG IN THE NUMBER OF ZONES FOR SPECIFIC STORAGE AND/OR
    % RECONSTRUCTION OF UNSATURATED PARAMETERS: PLUG THIS IN, IF
    % zones_from_file == false
    zone_for_element = 1*ones(no_elem,1); % Cange this for manual zone allocation accordingly!!!
end


%% Porosity
porosity_from_bin_file = true; % set this also true, if possible

if porosity_from_bin_file
    filenameP = strcat(prefix,'o.ElemPor_pm.0001');
    [n] = read_HGS_binary_porosity(dim,no_elem,filenameP);
end

% Manual list of unsaturated parameters: ("nan" translates to "always
% saturated")

theta_from_file = false;

if ~theta_from_file
    Sw_r_s = [0.15 0.15 0.15 0.15  0.15 0.15 nan nan 0.17 0.18 0.25 0.25];
    alpha_s = [1.235001 1.235001 1.235001 1.235001 1.235001 1.235001 nan nan 0.012178 11.512762 0.513486 0.513486];
    N_s = [3.805919 3.805919 3.805919 3.805919 3.805919 3.805919 nan nan 1.408750 2.195659 1.798006 1.798006];
end

if transient
    S_sp_s = ones(1,12)*1e-5; % Change this, if needed!
end

if ~porosity_from_bin_file
    n_v = [0.4 0.4 0.4 0.4 0.4 0.4 0.4];
    n = nan(no_elem,1);
end


if ~K_from_file
    k_xx_s = [0.000005 0.00018 1e-7 0.00018 0.000003 0.000001 0.00003];
    k_yy_s = [0.000005 0.00018 1e-7 0.00018 0.000003 0.000001 0.00003];
    k_zz_s = [0.0000005 0.000018 1e-7 0.000018 0.000003 0.0000001 0.000003];
    
    k_xx = nan(no_elem,1);
    k_yy = nan(no_elem,1);
    k_zz = nan(no_elem,1);
    
end

if ~theta_from_file
    alpha = nan(no_elem,1);
    N = nan(no_elem,1);
    Sw_r = nan(no_elem,1);
end

if transient
    S_sp = nan(no_elem,1);
end


for z = 1:no_of_zones
    
    if ~porosity_from_bin_file
        n(zone_for_element==z) = n_v(z);
    end
    
    if ~theta_from_file
        if ~isnan(Sw_r_s(z))
            alpha(zone_for_element==z) = alpha_s(z);
            N(zone_for_element==z) = N_s(z);
            Sw_r(zone_for_element==z) = Sw_r_s(z);
        else
            alpha(zone_for_element==z) = 1;
            N(zone_for_element==z) = 1;
            Sw_r(zone_for_element==z) = n(zone_for_element==z);
        end
    end
    
    if transient
        S_sp(zone_for_element==z) = S_sp_s(z);
    end
    
    
    if ~K_from_file
        k_xx(zone_for_element==z) = k_xx_s(z);
        k_yy(zone_for_element==z) = k_yy_s(z);
        k_zz(zone_for_element==z) = k_zz_s(z);
    end
    
end

%% load grid geometry

corr_fact_orth_name = strcat('corr_fact_orth_save_',prefix,'.mat');
A_f_all_name = strcat('A_f_all_',prefix,'.mat');
l_elem_all_name = strcat('l_elem_all_',prefix,'.mat');
l_neighb_all_name = strcat('l_neighb_all_',prefix,'.mat');
A_horiz_name = strcat('A_horiz_',prefix,'.mat');
A_top_name = strcat('A_top_',prefix,'.mat');
all_neighb_indices_name = strcat('all_neighb_indices_',prefix,'.mat');
hQ_all_name = strcat('hQ_all_',prefix,'.mat');
face_index_mat_name = strcat('face_index_mat_',prefix,'.mat');
Dirichlet_elements_mult_name = strcat('Dirichlet_elements_mult_',prefix,'.mat');
face_index_mat_multD_name = strcat('face_index_mat_multD_',prefix,'.mat');
Neumann_elements_mult_name = strcat('Neumann_elements_mult_',prefix,'.mat');
face_index_mat_multN_name = strcat('face_index_mat_mult_',prefix,'.mat');
face_flag_D_N_name = strcat('face_flag_D_N_',prefix,'.mat');
Vol_prisms_name = strcat('Vol_prisms_',prefix,'.mat');

corr_fact_orth_save = importdata(corr_fact_orth_name);
A_f_all = importdata(A_f_all_name);
l_elem_all = importdata(l_elem_all_name);
l_neighb_all = importdata(l_neighb_all_name);
face_index_mat = importdata(face_index_mat_name);
vector_of_d_bounds_mult = importdata(Dirichlet_elements_mult_name);
face_index_mat_multD = importdata(face_index_mat_multD_name);
vector_of_n_bounds_mult = importdata(Neumann_elements_mult_name);
face_index_mat_multN = importdata(face_index_mat_multN_name);
A_horiz = importdata(A_horiz_name);
A_top = importdata(A_top_name);
all_neighb_indices = importdata(all_neighb_indices_name);
face_flag_D_N = importdata(face_flag_D_N_name);
Vol_prisms = importdata(Vol_prisms_name);

if plotting
    centro_of_elem_name = strcat('centro_of_elem_',prefix,'.mat');
    centro_of_elem = importdata(centro_of_elem_name);
end


% Transient settings

FVM_corr_target_time_IDs = [2,3,4,5];% steady state = [1] % Minimum should be 2, if transient == true/ for steady-state simulations this should stay "1"

if transient
    target_times = [3e7,3e8,3e9,3e11,3e12]; % steady state: 1 target time
end


for target_time_step = FVM_corr_target_time_IDs
    
    % automatize strings
    
    if target_time_step<=9
        target_string_no = '000';
        preceding_string_no = '000';
    elseif target_time_step==10
        target_string_no = '00';
        preceding_string_no = '000';
    elseif target_time_step>= 11 && target_time_step <= 99
        target_string_no = '00';
        preceding_string_no = '00';
    elseif target_time_step == 100
        target_string_no = '0';
        preceding_string_no = '00';
    elseif target_time_step>=101 && target_time_step <= 999
        target_string_no = '0';
        preceding_string_no = '0';
    elseif target_time_step == 1000
        target_string_no = '';
        preceding_string_no = '0';
    elseif target_time_step >= 1001 && target_time_step <= 9999
        target_string_no = '';
        preceding_string_no = '';
    elseif target_time_step>= 10000
        disp('These are too many time steps. Code runs with time step == 9999')
        target_string_no = '';
        preceding_string_no = '';
        target_time_step = 9999;
    end
    
    
    if transient
        target_time_step_suffix = strcat(target_string_no,num2str(target_time_step));
        precedent_time_step_suffix = strcat(preceding_string_no,num2str(target_time_step-1));
    else
        target_time_step_suffix = strcat(target_string_no,num2str(target_time_step));
    end
    
    
    % read the binary heads
    filename3 = strcat(prefix,'o.head_pm.',target_time_step_suffix);
    [head] = read_HGS_binary_heads(dim,no_nodes,filename3);
    
    if transient
        
        filename3_2 = strcat(prefix,'o.head_pm.',precedent_time_step_suffix);
        
        [head_2] = read_HGS_binary_heads(dim,no_nodes,filename3_2);
        mean_head = (head + head_2)/2;
        
    end
    
    pressure_head = head - coordinates(:,3); % NOT NEEDED for flux reconstruction..
    
    %% NEW Surface domain, be carefull when using because it is still not mimicking the HGS primal solution in all details
    
    if surface_domain
        
        if transient
            filename_exflu1 = strcat(prefix,'o.ExchFlux_olf.',precedent_time_step_suffix);
            [exchange_flux1] = read_HGS_binary_exchangefluxes(dim,no_top_nodes,filename_exflu1);
            
            filename_exflu2 = strcat(prefix,'o.ExchFlux_olf.',target_time_step_suffix);
            [exchange_flux2] = read_HGS_binary_exchangefluxes(dim,no_top_nodes,filename_exflu2);
            
            exchange_flux = (exchange_flux1 + exchange_flux2)/2;
            
            filename_ETTot1 = strcat(prefix,'o.ETTotal_olf.',precedent_time_step_suffix);
            [ETTotal1] =  read_HGS_binary_totalET(dim,no_top_nodes,filename_ETTot1);
            
            filename_ETTot2 = strcat(prefix,'o.ETTotal_olf.',target_time_step_suffix);
            [ETTotal2] =  read_HGS_binary_totalET(dim,no_top_nodes,filename_ETTot2);
            
            ETTotal = (ETTotal1 + ETTotal2)/2;
            
            
        else
            filename_exflu = strcat(prefix,'o.ExchFlux_olf.',target_time_step_suffix);
            [exchange_flux] = read_HGS_binary_exchangefluxes(dim,no_top_nodes,filename_exflu);
            
            filename_ETTot = strcat(prefix,'o.ETTotal_olf.',target_time_step_suffix);
            [ETTotal] =  read_HGS_binary_totalET(dim,no_top_nodes,filename_ETTot);
            
            filename_rain = strcat(prefix,'o.rain_olf.',target_time_step_suffix);
            [nodal_rain] = read_HGS_binary_rain(dim,no_top_nodes,filename_rain);
            
            filename_ETEvap_olf = strcat(prefix,'o.ETEvap_olf.',target_time_step_suffix);
            [ETEvap_olf] = read_HGS_binary_ETEvap_olf(dim,no_top_nodes,filename_ETEvap_olf);
            
            filename_ETPmEvap_olf = strcat(prefix,'o.ETPmEvap_olf.',target_time_step_suffix);
            [ETPmEvap_olf] = read_HGS_binary_ETPmEvap_olf(dim,no_top_nodes+1000,filename_ETPmEvap_olf);
            
            filename_ETPmTranspire_olf = strcat(prefix,'o.ETPmTranspire_olf.',target_time_step_suffix);
            [ETPmTranspire_olf] = read_HGS_binary_ETPmTranspire_olf(dim,no_top_nodes,filename_ETPmTranspire_olf);
            
            ET_porous_medium = (ETPmTranspire_olf + ETPmEvap_olf);
        end
        
    end
    
    
    %% Read relative permeabilities
    
    if unsaturated_flow
        
        filenameN = strcat(prefix,'o.krw_pm.',target_time_step_suffix);
        
        [k_rels_bin] = read_HGS_binary_krels(no_nodes_per_element,no_elem,filenameN);
        k_rel = k_rels_bin;
        
        if transient
            filenameN_2 = strcat(prefix,'o.krw_pm.',precedent_time_step_suffix);
            [k_rels_bin_2] = read_HGS_binary_krels(no_nodes_per_element,no_elem,filenameN_2);
            k_rel = (k_rel + k_rels_bin_2)/2;
        end
        
    else
        k_rel = ones(no_elem,no_nodes_per_element);
    end
    
    pressure_head_actual = pressure_head;
    
    % van Genuchten prerequisites
    
    if ~theta_from_file
        theta_r = Sw_r.*n;
        
        rr_before = reshape(node_ids,[1,no_elem*no_nodes_per_element]);
        pressure_head_actual_mat_v = pressure_head(rr_before);
        pressure_head_actual_mat = reshape(pressure_head_actual_mat_v,[no_elem,no_nodes_per_element]);
        
        pressure_head_actual_mat_all = pressure_head_actual_mat;
        pressure_head_actual_centroid = mean(pressure_head_actual_mat_all,2);
        
        pressure_head_actual_mat(pressure_head_actual_mat>0) = 0;
        
        Se_mat_actual = (1 + (alpha.* (-1*pressure_head_actual_mat)).^N).^((1-N)./N);
        theta_actual_mat_recon_hc = (n-theta_r).*Se_mat_actual + theta_r;
        
        theta_actual = mean(theta_actual_mat_recon_hc,2);
    else
        filenameSW = strcat(prefix,'o.sw_pm.',target_time_step_suffix);
        [nodal_sat_bin_actual] = read_HGS_binary_nodal_sw(no_nodes_per_element,no_elem,filenameSW);
        theta_actual = nodal_sat_bin_actual.*n;
    end
    
    mean_theta = theta_actual;
    %% Dirichlet boundaries, TO ADAPT SECTION, Especially if transient==true
    
    % extract all nodes which are part of the Dirichlet boundary
    filename_lbc1 = strcat(prefix,'o.node_set.bc_ammer_in');
    filename_lbc2 = strcat(prefix,'o.node_set.bc_ammer_out');
    filename_lbc3 = strcat(prefix,'o.node_set.bc_neckar');
    
    inds_lbc1 = importdata(filename_lbc1);
    inds_lbc1 = inds_lbc1(2:end);
    
    inds_lbc2 = importdata(filename_lbc2);
    inds_lbc2 = inds_lbc2(2:end);
    
    inds_lbc3 = importdata(filename_lbc3);
    inds_lbc3 = inds_lbc3(2:end);
    
    diri_nodes = [inds_lbc1;inds_lbc2;inds_lbc3];
    diri_nodes = unique(diri_nodes);
    
    
    %% Neumann unlike zero BC
    
    neumann_unlike0_exists = false;
    
    if neumann_unlike0_exists
        q_Neumann_r = -1*10^-6*ones(6,1); % Normal flux on boundary face in ascending order in analogy to the element ID
        % negative flux = flux leaving the domain; positive flux = flux
        % entering the domain
    end
    
    %% RECHARGE, TO ADAPT
    
    recharge = true;
    
    recharge_as_rain = false; % TO ADAPT LINE, TYPE OF BC
    recharge_as_flux = true; % TO ADAPT LINE, TYPE OF BC
    
    inds_recharge = [];
    
    recharge_surf_q = zeros(no_elem,1);
    
    if recharge
        
        filename_recharge1 = strcat(prefix,'o.face_set.rr_cropland');
        inds_recharge1 = load(filename_recharge1);
        inds_recharge1 = inds_recharge1(:,1);
        recharge_surf_q(inds_recharge1) = 1.777660e-09;
        
        filename_recharge2 = strcat(prefix,'o.face_set.rr_floodplain');
        inds_recharge2 = importdata(filename_recharge2);
        inds_recharge2 = inds_recharge2(:,1);
        recharge_surf_q(inds_recharge2) = 1.187748e-09;
        %
        filename_recharge3 = strcat(prefix,'o.face_set.rr_sandstone');
        inds_recharge3 = importdata(filename_recharge3);
        inds_recharge3 = inds_recharge3(:,1);
        recharge_surf_q(inds_recharge3) = 5.159632e-10;
        
        filename_recharge4 = strcat(prefix,'o.face_set.rr_urban');
        inds_recharge4 = importdata(filename_recharge4);
        inds_recharge4 = inds_recharge4(:,1);
        recharge_surf_q(inds_recharge4) = 1.489827e-09;
        
        inds_recharge = [inds_recharge1;inds_recharge2;inds_recharge3;inds_recharge4]; % list of all elements exhibiting recharge
    end
    
    
    %% injection/extraction wells:  Well-Aenderung
    
    wells_exist = false;
    
    if wells_exist
        wells_name = strcat(prefix,'o.node_set.well');
        wells_nodes_v = load(wells_name);
        
        wells_name2 = strcat(prefix,'o.node_set.well2');
        wells_nodes_v2 = load(wells_name2);
        
        wells_nodes = [wells_nodes_v(2:end) wells_nodes_v2(2:end)];
        
        list_of_total_well_fluxes = [-1.0e-4 +1.0e-4]; % [m^3/s]
        
        Vol_prisms_name = strcat('Vol_prisms_',prefix,'.mat');
        Vol_prisms = importdata(Vol_prisms_name);
        
    end
    
    
    if transient && target_time_step>=2
        
        endtime = target_times(target_time_step);
        beforetime = target_times(target_time_step-1);
        delta_t = endtime-beforetime;
        
        % pressure head before:
        
        filenameheadbefore = strcat(prefix,'o.head_pm.',precedent_time_step_suffix);
        dim = 1;
        [head_before] = read_HGS_binary_heads(dim,no_nodes,filenameheadbefore);
        pressure_head_before = head_before - coordinates(:,3);
        
        if ~theta_from_file
            pressure_head_before_mat_v = pressure_head_before(rr_before);
            pressure_head_before_mat = reshape(pressure_head_before_mat_v,[no_elem,no_nodes_per_element]);
            
            %%
            pressure_head_before_mat_all = pressure_head_before_mat;
            pressure_head_before_centroid = mean(pressure_head_before_mat_all,2);
            pressure_head_before_mat(pressure_head_before_mat>0) = 0;
            
            Se_mat_before = (1 + (alpha.* (-1*pressure_head_before_mat)).^N).^((1-N)./N);
            theta_before_mat_recon_hc = (n-theta_r).*Se_mat_before + theta_r;
            theta_before = mean(theta_before_mat_recon_hc,2);
        else
            filenameswbefore = strcat(prefix,'o.sw_pm.',precedent_time_step_suffix);
            [nodal_sat_bin_before] = read_HGS_binary_nodal_sw(no_nodes_per_element,no_elem,filenameswbefore);
            theta_before = nodal_sat_bin_before.*n;
            
        end
        
        diff_thetas = theta_actual - theta_before;
        diff_pressure_heads = pressure_head_actual_centroid - pressure_head_before_centroid;
        mean_theta = (theta_actual + theta_before)/2;
        
        
        lump_transflow_inrhs =  -(   (diff_thetas./delta_t)   +  (((((mean_theta))./n).*S_sp ).*(diff_pressure_heads./delta_t))      ).*Vol_prisms; % NEW CHANGE 12.01.2020, replaced theta_actual by mean_theta
        
    end
    
    
    %% DITCHES, TO ADAPT SECTION
    
    ditches_exist = true;
    
    elem_list_ditch = [];
    
    if ditches_exist
        
        filename_ditch1_nodes = strcat(prefix,'o.node_set.bc_ditches_and_arbach');
        n_inds_ditch1 = importdata(filename_ditch1_nodes);
        n_inds_ditch1 = n_inds_ditch1(2:end);
        
        crit_head_ditch1 = 0.001;
        conductance_ditch1 = 1e-1;
        %
        %     filename_ditch2_nodes = strcat(prefix,'o.node_set.riverStretch_2');
        %     n_inds_ditch2 = importdata(filename_ditch2_nodes);
        %     n_inds_ditch2 = n_inds_ditch2(2:end);
        %
        %     crit_head_ditch2 = 0.16852;
        %     conductance_ditch2 = 1e-2;
        %
        %     filename_ditch3_nodes = strcat(prefix,'o.node_set.riverStretch_3');
        %     n_inds_ditch3 = importdata(filename_ditch3_nodes);
        %     n_inds_ditch3 = n_inds_ditch3(2:end);
        %
        %     crit_head_ditch3 = 0.16852;
        %     conductance_ditch3 = 1e-2;
        
        len_inds_ditches = [length(n_inds_ditch1)];%,length(n_inds_ditch2),length(n_inds_ditch3)];
        
        crit_heads_ditch_all = [crit_head_ditch1];%,crit_head_ditch2,crit_head_ditch3];
        n_inds_ditches_all = {n_inds_ditch1};%,n_inds_ditch2,n_inds_ditch3};
        conduct_ditches_all = [conductance_ditch1];%,conductance_ditch2,conductance_ditch3];
    else
        len_inds_ditches = [];
        crit_heads_ditch_all = [];
        n_inds_ditches_all = {};
        conduct_ditches_all = [];
    end
    
    %% Leakage BC (i.e., third-type, e.g., rivers), only On TOP of domain, TO ADAPT SECTION
    filename_river = strcat(prefix,'o.node_set.river_set');
    n_inds_leak1 = importdata(filename_river);
    n_inds_leak1 = n_inds_leak1(2:end);
    
    river_data = load('bc_river.lst');
    
    ref_head_leak1 = river_data(:,1); % vector
    conductance_leak1 = river_data(:,2); % vector
    
    n_inds_leak2 = []; % vector
    
    ref_head_leak2 = []; % vector
    conductance_leak2 = []; % vector
    
    len_inds_leak = [length(n_inds_leak1)];%,length(n_inds_leak2)];
    
    ref_heads_leak_all = {ref_head_leak1};%,ref_head_leak2};
    n_inds_leak_all = {n_inds_leak1};%,n_inds_leak2};
    conduct_leak_all = {conductance_leak1};%,conductance_leak2};
    
    %% Lateral Leakage BC (i.e., third type, hydrosystem fluid interaction/transfer), only on lateral faces of domain, TO ADAPT SECTION
    n_inds_leak1_lat = [];%importdata('southRiver.asc.txt');
    elem_inds_leak1_lat = [];%importdata('KaesBachLite_v1o.face_set.outlet_nodes');
    elem_inds_leak1_lat = [];%sort(elem_inds_leak1_lat(:,1));
    
    A_lat_robin = [];%zeros(length(elem_inds_leak1_lat),1);
    ref_head_leak1_lat = [];%352.2306; % m.a.s.l.
    k_leak1_lat = [];%1e-5; % [m/s]
    dist_leak1_lat = [];%1e3; % [m]
    fact_leak1_lat = (k_leak1_lat/dist_leak1_lat);
    
    
    n_inds_leak2_lat = []; % vector
    elem_inds_leak2_lat = [];
    
    ref_head_leak2_lat = []; % vector
    k_leak2_lat = []; % [m/s]
    dist_leak2_lat = []; % [m]
    fact_leak2_lat = k_leak2_lat/dist_leak2_lat;
    
    len_elem_inds_leak_lat = [length(elem_inds_leak1_lat);length(elem_inds_leak2_lat)];
    len_inds_leak_lat = [length(n_inds_leak1_lat),length(n_inds_leak2_lat)];
    
    ref_heads_leak_all_lat = {ref_head_leak1_lat,ref_head_leak2_lat};
    
    ref_heads_leak_all_lat_vec = {ref_head_leak1_lat,ref_head_leak2_lat};
    k_leak1_lat_vec = {k_leak1_lat,k_leak2_lat};
    dist_leak_lat_robin_vec = {dist_leak1_lat,dist_leak2_lat};
    
    n_inds_leak_all_lat = {n_inds_leak1_lat,n_inds_leak2_lat};
    elem_inds_leak_all_lat = {elem_inds_leak1_lat,elem_inds_leak2_lat};
    fact_leak_all_lat = {fact_leak1_lat,fact_leak2_lat};
    
    
    %% Initialization
    
    elem_list_robin =zeros((sum(len_inds_ditches) + sum(len_inds_leak))*12 + sum(len_elem_inds_leak_lat),1);
    rhs_elems_robin = zeros((sum(len_inds_ditches) + sum(len_inds_leak))*12 + sum(len_elem_inds_leak_lat),1);
    mob_coeff_elems_robin = zeros((sum(len_inds_ditches) + sum(len_inds_leak))*12 + sum(len_elem_inds_leak_lat),1);
    
    robin_elem_counter = 1;
    
    %% Lateral Leakage BC (i.e., third type, hydrosystem fluid interaction/transfer), only on lateral faces of domain
    
    right_areas = zeros(no_nodes,1);
    
    for dd = 1:length(len_inds_leak_lat)
        n_inds_leak_current =  n_inds_leak_all_lat{dd};
        elem_inds_leak_current = elem_inds_leak_all_lat{dd};
        
        ref_head_current = ref_heads_leak_all_lat{dd};
        leak_fact_current = fact_leak_all_lat{dd};
        
        for kk = 1:length(elem_inds_leak_all_lat{dd})
            current_element = elem_inds_leak_current(kk);
            
            node_id1 = node_ids(current_element,1);
            node_id2 = node_ids(current_element,2);
            node_id3 = node_ids(current_element,3);
            node_id4 = node_ids(current_element,4);
            node_id5 = node_ids(current_element,5);
            node_id6 = node_ids(current_element,6);
            
            find_node_id1 = find(node_id1==n_inds_leak_current);
            find_node_id2 = find(node_id2==n_inds_leak_current);
            find_node_id3 = find(node_id3==n_inds_leak_current);
            
            if isempty(find_node_id1)
                
                x1 = coordinates(node_id2,1); x2 = coordinates(node_id3,1);
                y1 = coordinates(node_id2,2); y2 = coordinates(node_id3,2);
                z1 = coordinates(node_id2,3); z2 = coordinates(node_id3,3);
                z3 = coordinates(node_id5,3); z4 = coordinates(node_id6,3);
                
            elseif isempty(find_node_id2)
                
                x1 = coordinates(node_id1,1); x2 = coordinates(node_id3,1);
                y1 = coordinates(node_id1,2); y2 = coordinates(node_id3,2);
                z1 = coordinates(node_id1,3); z2 = coordinates(node_id3,3);
                z3 = coordinates(node_id4,3); z4 = coordinates(node_id6,3);
                
            elseif isempty(find_node_id3)
                
                x1 = coordinates(node_id1,1); x2 = coordinates(node_id2,1);
                y1 = coordinates(node_id1,2); y2 = coordinates(node_id2,2);
                z1 = coordinates(node_id1,3); z2 = coordinates(node_id2,3);
                z3 = coordinates(node_id4,3); z4 = coordinates(node_id5,3);
                
            end
            
            LQ = sqrt((x1 - x2)^2 + (y1 - y2)^2);
            hQ_d = (z1 + z2)/2;
            hQ_u = (z3 + z4)/2;
            hQ = sqrt((hQ_d - hQ_u)^2);
            A_elem_leak = LQ*hQ;
            
            elem_list_robin(robin_elem_counter) = current_element;
            rhs_elems_robin(robin_elem_counter)  =  A_elem_leak*leak_fact_current*ref_head_current;
            mob_coeff_elems_robin(robin_elem_counter) = A_elem_leak*leak_fact_current;
            
            robin_elem_counter = robin_elem_counter + 1;
        end
        
        
    end
    
    
    
    for dd = 1:length(len_inds_ditches)
        n_inds_ditch_current =  n_inds_ditches_all{dd};
        
        crit_head_ditch_cur = crit_heads_ditch_all(dd);
        conductance_ditch_cur = conduct_ditches_all(dd);
        
        
        for kk = 1:len_inds_ditches(dd)
            
            
            cur_ditch_node = n_inds_ditch_current(kk);
            
            %%%if ~enforce_unsatflow_correctness
            
            if (head(cur_ditch_node)>=coordinates(cur_ditch_node,3)+crit_head_ditch_cur)
                find_topele_node4 = find(cur_ditch_node==node_ids((no_layers-1)*elem_per_layer+1:(no_layers)*elem_per_layer,4));
                find_topele_node5 = find(cur_ditch_node==node_ids((no_layers-1)*elem_per_layer+1:(no_layers)*elem_per_layer,5));
                find_topele_node6 = find(cur_ditch_node==node_ids((no_layers-1)*elem_per_layer+1:(no_layers)*elem_per_layer,6));
                
                all_elem_ids_ditchn = [find_topele_node4; find_topele_node5; find_topele_node6]+(no_layers-1)*elem_per_layer;
                elem_list_robin(robin_elem_counter:robin_elem_counter+length(all_elem_ids_ditchn)-1) = all_elem_ids_ditchn;
                
                A_horiz_relevant = A_horiz(all_elem_ids_ditchn);
                w_A_horiz = A_horiz_relevant/sum(A_horiz_relevant);
                
                rhs_elems_ditch_v = (conductance_ditch_cur*(coordinates(cur_ditch_node,3)+crit_head_ditch_cur))*w_A_horiz;
                
                rhs_elems_robin(robin_elem_counter:robin_elem_counter+length(all_elem_ids_ditchn)-1)  =  rhs_elems_ditch_v;
                
                mob_coeff_elems_ditch_v = conductance_ditch_cur*w_A_horiz;
                mob_coeff_elems_robin(robin_elem_counter:robin_elem_counter+length(all_elem_ids_ditchn)-1) = mob_coeff_elems_ditch_v;
                
                robin_elem_counter = robin_elem_counter + length(all_elem_ids_ditchn);
            end
            
            %%%elseif enforce_unsatflow_correctness
            % % %             if (head(cur_ditch_node)>=coordinates(cur_ditch_node,3)+crit_head_ditch_cur) || dd==length(len_inds_ditches)
            % % %                 find_topele_node4 = find(cur_ditch_node==node_ids((no_layers-1)*elem_per_layer+1:(no_layers)*elem_per_layer,4));
            % % %                 find_topele_node5 = find(cur_ditch_node==node_ids((no_layers-1)*elem_per_layer+1:(no_layers)*elem_per_layer,5));
            % % %                 find_topele_node6 = find(cur_ditch_node==node_ids((no_layers-1)*elem_per_layer+1:(no_layers)*elem_per_layer,6));
            % % %
            % % %                 all_elem_ids_ditchn = [find_topele_node4; find_topele_node5; find_topele_node6]+(no_layers-1)*elem_per_layer;
            % % %                 elem_list_robin(robin_elem_counter:robin_elem_counter+length(all_elem_ids_ditchn)-1) = all_elem_ids_ditchn;
            % % %
            % % %                 A_horiz_relevant = A_horiz(all_elem_ids_ditchn);
            % % %                 w_A_horiz = A_horiz_relevant/sum(A_horiz_relevant);
            % % %
            % % %                 if dd==length(len_inds_ditches)
            % % %                     rhs_elems_ditch_v = (conductance_ditch_cur*(head(cur_ditch_node)))*w_A_horiz;
            % % %                 else
            % % %                     rhs_elems_ditch_v = (conductance_ditch_cur*(coordinates(cur_ditch_node,3)+crit_head_ditch_cur))*w_A_horiz;
            % % %                 end
            % % %
            % % %                 rhs_elems_robin(robin_elem_counter:robin_elem_counter+length(all_elem_ids_ditchn)-1)  =  rhs_elems_ditch_v;
            % % %
            % % %                 mob_coeff_elems_ditch_v = conductance_ditch_cur*w_A_horiz;
            % % %                 mob_coeff_elems_robin(robin_elem_counter:robin_elem_counter+length(all_elem_ids_ditchn)-1) = mob_coeff_elems_ditch_v;
            % % %
            % % %                 robin_elem_counter = robin_elem_counter + length(all_elem_ids_ditchn);
            % % %             end
            %%%end
            
        end
    end
    
    %% Leakage BC (i.e., third-type, e.g., rivers), only on top of domain
    
    for dd = 1:length(len_inds_leak)
        n_inds_leak_current = n_inds_leak_all{dd};
        
        crit_head_leak_cur = ref_heads_leak_all{dd};
        conductance_leak_cur = conduct_leak_all{dd};
        for kk = 1:len_inds_leak(dd)
            
            
            cur_leak_node = n_inds_leak_current(kk);
            
            
            find_topele_node4 = find(cur_leak_node==node_ids((no_layers-1)*elem_per_layer+1:(no_layers)*elem_per_layer,4));
            find_topele_node5 = find(cur_leak_node==node_ids((no_layers-1)*elem_per_layer+1:(no_layers)*elem_per_layer,5));
            find_topele_node6 = find(cur_leak_node==node_ids((no_layers-1)*elem_per_layer+1:(no_layers)*elem_per_layer,6));
            
            all_elem_ids_leakn = [find_topele_node4; find_topele_node5; find_topele_node6]+(no_layers-1)*elem_per_layer;
            elem_list_robin(robin_elem_counter:robin_elem_counter+length(all_elem_ids_leakn)-1) = all_elem_ids_leakn;
            
            A_horiz_relevant = A_horiz(all_elem_ids_leakn);
            w_A_horiz = A_horiz_relevant/sum(A_horiz_relevant);
            
            rhs_elems_leak_v = (conductance_leak_cur(kk)*(crit_head_leak_cur(kk)))*w_A_horiz;
            
            rhs_elems_robin(robin_elem_counter:robin_elem_counter+length(all_elem_ids_leakn)-1)  =  rhs_elems_leak_v;
            
            mob_coeff_elems_leak_v = conductance_leak_cur(kk)*w_A_horiz;
            mob_coeff_elems_robin(robin_elem_counter:robin_elem_counter+length(all_elem_ids_leakn)-1) = mob_coeff_elems_leak_v;
            
            robin_elem_counter = robin_elem_counter + length(all_elem_ids_leakn);
            
            
            
        end
    end
    
    
    if surface_domain
        Q_inex_olf_pm = zeros(no_elem,1);
        n_inds_top_IDs =  [no_nodes - no_top_nodes + 1:1:no_nodes]';
        
        for curEL = no_elem-elem_per_layer +1:no_elem
            node_ids_el = node_ids(curEL,:);
            rel_ID1 = node_ids_el(4) - (no_nodes - no_top_nodes);
            rel_ID2 = node_ids_el(5) - (no_nodes - no_top_nodes);
            rel_ID3 = node_ids_el(6) - (no_nodes - no_top_nodes);
            
            Q_inex_olf_pm(curEL) = (((-1*((exchange_flux(rel_ID1) + exchange_flux(rel_ID2) + exchange_flux(rel_ID3))/3))+((ET_porous_medium(rel_ID1) + ET_porous_medium(rel_ID2) + ET_porous_medium(rel_ID3))/3 )) )*A_horiz(curEL);
        end
        
    end
    
    
    %% Finalize it
    
    elem_list_robin(elem_list_robin==0) = [];
    rhs_elems_robin(length(elem_list_robin)+1:end) = [];
    mob_coeff_elems_robin(length(elem_list_robin)+1:end) = [];
    
    disp('Prestuff done!')
    
    %% Flux correction itself
    no_Diri_elem = length(vector_of_d_bounds_mult);
    elem_diri_bound = vector_of_d_bounds_mult;
    counter_D_val_bounds = 1;
    
    no_Neumann_elem = length(vector_of_n_bounds_mult);
    elem_neumann_bound = vector_of_n_bounds_mult;
    counter_N_val_bounds = 1;
    
    no_poss_neighbs = 5;
    
    K_aniso = zeros(no_elem,5);
    
    vector_of_D_val_bounds = zeros(no_elem,1);
    rhs_vector_FVM = zeros(no_elem + no_Diri_elem + no_Neumann_elem,1);
    only_recharge = zeros(no_elem + no_Diri_elem + no_Neumann_elem,1);
    
    centro_head = zeros(no_elem,1);
    
    diri_neum_counter = 1;
    diri_counter = 1;
    neum_counter = 1;
    rob = 1;
    
    for j = 1:no_elem
        
        % Dirichlet Assignment
        
        heads_Diri = zeros(1,6);
        
        node_ids1 = node_ids(j,1);
        node_ids2 = node_ids(j,2);
        node_ids3 = node_ids(j,3);
        node_ids4 = node_ids(j,4);
        node_ids5 = node_ids(j,5);
        node_ids6 = node_ids(j,6);
        
        if ~transient
            centro_head(j) = (head(node_ids1) + head(node_ids2) + head(node_ids3) + head(node_ids4) + head(node_ids5) + head(node_ids6))/6;
        else
            centro_head(j) = (mean_head(node_ids1) + mean_head(node_ids2) + mean_head(node_ids3) + mean_head(node_ids4) + mean_head(node_ids5) + mean_head(node_ids6))/6;
        end
        
        
        d_nodes = intersect(node_ids(j,:),diri_nodes);
        
        
        
        find_jD = find(j==elem_diri_bound);
        
        if ~isempty(find_jD)
            
            vector_of_D_val_bounds_intmed = zeros(1,5);
            no_dirifaces_j = face_index_mat(j,:);
            no_dirifaces_j(no_dirifaces_j==0) = [];
            
            for kk = 1:length(no_dirifaces_j)
                
                if face_flag_D_N(j,kk)==1
                    
                    act_diri_face = no_dirifaces_j(kk);
                    
                    if act_diri_face==1
                        d_nodes_act = [node_ids2;node_ids3;node_ids5;node_ids6];
                        vector_of_D_val_bounds(counter_D_val_bounds) = mean(head(d_nodes_act));
                        vector_of_D_val_bounds_intmed(kk) = vector_of_D_val_bounds(counter_D_val_bounds);
                        counter_D_val_bounds = counter_D_val_bounds + 1;
                        %                 k_rel(j,[2,3,5,6]) = 1;
                    elseif act_diri_face==2
                        d_nodes_act = [node_ids1;node_ids3;node_ids4;node_ids6];
                        vector_of_D_val_bounds(counter_D_val_bounds) = mean(head(d_nodes_act));
                        vector_of_D_val_bounds_intmed(kk) = vector_of_D_val_bounds(counter_D_val_bounds);
                        counter_D_val_bounds = counter_D_val_bounds + 1;
                        %                 k_rel(j,[1,3,4,6]) = 1;
                    elseif act_diri_face==3
                        d_nodes_act = [node_ids1;node_ids2;node_ids4;node_ids5];
                        vector_of_D_val_bounds(counter_D_val_bounds) = mean(head(d_nodes_act));
                        vector_of_D_val_bounds_intmed(kk) = vector_of_D_val_bounds(counter_D_val_bounds);
                        counter_D_val_bounds = counter_D_val_bounds + 1;
                        %                 k_rel(j,[1,2,4,5]) = 1;
                    elseif act_diri_face==4
                        d_nodes_act = [node_ids1;node_ids2;node_ids3];
                        vector_of_D_val_bounds(counter_D_val_bounds) = mean(head(d_nodes_act));
                        vector_of_D_val_bounds_intmed(kk) = vector_of_D_val_bounds(counter_D_val_bounds);
                        counter_D_val_bounds = counter_D_val_bounds + 1;
                        %                 k_rel(j,[1,2,3]) = 1;
                    elseif act_diri_face==5
                        d_nodes_act = [node_ids4;node_ids5;node_ids6];
                        vector_of_D_val_bounds(counter_D_val_bounds) = mean(head(d_nodes_act));
                        vector_of_D_val_bounds_intmed(kk) = vector_of_D_val_bounds(counter_D_val_bounds);
                        counter_D_val_bounds = counter_D_val_bounds + 1;
                        %                 k_rel(j,[4,5,6]) = 1;
                    end
                end
                
                
            end
            vector_of_D_val_bounds_intmed(vector_of_D_val_bounds_intmed==0) = [];
            
        end
        
        % k_rel-Aederung
        is_recharge_ele = find(j==inds_recharge);
        
        if ~isempty(is_recharge_ele)
            k_rel(j,4:6) = [1,1,1];
        end
        
        
        %% Assigning hydraulic conductivity values
        
        for k = 1:no_poss_neighbs
            
            K_main_dim = true; % SET THIS TO TRUE. Anisotropy is not yet included in the code. However, HGS is also
            % doing some approximations.
            if K_main_dim
                if k<=3
                    K_aniso(j,k) = k_xx(j);
                else
                    K_aniso(j,k) = k_zz(j);
                end
                
            end
            
        end
        %% Recharge and Robin-BC
        
        robin_elems = find(j==elem_list_robin);
        rhs_elems_robin_j = 0;
        
        if ~isempty(robin_elems)
            rhs_elems_robin_j = rhs_elems_robin(robin_elems);
            rhs_vector_FVM(j) = sum(rhs_elems_robin_j);
            rob = rob + 1;
        end
        
        
        if j >=(no_layers-1)*elem_per_layer+1
            
            is_j_recharge_elem = find(j==inds_recharge);
            
            
            if ~isempty(is_j_recharge_elem)
                if recharge_as_rain
                    rhs_vector_FVM(j) = rhs_vector_FVM(j) + recharge_surf_q(j)*A_horiz(j);
                    only_recharge(j) = recharge_surf_q(j)*A_horiz(j);
                elseif recharge_as_flux
                    rhs_vector_FVM(j) = rhs_vector_FVM(j) + recharge_surf_q(j)*A_top(j);
                    only_recharge(j) = recharge_surf_q(j)*A_top(j);
                end
                
            end
            
            % SURFACE DOMAIN COUPLING
            if surface_domain
                rhs_vector_FVM(j) = rhs_vector_FVM(j) + Q_inex_olf_pm(j);
            end
            
        end
        %% Neumann flux right-hand side
        find_jN = find(j==elem_neumann_bound);
        
        if ~isempty(find_jN)
            
            no_neumannfaces_j = face_index_mat(j,:);
            no_neumannfaces_j(no_neumannfaces_j==0) = [];
            for kk = 1:length(no_neumannfaces_j)
                if face_flag_D_N(j,kk)==2
                    ind_Neum_face = find(all_neighb_indices(j,:)>no_elem);
                    A_f_Neumann = A_f_all(j,ind_Neum_face);
                    
                    rhs_vector_FVM(no_elem + diri_neum_counter) = A_f_Neumann*q_Neumann_r(neum_counter);
                    neum_counter = neum_counter + 1;
                    diri_neum_counter = diri_neum_counter + 1;
                end
            end
            
        end
        
        if ~isempty(find_jD)
            for kk = 1:length(no_dirifaces_j)
                if face_flag_D_N(j,kk)==1
                    rhs_vector_FVM(no_elem + diri_neum_counter) = vector_of_D_val_bounds_intmed(kk);
                    diri_neum_counter = diri_neum_counter + 1;
                    diri_counter = diri_counter + 1;
                end
            end
        end
    end
    
    vector_of_D_val_bounds(vector_of_D_val_bounds==0) = [];
    
    %% injection/extraction wells: Well-Aenderung
    
    if wells_exist
        
        ele_list_wells = cell(1,length(wells_nodes));
        
        for ww = 1:length(wells_nodes)
            
            well_node = wells_nodes(ww);
            
            find1 = find(well_node==node_ids(:,1));
            find2 = find(well_node==node_ids(:,2));
            find3 = find(well_node==node_ids(:,3));
            find4 = find(well_node==node_ids(:,4));
            find5 = find(well_node==node_ids(:,5));
            find6 = find(well_node==node_ids(:,6));
            
            ele_list_well = [find1;find2;find3;find4;find5;find6];
            ele_list_wells{ww} = ele_list_well;
            
            prismatic_volumes = Vol_prisms(ele_list_well);
            
            vol_weight_facs = prismatic_volumes/sum(prismatic_volumes);
            
            Q_well = list_of_total_well_fluxes(ww);
            tot_elem_well_discharge = Q_well*vol_weight_facs;
            
            rhs_vector_FVM(ele_list_well) = rhs_vector_FVM(ele_list_well) + tot_elem_well_discharge;
            
        end
        
        ele_list_wells_name = strcat('ele_list_wells_',prefix,'.mat');
        save(ele_list_wells_name,'ele_list_wells')
        
    end
    
    
    mob_a = zeros(no_elem,no_poss_neighbs);
    
    k_rel_F_all = zeros(no_elem,no_poss_neighbs);
    
    disp('Boundary computations done!')
    
    %% mobility coefficients for the FVM scheme:
    Neumann_counter_place = 1;
    Bound_counter_place = 1;
    
    all_bounds_DiriNeum = sort([elem_diri_bound;elem_neumann_bound]);
    
    mob_matrix_diag = NaN(no_elem,1);
    
    row_mob_matrix_offdiag = NaN(no_elem*no_poss_neighbs,1);
    col_mob_matrix_offdiag = NaN(no_elem*no_poss_neighbs,1);
    mob_matrix_offdiag = NaN(no_elem*no_poss_neighbs,1);
    off_diag_counter = 1;
    
    row_mob_matrix_diri = NaN(no_Diri_elem,1);
    col_mob_matrix_diri = NaN(no_Diri_elem,1);
    diri_mat_counter = 1;
    
    row_mob_matrix_neum_diag = NaN(no_Neumann_elem,1);
    col_mob_matrix_neum_diag = NaN(no_Neumann_elem,1);
    row_mob_matrix_neum_offdiag = NaN(no_Neumann_elem,1);
    col_mob_matrix_neum_offdiag = NaN(no_Neumann_elem,1);
    neum_mat_counter = 1;
    
    starting_vec_solver = ones(no_elem + no_Diri_elem + no_Neumann_elem,1);
    
    for i = 1:no_elem + no_Diri_elem + no_Neumann_elem
        
        if i <=no_elem
            
            starting_vec_solver(i) = centro_head(i);
            
            %% facial upstream of relative permeability
            if unsaturated_flow
                
                node_ids1 = node_ids(i,1);
                node_ids2 = node_ids(i,2);
                node_ids3 = node_ids(i,3);
                node_ids4 = node_ids(i,4);
                node_ids5 = node_ids(i,5);
                node_ids6 = node_ids(i,6);
                
                k_relF1_v1 = [k_rel(i,2) , k_rel(i,3) , k_rel(i,5) , k_rel(i,6)];
                k_relF2_v1 = [k_rel(i,1) , k_rel(i,3) , k_rel(i,4) , k_rel(i,6)];
                k_relF3_v1 = [k_rel(i,1) , k_rel(i,2) , k_rel(i,4) , k_rel(i,5)];
                k_relF4_v1 = [k_rel(i,1) , k_rel(i,2) , k_rel(i,3)];
                k_relF5_v1 = [k_rel(i,4) , k_rel(i,5) , k_rel(i,6)];
                
                %% k_rels neighbour
                neighb_elem_inds = all_neighb_indices(i,:);
                % elem_inds_neighb
                
                % node_inds neighb1: 2,3,5,6 of i shared by neighb 1
                act_neighb = neighb_elem_inds(1);
                if act_neighb~=0 && act_neighb<=no_elem
                    find2 = find(node_ids2==node_ids(act_neighb,:));
                    find3 = find(node_ids3==node_ids(act_neighb,:));
                    find5 = find2 + 3;
                    find6 = find3 + 3;
                    k_relF1_v2 = [k_rel(act_neighb,find2) , k_rel(act_neighb,find3) , k_rel(act_neighb,find5) , k_rel(act_neighb,find6)];
                else
                    k_relF1_v2 = [k_rel(i,2) , k_rel(i,3) , k_rel(i,5) , k_rel(i,6)];
                end
                % node_inds neighb2: 1,3,4,6 of i shared by neighb 2
                act_neighb = neighb_elem_inds(2);
                if act_neighb~=0 && act_neighb<=no_elem
                    find1 = find(node_ids1==node_ids(act_neighb,:));
                    find3 = find(node_ids3==node_ids(act_neighb,:));
                    find4 = find1 + 3;
                    find6 = find3 + 3;
                    k_relF2_v2 = [k_rel(act_neighb,find1) , k_rel(act_neighb,find3) , k_rel(act_neighb,find4) , k_rel(act_neighb,find6)];
                else
                    k_relF2_v2 = [k_rel(i,1) , k_rel(i,3) , k_rel(i,4) , k_rel(i,6)];
                end
                % node_inds neighb3: 1,2,4,5 of i shared by neighb 3
                act_neighb = neighb_elem_inds(3);
                if act_neighb~=0 && act_neighb<=no_elem
                    find1 = find(node_ids1==node_ids(act_neighb,:));
                    find2 = find(node_ids2==node_ids(act_neighb,:));
                    find4 = find1 + 3;
                    find5 = find2 + 3;
                    k_relF3_v2 = [k_rel(act_neighb,find1) , k_rel(act_neighb,find2) , k_rel(act_neighb,find4) , k_rel(act_neighb,find5)];
                else
                    k_relF3_v2 = [k_rel(i,1) , k_rel(i,2) , k_rel(i,4) , k_rel(i,5)];
                end
                % node_inds neighb4: 1,2,3 of i shared by neighb 4
                act_neighb = neighb_elem_inds(4);
                if act_neighb~=0 && act_neighb<=no_elem
                    k_relF4_v2 = [k_rel(act_neighb,4) , k_rel(act_neighb,5) , k_rel(act_neighb,6)];
                else
                    k_relF4_v2 = [k_rel(i,1) , k_rel(i,2) , k_rel(i,3)];
                end
                % node_inds neighb5: 4,5,6 of i shared by neighb 5
                act_neighb = neighb_elem_inds(5);
                if act_neighb~=0 && act_neighb<=no_elem
                    k_relF5_v2 = [k_rel(act_neighb,1) , k_rel(act_neighb,2) , k_rel(act_neighb,3)];
                else
                    k_relF5_v2 = [k_rel(i,4) , k_rel(i,5) , k_rel(i,6)];
                end
                
                
                
                
                %% weighting of relative permeability on the face
                if lincentro_upface
                    mean_F1_v1 = mean(k_relF1_v1);
                    mean_F2_v1 = mean(k_relF2_v1);
                    mean_F3_v1 = mean(k_relF3_v1);
                    mean_F4_v1 = mean(k_relF4_v1);
                    mean_F5_v1 = mean(k_relF5_v1);
                    
                    mean_F1_v2 = mean(k_relF1_v2);
                    mean_F2_v2 = mean(k_relF2_v2);
                    mean_F3_v2 = mean(k_relF3_v2);
                    mean_F4_v2 = mean(k_relF4_v2);
                    mean_F5_v2 = mean(k_relF5_v2);
                    
                    if neighb_elem_inds(1)~=0 && neighb_elem_inds(1)<=no_elem
                        
                        if centro_head(i)>=centro_head(neighb_elem_inds(1))
                            k_relF1 = mean_F1_v1;
                        else
                            k_relF1 = mean_F1_v2;
                        end
                        
                    else
                        k_relF1 = mean_F1_v1;
                    end
                    
                    if neighb_elem_inds(2)~=0 && neighb_elem_inds(2)<=no_elem
                        
                        if centro_head(i)>=centro_head(neighb_elem_inds(2))
                            k_relF2 = mean_F2_v1;
                        else
                            k_relF2 = mean_F2_v2;
                        end
                        
                    else
                        k_relF2 = mean_F2_v1;
                    end
                    
                    if neighb_elem_inds(3)~=0 && neighb_elem_inds(3)<=no_elem
                        
                        if centro_head(i)>=centro_head(neighb_elem_inds(3))
                            k_relF3 = mean_F3_v1;
                        else
                            k_relF3 = mean_F3_v2;
                        end
                        
                    else
                        k_relF3 = mean_F3_v1;
                    end
                    
                    if neighb_elem_inds(4)~=0 && neighb_elem_inds(4)<=no_elem
                        
                        if centro_head(i)>=centro_head(neighb_elem_inds(4))
                            k_relF4 = mean_F4_v1;
                        else
                            k_relF4 = mean_F4_v2;
                        end
                        
                    else
                        k_relF4 = mean_F4_v1;
                    end
                    
                    if neighb_elem_inds(5)~=0 && neighb_elem_inds(5)<=no_elem
                        
                        if centro_head(i)>=centro_head(neighb_elem_inds(5))
                            k_relF5 = mean_F5_v1;
                        else
                            k_relF5 = mean_F5_v2;
                        end
                        
                    else
                        k_relF5 = mean_F5_v1;
                    end
                    
                end
                
                k_rel_F_ord = [k_relF1,k_relF2,k_relF3,k_relF4,k_relF5];
                
                k_rel_F_all(i,:) = k_rel_F_ord;
                
            else
                k_rel_F_all(i,:) = ones(1,5);
            end
            
            %% mobility coefficients
            
            new_counter = 1;
            
            for k = 1:no_poss_neighbs
                
                if (all_neighb_indices(i,k)~=0) && all_neighb_indices(i,k)<=no_elem
                    
                    act_neighb = all_neighb_indices(i,k);
                    neighbs_act_neighb = all_neighb_indices(act_neighb,:);
                    ind_j_in_neighb = find(i==neighbs_act_neighb);
                    K_aniso_neighb = K_aniso(act_neighb,ind_j_in_neighb);
                    mob_a(i,k) = (k_rel_F_all(i,k)*corr_fact_orth_save(i,k)*A_f_all(i,k)* (((l_elem_all(i,k)*K_aniso(i,k) + l_neighb_all(i,k)*K_aniso_neighb))/(l_elem_all(i,k) + l_neighb_all(i,k))  )   ) /(l_elem_all(i,k) + l_neighb_all(i,k)) ;% It's all about mobility..
                    
                elseif (all_neighb_indices(i,k)~=0) && all_neighb_indices(i,k)>no_elem
                    mob_a(i,k) = k_rel_F_all(i,k)*(K_aniso(i,k)/l_elem_all(i,k))*A_f_all(i,k);
                    find_jN = find(i==elem_neumann_bound);
                    
                    if ~isempty(find_jN)
                        mob_a(i,k) = 1;
                    end
                    
                    new_counter = new_counter + 1;
                elseif (all_neighb_indices(i,k)==0)
                    mob_a(i,k) = 0;
                end
                
            end
            
            %% Robin mob-coeff:
            mob_robin = 0;
            
            robin_elems = find(i==elem_list_robin);
            
            if ~isempty(robin_elems)
                mob_coeff_elems_robin_i = mob_coeff_elems_robin(robin_elems);
                mob_robin = sum(mob_coeff_elems_robin_i);
            end
            
            %% Place them all..
            
            mob_matrix_diag(i) = mob_a(i,1) + mob_a(i,2) + mob_a(i,3) + mob_a(i,4) + mob_a(i,5) + mob_robin;
            
            for k = 1:no_poss_neighbs
                if (all_neighb_indices(i,k)~=0)
                    
                    row_mob_matrix_offdiag(off_diag_counter) = i;
                    col_mob_matrix_offdiag(off_diag_counter) = all_neighb_indices(i,k);
                    mob_matrix_offdiag(off_diag_counter) = - mob_a(i,k);
                    
                    off_diag_counter = off_diag_counter + 1;
                end
            end
        end
        
        
        if i>no_elem
            
            find_jD_place = find(all_bounds_DiriNeum(Bound_counter_place)==elem_diri_bound);
            find_jN_place = find(all_bounds_DiriNeum(Bound_counter_place)==elem_neumann_bound);
            
            if ~isempty(find_jD_place) && isempty(find_jN_place)
                
                row_mob_matrix_diri(diri_mat_counter) = i;
                col_mob_matrix_diri(diri_mat_counter) = i;
                
                starting_vec_solver(i) = vector_of_D_val_bounds(diri_mat_counter); % bicgstab-Aenderung
                
                diri_mat_counter = diri_mat_counter + 1;
                
            elseif ~isempty(find_jN_place)
                
                Neumann_counter_place = Neumann_counter_place + 1;
                
                row_mob_matrix_neum_diag(neum_mat_counter) = i;
                col_mob_matrix_neum_diag(neum_mat_counter) = i;
                row_mob_matrix_neum_offdiag(neum_mat_counter) = i;
                col_mob_matrix_neum_offdiag(neum_mat_counter) = elem_neumann_bound(Neumann_counter_place);
                
                neum_mat_counter = neum_mat_counter + 1;
                
                starting_vec_solver(i) = 1;
                
                
            end
            Bound_counter_place = Bound_counter_place  + 1;
            
            
        end
        
    end
    
    row_mob_matrix_diag = [1:1:no_elem]';
    col_mob_matrix_diag = [1:1:no_elem]';
    
    row_mob_matrix_offdiag(isnan(row_mob_matrix_offdiag)) = [];
    col_mob_matrix_offdiag(isnan(col_mob_matrix_offdiag)) = [];
    mob_matrix_offdiag(isnan(mob_matrix_offdiag)) = [];
    
    mob_matrix_diri = ones(no_Diri_elem,1);
    
    mob_matrix_neum_diag = ones(no_Neumann_elem,1);
    
    mob_matrix_neum_offdiag = ones(no_Neumann_elem,1)*-1;
    
    all_row = [row_mob_matrix_diag;row_mob_matrix_offdiag;row_mob_matrix_diri;row_mob_matrix_neum_diag;row_mob_matrix_neum_offdiag];
    all_col = [col_mob_matrix_diag;col_mob_matrix_offdiag;col_mob_matrix_diri;col_mob_matrix_neum_diag;col_mob_matrix_neum_offdiag];
    all_val = [mob_matrix_diag;mob_matrix_offdiag;mob_matrix_diri;mob_matrix_neum_diag;mob_matrix_neum_offdiag];
    
    B_mob_a = sparse(all_row,all_col,all_val);
    
    disp('Matrix assembly done. Wait for solving the system..')
    
    %% SOLVE ME: DETERMINE THE HEADS AT THE CENTROIDS
    
    if transient
        add_zeros_no = length(rhs_vector_FVM) - length(lump_transflow_inrhs);
        rhs_vector_FVM = rhs_vector_FVM + [lump_transflow_inrhs;zeros(add_zeros_no,1)];
    end
    
    [L,U] = ilu(B_mob_a);
    heads_centro_FVM_v = bicgstab(B_mob_a,rhs_vector_FVM,tol,maxit,L,U,starting_vec_solver);
    heads_centro_FVM = heads_centro_FVM_v;
    heads_centro_FVM(no_elem+1:end) = [];
    
    disp('System of equations solved!')
    %% Recover fluxes on faces
    Q_FV_fluxes = zeros(no_elem,5);
    Q_Neumann_int = 0;
    for i = 1:no_elem
        
        for k = 1:no_poss_neighbs
            if (all_neighb_indices(i,k)~=0)
                Q_FV_fluxes(i,k) = mob_a(i,k)*(heads_centro_FVM_v(i) - heads_centro_FVM_v(all_neighb_indices(i,k)));
            end
            
        end
        
        is_elem_robin = find(i==elem_list_robin);
        is_i_recharge_elem = find(i==inds_recharge);
        
        elms_lat_robin = elem_inds_leak_all_lat{:};
        is_elem_lat_robin = find(i==elms_lat_robin);
        
        ref_head_robin = ref_heads_leak_all_lat_vec{1};
        k_leak_lat_robin = k_leak1_lat_vec{1};
        dist_leak_lat_robin = dist_leak_lat_robin_vec{1};
        
        % Look at this in more detail later: this part still does not cover all cases!!!!
        
        if i>=(no_layers-1)*elem_per_layer+1 && surface_domain
            Q_FV_fluxes(i,5) = -Q_inex_olf_pm(i);
        end
        
        if i>=(no_layers-1)*elem_per_layer+1 && isempty(is_elem_lat_robin) && ~isempty(is_elem_robin)
            Q_FV_fluxes(i,5) = Q_FV_fluxes(i,5) - (sum(Q_FV_fluxes(i,:)));
            
            
        elseif i>=(no_layers-1)*elem_per_layer+1 && isempty(is_elem_lat_robin) && isempty(is_elem_robin)
            Q_FV_fluxes(i,5) = Q_FV_fluxes(i,5) - only_recharge(i);
            
        elseif i>=(no_layers-1)*elem_per_layer+1 && ~isempty(is_elem_lat_robin)
            % 1st face: 2,3,5,6
            % 2nd face: 1,3,4,6
            % 3rd face: 1,2,4,5
            
            % all lat Robin-BC nodes: n_inds_leak1_lat
            node_id1 = node_ids(i,1);
            node_id2 = node_ids(i,2);
            node_id3 = node_ids(i,3);
            node_id4 = node_ids(i,4);
            node_id5 = node_ids(i,5);
            node_id6 = node_ids(i,6);
            
            find_n1 = find(node_id1==n_inds_leak1_lat);
            find_n2 = find(node_id2==n_inds_leak1_lat);
            find_n3 = find(node_id3==n_inds_leak1_lat);
            find_n4 = find(node_id4==n_inds_leak1_lat);
            find_n5 = find(node_id5==n_inds_leak1_lat);
            find_n6 = find(node_id6==n_inds_leak1_lat);
            
            if ~isempty(find_n2) && ~isempty(find_n3) && ~isempty(find_n5) && ~isempty(find_n6)
                
                Q_FV_fluxes(i,5) = Q_FV_fluxes(i,5) - only_recharge(i);
                Q_FV_fluxes(i,1) = - (sum(Q_FV_fluxes(i,:)));
                
            elseif ~isempty(find_n1) && ~isempty(find_n3) && ~isempty(find_n4) && ~isempty(find_n6)
                
                Q_FV_fluxes(i,2) = Q_FV_fluxes(i,5) - only_recharge(i);
                Q_FV_fluxes(i,5) = - (sum(Q_FV_fluxes(i,:)));
            elseif ~isempty(find_n1) && ~isempty(find_n2) && ~isempty(find_n4) && ~isempty(find_n5)
                
                Q_FV_fluxes(i,3) = Q_FV_fluxes(i,5) - only_recharge(i);
                Q_FV_fluxes(i,5) = - (sum(Q_FV_fluxes(i,:)));
            end
            
            
        elseif ~(i>=(no_layers-1)*elem_per_layer+1) && ~isempty(is_elem_lat_robin)
            % 1st face: 2,3,5,6
            % 2nd face: 1,3,4,6
            % 3rd face: 1,2,4,5
            
            % all lat Robin-BC nodes: n_inds_leak1_lat
            node_id1 = node_ids(i,1);
            node_id2 = node_ids(i,2);
            node_id3 = node_ids(i,3);
            node_id4 = node_ids(i,4);
            node_id5 = node_ids(i,5);
            node_id6 = node_ids(i,6);
            
            find_n1 = find(node_id1==n_inds_leak1_lat);
            find_n2 = find(node_id2==n_inds_leak1_lat);
            find_n3 = find(node_id3==n_inds_leak1_lat);
            find_n4 = find(node_id4==n_inds_leak1_lat);
            find_n5 = find(node_id5==n_inds_leak1_lat);
            find_n6 = find(node_id6==n_inds_leak1_lat);
            
            if ~isempty(find_n2) && ~isempty(find_n3) && ~isempty(find_n5) && ~isempty(find_n6)
                
                Q_FV_fluxes(i,1) = - sum(Q_FV_fluxes(i,:));
                
            elseif ~isempty(find_n1) && ~isempty(find_n3) && ~isempty(find_n4) && ~isempty(find_n6)
                
                Q_FV_fluxes(i,2) = - sum(Q_FV_fluxes(i,:));
                
            elseif ~isempty(find_n1) && ~isempty(find_n2) && ~isempty(find_n4) && ~isempty(find_n5)
                
                Q_FV_fluxes(i,3) = - sum(Q_FV_fluxes(i,:));
            end
            
            
            
        end
        
        
    end
    
    disp('Fluxes are recovered!')
    %% Save fluxes:
    % mass conservation test:
    mass_cons = sum(Q_FV_fluxes,2);
    
    q_FV_fluxes = Q_FV_fluxes./A_f_all;
    q_FV_fluxes(isnan(q_FV_fluxes)) = 0;
    
    Q_div_zero_name = strcat('Q_div_zero_',prefix,'.',target_time_step_suffix,'.mat');
    save(Q_div_zero_name,'Q_FV_fluxes')
    
    q_locn_name = strcat('q_locn_',prefix,'.',target_time_step_suffix,'.mat');
    save(q_locn_name,'q_FV_fluxes')
    
    
    %% VELOCITY FIELD COMPUTATIONS
    
    Q_div0 = Q_FV_fluxes;
    q_locn = q_FV_fluxes;
    
    hQ_all = importdata(hQ_all_name);
    
    a_x = zeros(no_elem,1);
    a_y = zeros(no_elem,1);
    a_z = zeros(no_elem,1);
    b_div_horiz = zeros(no_elem,1);
    b_div_vert = zeros(no_elem,1);
    
    a_zn = zeros(no_elem,1);
    b_div_vertn = zeros(no_elem,1);
    
    detJ_centro = zeros(no_elem,1);
    
    for j = 1:no_elem
        
        node_ids1 = node_ids(j,1);
        node_ids2 = node_ids(j,2);
        node_ids3 = node_ids(j,3);
        node_ids4 = node_ids(j,4);
        node_ids5 = node_ids(j,5);
        node_ids6 = node_ids(j,6);
        
        xA = coordinates(node_ids1,1);
        yA = coordinates(node_ids1,2);
        zA = coordinates(node_ids1,3);
        
        xB = coordinates(node_ids2,1);
        yB = coordinates(node_ids2,2);
        zB = coordinates(node_ids2,3);
        
        xC = coordinates(node_ids3,1);
        yC = coordinates(node_ids3,2);
        zC = coordinates(node_ids3,3);
        
        xD = coordinates(node_ids4,1);
        yD = coordinates(node_ids4,2);
        zD = coordinates(node_ids4,3);
        
        xE = coordinates(node_ids5,1);
        yE = coordinates(node_ids5,2);
        zE = coordinates(node_ids5,3);
        
        xF = coordinates(node_ids6,1);
        yF = coordinates(node_ids6,2);
        zF = coordinates(node_ids6,3);
        
        Xpri = [xA yA zA;xB yB zB;xC yC zC;xD yD zD;xE yE zE;xF yF zF];
        
        
        Qpri = Q_div0(j,1:3);
        hQ = hQ_all(j,:);
        hQ(hQ==0) = 1; % heuristic, for no-flow boundary conditions, be aware in case of code changes that all boundaries are covered
        
        a_x(j) = (-1/(2*A_horiz(j)))*( (Qpri(1)/hQ(1))*Xpri(1,1) + (Qpri(2)/hQ(2))*Xpri(2,1) + (Qpri(3)/hQ(3))*Xpri(3,1));
        a_y(j) = (-1/(2*A_horiz(j)))*( (Qpri(1)/hQ(1))*Xpri(1,2) + (Qpri(2)/hQ(2))*Xpri(2,2) + (Qpri(3)/hQ(3))*Xpri(3,2));
        b_div_horiz(j) = (1/(2*A_horiz(j)))*((Qpri(1)/hQ(1)) + (Qpri(2)/hQ(2)) + (Qpri(3)/hQ(3)));
        
        xloc_centro = [(xA + xB + xC)/3 ,(yA + yB + yC)/3, 0.5];
        A_n = A_horiz(j);
        
        gradN_locmix = [ ((yB-yC)/(2*A_n))*((1-xloc_centro(3))), ((yC-yA)/(2*A_n))*(1-xloc_centro(3)), ((yA-yB)/(2*A_n))*(1-xloc_centro(3)), ((yB-yC)/(2*A_n))* (xloc_centro(3)), ((yC-yA)/(2*A_n))* (xloc_centro(3)), ((yA-yB)/(2*A_n))* (xloc_centro(3));...
            ((xC-xB)/(2*A_n))*(1-xloc_centro(3)), ((xA-xC)/(2*A_n))*(1-xloc_centro(3)),((xB-xA)/(2*A_n))*(1-xloc_centro(3)), ((xC-xB)/(2*A_n))* (xloc_centro(3)), ((xA-xC)/(2*A_n))* (xloc_centro(3)),((xB-xA)/(2*A_n))* (xloc_centro(3));...
            -1*((xB*yC - xC*yB)/(2*A_n) + ((yB-yC)/(2*A_n))*xloc_centro(1) + ((xC-xB)/(2*A_n))*xloc_centro(2)), -1*((xC*yA - xA*yC)/(2*A_n) + ((yC-yA)/(2*A_n))*xloc_centro(1) + ((xA-xC)/(2*A_n))*xloc_centro(2)),...
            -1*((xA*yB - xB*yA)/(2*A_n) + ((yA-yB)/(2*A_n))*xloc_centro(1) + ((xB-xA)/(2*A_n))*xloc_centro(2)),((xB*yC - xC*yB)/(2*A_n) + ((yB-yC)/(2*A_n))*xloc_centro(1) + ((xC-xB)/(2*A_n))*xloc_centro(2))...
            ,((xC*yA - xA*yC)/(2*A_n) + ((yC-yA)/(2*A_n))*xloc_centro(1) + ((xA-xC)/(2*A_n))*xloc_centro(2)),((xA*yB - xB*yA)/(2*A_n) + ((yA-yB)/(2*A_n))*xloc_centro(1) + ((xB-xA)/(2*A_n))*xloc_centro(2))];
        
        J_loc_centro = gradN_locmix*Xpri;
        
        % making use of mixed coordinates:
        detJ_centro(j) = J_loc_centro(3,3);
        
        scaling_fac_bottom = A_f_all(j,4)/A_n;
        scaling_fac_top = A_f_all(j,5)/A_n;
        
        a_z(j) = (- q_locn(j,4)*scaling_fac_bottom/detJ_centro(j));
        b_div_vert(j) = (((q_locn(j,4)*scaling_fac_bottom) + (q_locn(j,5))* scaling_fac_top)/detJ_centro(j));
        
        
    end
    
    a_x_name = strcat('a_x_',prefix,'.',target_time_step_suffix,'.mat');
    save(a_x_name,'a_x')
    
    a_y_name = strcat('a_y_',prefix,'.',target_time_step_suffix,'.mat');
    save(a_y_name,'a_y')
    
    a_z_name = strcat('a_z_',prefix,'.',target_time_step_suffix,'.mat');
    save(a_z_name,'a_z')
    
    b_div_horiz_name = strcat('b_div_horiz_',prefix,'.',target_time_step_suffix,'.mat');
    save(b_div_horiz_name,'b_div_horiz')
    
    b_div_vert_name = strcat('b_div_vert_',prefix,'.',target_time_step_suffix,'.mat');
    save(b_div_vert_name,'b_div_vert')
    
    theta_actual_name = strcat('theta_',prefix,'.',target_time_step_suffix,'.mat');
    save(theta_actual_name,'mean_theta')
    
    disp('Attributes of RTN_0 velocity field are computed!')
end

detJ_centro_name = strcat('detJ_centro_',prefix,'.mat');
save(detJ_centro_name,'detJ_centro')

n_name = strcat('n_',prefix,'.mat');
save(n_name,'n')

zone_for_elements_name = strcat('zone_for_elements_',prefix,'.mat');
save(zone_for_elements_name,'zone_for_element')

%% Plotting:
toc
if plotting
    
    down = 0; % TO ADAPT LINE, go down the layers: down = 0 == top layer, down = +1, one layer beneath, and then iterate it: +=1
    head_top = heads_centro_FVM((no_layers-1-down)*elem_per_layer+1:(no_layers-down)*elem_per_layer);
    
    heads_max = heads_centro_FVM((no_layers-1-down)*elem_per_layer+1:(no_layers-down)*elem_per_layer);
    
    inds_crit = (no_layers-1-down)*elem_per_layer+1:(no_layers-down)*elem_per_layer;
    
    elms_crit = (no_layers-1-down)*elem_per_layer+1:(no_layers-down)*elem_per_layer;
    resp_inds = zeros(length(elms_crit),1);
    
    n_grid = 20;
    min_x_grid = floor(min(coordinates(:,1)));
    max_x_grid = ceil(max(coordinates(:,1)));
    min_y_grid = floor(min(coordinates(:,2)));
    max_y_grid = ceil(max(coordinates(:,2)));
    x_grid_step = (max_x_grid - min_x_grid)/n_grid;
    y_grid_step = (max_y_grid - min_y_grid)/n_grid;
    
    [xq,yq]=meshgrid(min_x_grid:x_grid_step:max_x_grid,min_y_grid:y_grid_step:max_y_grid);
    Vq = griddata(centro_of_elem(1:elem_per_layer,1),centro_of_elem(1:elem_per_layer,2),full(head_top),xq,yq);
    
    fig1=figure(1);
    hold on
    
    for i = 1:elem_per_layer
        
        node_ids_sh1 = node_ids(i,1);
        node_ids_sh2 = node_ids(i,2);
        node_ids_sh3 = node_ids(i,3);
        
        x1 = coordinates(node_ids_sh1,1);
        y1 = coordinates(node_ids_sh1,2);
        x2 = coordinates(node_ids_sh2,1);
        y2 = coordinates(node_ids_sh2,2);
        x3 = coordinates(node_ids_sh3,1);
        y3 = coordinates(node_ids_sh3,2);
        
        
        
        line([x1 x2 x3 x1],[y1 y2 y3 y1],'Color',[0.6 0.6 0.6])
    end
    
    xlabel('x [m]');ylabel('y [m]');
    
    frame_plot_x = (max(coordinates(:,1)) - min(coordinates(:,1)))/20;
    frame_plot_y = (max(coordinates(:,2)) - min(coordinates(:,2)))/20;
    
    axis([min(coordinates(:,1))-frame_plot_x max(coordinates(:,1))+frame_plot_x min(coordinates(:,2))-frame_plot_y max(coordinates(:,2))+frame_plot_y]);
    daspect([1 1 1]);
    
    contour(xq,yq,Vq,'-b','LineWidth',1.5,'ShowText','on');
    daspect([1 1 1]);
    
    set(gca,'layer','top');
    set(gcf,'PaperPositionMode','auto')
    
    set(gca,'FontSize',12)
    
    Pos = get(gcf,'Position');
    set(gcf,'Position',[Pos(1),Pos(2),4*156,3*156])
    
    hold off
    title('Pseudopotential: hydraulic head [m]')
end

% print(fig1,'a_floodplain_pot','-dpng','-r300');

disp('ready')



