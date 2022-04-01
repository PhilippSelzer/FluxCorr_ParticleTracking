%% Metainformation:
%
% Author:           Philipp Selzer
% Version:          2.0
% Date:             20.11.2019
% Last update:      24.03.2022
% Purpose:          This code is the actual particle tracking scheme
%                   applicable to variably saturated porous media flow on
%                   deformed triangular prisms.

% Correspondance:   philipp.selzer@gmx.net
% Copyright(c) 2019: Philipp Selzer
% License: An extended version of the LGPLv3 (http://www.gnu.org/licenses/lgpl-3.0.en.html), see
% attached license-file

% Acknowledgement: The author thanks Jonas Allgeier and Daniel Erdal for
% their valuable comments and suggestions, as well as for testing my codes
% and their help in the debugging process. The author thanks his PhD
% supervisors Olaf A. Cirpka and René Therrien for their valuable 
% contributions and their advise.

%% Parallel Particle Tracking on CPU or GPU
%%
clear all
close all
clc

format long

prefix = 'floodplain_model'; % TO ADAPT
plotting = true;

reverse_particle_tracking = false;

crit_div_horiz = 1e-19; % Switch between divergent and non-divergent flow in the horizontal, equal to critical divergence term (i.e., "b") in RTN_0 velocity space
crit_div_vert = 1e-19; % Switch between divergent and non-divergent flow in the vertical, equal to critical divergence term (i.e., "b") in RTN_0 velocity space

%% GPU-flag

on_gpu = false; % compute particle tracking on the GPU
gpu_slim = false; % can only be true, if on_gpu == true, this flag puts
% only very little data on the gpu

%% ===================================================================
% Read all relevant HGS output files
%% ===================================================================

% read the binary node ids of the elements
n_dim = 3; % number of dimensions (not to alter)
filename1 = strcat(prefix,'o.elements_pm');

no_nodes_per_element = 6;

no_elem =138565;% TO ADAPT
no_layers = 35;% TO ADAPT
elem_per_layer = no_elem/no_layers;

% read the binary elements
[node_ids] = read_HGS_binary_elements(no_nodes_per_element,no_elem,filename1);

node_ids = uint32(node_ids);

% read the binary coordinates
filename2 = strcat(prefix,'o.coordinates_pm');
xyz = 3;
no_nodes = 74412;% TO ADAPT

[coordinates] = read_HGS_binary_coordinates(xyz,no_nodes,filename2);

preallocation_no_iterations = round(elem_per_layer/10);

%% extraction/injection wells
wells_exist = false;

if wells_exist
    ele_list_wells_name = strcat('ele_list_wells_',prefix,'.mat');
    ele_list_wells = importdata(ele_list_wells_name);
    
    ele_list_wells_all = NaN(no_elem,1);
    counter_w = 1;
    for ww = 1:length(ele_list_wells)
        ele_list_wells_act = ele_list_wells{ww};
        
        ele_list_wells_all(counter_w:counter_w+ length(ele_list_wells_act)-1) = ele_list_wells_act;
        
        counter_w = counter_w + length(ele_list_wells_act);
        
    end
    
    ele_list_wells_all(isnan(ele_list_wells_all)) = [];
    
    if on_gpu
        ele_list_wells_all = gpuarray(ele_list_wells_all);
    end
    
end
%% Initial particle locations

starting_elements_are_known = true; % Starting elements are equal to elem_ID_parts later
% to be defined, if true

if ~starting_elements_are_known
    tetrahedral_nodes_name = strcat('tetrahedral_nodes_',prefix,'.mat');
    tetrahedral_nodes = importdata(tetrahedral_nodes_name);
end

% Example for starting locations, define them as needed. You may also not define
% element IDs by starting_elements_are_known = false, in this case the
% elements IDs are searched vie the input coordinates within the grid.
% However, you still need to give initial coordinates for the particles.

elem_ID_parts = no_elem-elem_per_layer+1:no_elem;

elem_ID_parts_vv = repmat(elem_ID_parts,[1,1]);
elem_ID_parts = elem_ID_parts_vv;
no_particles = length(elem_ID_parts);

vec_of_part_N = 1:1:no_particles;

if on_gpu   
    vec_of_part_N = gpuArray(vec_of_part_N );
end

centro_of_elem_name = strcat('centro_of_elem_',prefix,'.mat');
centro_of_elem = importdata(centro_of_elem_name);

xloc_part = centro_of_elem(elem_ID_parts,1);
yloc_part = centro_of_elem(elem_ID_parts,2);
zloc_part_glob = centro_of_elem(elem_ID_parts,3);

max_simulation_time= Inf; % default should be: Inf

%% Load needed variables

transient = true; % if transient == false FVM_corr_target_time_IDs should be FVM_corr_target_time_IDs=1;
% and target_times can be any value, e.g. "Inf"

% adaptation of strings

FVM_corr_target_time_IDs = [2,3,4,5];% Minimum value should be 2, and minmum two values, if transient == true/ for steady-state simulations this should stay "1"
% FVM_corr_target_time_IDs = [2]; % steady state

no_of_tracking_time_steps = length(FVM_corr_target_time_IDs);
target_times =[3e7,3e8,3e9,3e11,3e12];

steady_state_time = target_times(end);

if transient
    previous_time  = 0;
    next_time_update = target_times(FVM_corr_target_time_IDs(1));
end

target_time_step = FVM_corr_target_time_IDs(1);

if on_gpu
    FVM_corr_target_time_IDs = gpuArray(FVM_corr_target_time_IDs);
    target_times = gpuArray(target_times);
    steady_state_time = gpuArray(steady_state_time);
end

target_time_step_suffix = strcat('000',num2str(target_time_step));

if ~transient
    
    q_locn_name = strcat('q_locn_',prefix,'.',target_time_step_suffix,'.mat');
    q_FV_fluxes = importdata(q_locn_name);
    q_locn = q_FV_fluxes;
    
    a_x_name = strcat('a_x_',prefix,'.',target_time_step_suffix,'.mat');
    a_y_name = strcat('a_y_',prefix,'.',target_time_step_suffix,'.mat');
    a_z_name = strcat('a_z_',prefix,'.',target_time_step_suffix,'.mat');
    b_div_horiz_name = strcat('b_div_horiz_',prefix,'.',target_time_step_suffix,'.mat');
    b_div_vert_name = strcat('b_div_vert_',prefix,'.',target_time_step_suffix,'.mat');
    theta_name = strcat('theta_',prefix,'.',target_time_step_suffix,'.mat');
    
    a_x = importdata(a_x_name);
    a_y = importdata(a_y_name);
    a_z = importdata(a_z_name);
    b_div_horiz = importdata(b_div_horiz_name);
    b_div_vert = importdata(b_div_vert_name);
    theta = importdata(theta_name);
    
    a_x = a_x./theta;
    a_y = a_y./theta;
    a_z = a_z./theta;
    b_div_horiz = b_div_horiz./theta;
    b_div_vert = b_div_vert./theta;
    
    if reverse_particle_tracking
        q_locn = -q_locn;
        a_x = -a_x;
        a_y = -a_y;
        a_z = -a_z;
        b_div_horiz = -b_div_horiz;
        b_div_vert = -b_div_vert;
    end
    
    a_xmat = a_x;
    a_ymat = a_y;
    a_zmat = a_z;
    b_div_horiz_mat = b_div_horiz;
    b_div_vert_mat = b_div_vert;
    q_locn_mat = q_locn;
    
else
    
    v_field_transient_counter = 1;
    
    a_xmat = NaN(no_elem,no_of_tracking_time_steps);
    a_ymat = NaN(no_elem,no_of_tracking_time_steps);
    a_zmat = NaN(no_elem,no_of_tracking_time_steps);
    b_div_horiz_mat = NaN(no_elem,no_of_tracking_time_steps);
    b_div_vert_mat = NaN(no_elem,no_of_tracking_time_steps);
    q_locn_mat = NaN(no_elem,5,no_of_tracking_time_steps);
    
    for target_time_step = FVM_corr_target_time_IDs
        
        % target_time_step suffix
        
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
        
        target_time_step_suffix = strcat(target_string_no,num2str(target_time_step));
        
        q_locn_name = strcat('q_locn_',prefix,'.',target_time_step_suffix,'.mat');
        q_FV_fluxes = importdata(q_locn_name);
        
        a_x_name = strcat('a_x_',prefix,'.',target_time_step_suffix,'.mat');
        a_y_name = strcat('a_y_',prefix,'.',target_time_step_suffix,'.mat');
        a_z_name = strcat('a_z_',prefix,'.',target_time_step_suffix,'.mat');
        b_div_horiz_name = strcat('b_div_horiz_',prefix,'.',target_time_step_suffix,'.mat');
        b_div_vert_name = strcat('b_div_vert_',prefix,'.',target_time_step_suffix,'.mat');
        theta_name = strcat('theta_',prefix,'.',target_time_step_suffix,'.mat');
        
        a_xv = importdata(a_x_name);
        a_yv = importdata(a_y_name);
        a_zv = importdata(a_z_name);
        b_div_horizv = importdata(b_div_horiz_name);
        b_div_vertv = importdata(b_div_vert_name);
        thetav = importdata(theta_name);
        
        a_xv = a_xv./thetav;
        a_yv = a_yv./thetav;
        a_zv = a_zv./thetav;
        b_div_horizv = b_div_horizv./thetav;
        b_div_vertv = b_div_vertv./thetav;
        
        if reverse_particle_tracking
            a_xv = -a_xv;
            a_yv = -a_yv;
            a_zv = -a_zv;
            b_div_horizv = -b_div_horizv;
            b_div_vertv = -b_div_vertv;
            q_FV_fluxes = -q_FV_fluxes;
        end
        
        a_xmat(:,v_field_transient_counter) = a_xv;
        a_ymat(:,v_field_transient_counter) = a_yv;
        a_zmat(:,v_field_transient_counter) = a_zv;
        b_div_horiz_mat(:,v_field_transient_counter) = b_div_horizv;
        b_div_vert_mat(:,v_field_transient_counter) = b_div_vertv;
        q_locn_mat(:,:,v_field_transient_counter) = q_FV_fluxes;
        
        v_field_transient_counter = v_field_transient_counter + 1;
        
    end
    
    q_locn =  q_locn_mat(:,:,1);
    a_x = a_xmat(:,1);
    a_y = a_ymat(:,1);
    a_z = a_zmat(:,1);
    b_div_horiz = b_div_horiz_mat(:,1);
    b_div_vert =  b_div_vert_mat(:,1);
end

if on_gpu && ~gpu_slim
    q_locn_mat = gpuArray(q_locn_mat);
    a_xmat = gpuArray(a_xmat);
    a_ymat = gpuArray(a_ymat);
    a_zmat = gpuArray(a_zmat);
    b_div_horiz_mat = gpuArray(b_div_horiz_mat);
    b_div_vert_mat = gpuArray(b_div_vert_mat);
end

A_horiz_name = strcat('A_horiz_',prefix,'.mat');
A_horiz = importdata(A_horiz_name);

all_neighb_indices_name = strcat('all_neighb_indices_',prefix,'.mat');
all_neighb_indices = importdata(all_neighb_indices_name);
all_neighb_indices(all_neighb_indices>no_elem) = 0;

hQ_all_name = strcat('hQ_all_',prefix,'.mat');
hQ_all = importdata(hQ_all_name);

n_name = strcat('n_',prefix,'.mat');
zone_for_elements_name = strcat('zone_for_elements_',prefix,'.mat');

n = importdata(n_name);
zone_for_elements = importdata(zone_for_elements_name);

if on_gpu
    if ~gpu_slim
        A_horiz = gpuArray(A_horiz);
        all_neighb_indices = gpuArray(all_neighb_indices);
        hQ_all = gpuArray(hQ_all);
        n = gpuArray(n);
    end
    zone_for_elements = gpuArray(zone_for_elements);
end

%% Search algorithm for initial element IDs of particle locations:

% 1. Split every prism in three tetrahedra
% 2. Search in every tetrahedron for the coordinates by baraycentric coordinates
% 3. Determine via the ID of the tetrahedron inluding the coordinates the ID of th prism

starting_ele_part= zeros(1,no_particles);
if on_gpu
    starting_ele_part= gpuArray(starting_ele_part);
end

if ~starting_elements_are_known
    
    lambda = nan(4,no_particles);
    lambda_tri = nan(3,no_particles);
    lambda_pos = nan(4,no_particles);
    
    found_in_ele = false(no_particles,1);
    rt_tri = zeros(2,no_particles);
    rt = zeros(3,no_particles);
    
    
    
    while_iter = 0;
    while sum(found_in_ele)<no_particles
        
        % Grid needs to be stacked for this!!! (Same as for algorithm determining element neighbours)
        
        for kk = 1:elem_per_layer % search in triangular base areas, brute force, needs to be replaced by something smart
            
            node_ids1 = node_ids(kk,1);
            node_ids2 = node_ids(kk,2);
            node_ids3 = node_ids(kk,3);
            
            xA = coordinates(node_ids1,1);
            yA = coordinates(node_ids1,2);
            
            xB = coordinates(node_ids2,1);
            yB = coordinates(node_ids2,2);
            
            xC = coordinates(node_ids3,1);
            yC = coordinates(node_ids3,2);
            
            rt_tri(:,vec_of_part_N(~found_in_ele)) = [xloc_part(vec_of_part_N(~found_in_ele))-xC,yloc_part(vec_of_part_N(~found_in_ele))-yC]';
            
            At_tri = [(xA-xC) (xB-xC); (yA-yC) (yB-yC)];
            
            at = At_tri(1,1);
            bt = At_tri(1,2);
            ct = At_tri(2,1);
            dt = At_tri(2,2);
            
            At_det_tri = At_tri(1,1)*At_tri(2,2) - At_tri(2,1)*At_tri(1,2);
            At_m_tri = [At_tri(2,2) -At_tri(1,2);-At_tri(2,1) At_tri(1,1)];
            At_inv_tri = (1/At_det_tri)*At_m_tri;
            lambda_tri(1:2,vec_of_part_N(~found_in_ele)) = At_inv_tri*rt_tri(:,vec_of_part_N(~found_in_ele));
            
            lambda_tri(3,vec_of_part_N(~found_in_ele)) = 1-lambda_tri(1,vec_of_part_N(~found_in_ele))-lambda_tri(2,vec_of_part_N(~found_in_ele));
            
            
            condition_of_truth = ((lambda_tri(1,vec_of_part_N(~found_in_ele))>=0 & lambda_tri(1,vec_of_part_N(~found_in_ele))<=1) & (lambda_tri(2,vec_of_part_N(~found_in_ele))>=0 & lambda_tri(2,vec_of_part_N(~found_in_ele))<=1) & (lambda_tri(3,vec_of_part_N(~found_in_ele))>=-50*eps & lambda_tri(3,vec_of_part_N(~found_in_ele))<=(1+50*eps)))';
            
            if sum(condition_of_truth)>0
                
                sum_condition_of_truth_tetra = 0;
                inds_still_available = vec_of_part_N(~found_in_ele);
                inds_of_parts_found_of_avail = find(condition_of_truth==1);
                inds_of_parts_found = inds_still_available(inds_of_parts_found_of_avail);
                
                list_of_elements_to_search = [ones(1,no_layers)*kk + elem_per_layer*[0:1:no_layers-1]];
                
                %% Algorithm about relevant elements
                
                for j = list_of_elements_to_search
                    tetras_of_prism = [tetrahedral_nodes((j-1)*3 + 1,:);tetrahedral_nodes((j-1)*3 + 2,:);tetrahedral_nodes((j-1)*3 + 3,:)];
                    for jj = 1:3
                        
                        current_tetra = tetras_of_prism(jj,:);
                        
                        xA = coordinates(current_tetra(1),1);
                        yA = coordinates(current_tetra(1),2);
                        zA = coordinates(current_tetra(1),3);
                        
                        xB = coordinates(current_tetra(2),1);
                        yB = coordinates(current_tetra(2),2);
                        zB = coordinates(current_tetra(2),3);
                        
                        xC = coordinates(current_tetra(3),1);
                        yC = coordinates(current_tetra(3),2);
                        zC = coordinates(current_tetra(3),3);
                        
                        xD = coordinates(current_tetra(4),1);
                        yD = coordinates(current_tetra(4),2);
                        zD = coordinates(current_tetra(4),3);
                        
                        rt(:,inds_of_parts_found) = [xloc_part(inds_of_parts_found)-xD,yloc_part(inds_of_parts_found)-yD,zloc_part_glob(inds_of_parts_found)-zD]';%[xloc_part(ip)-xD;yloc_part(ip)-yD;zloc_part_glob(ip)-zD];
                        At = [(xA-xD) (xB-xD) (xC-xD); (yA-yD) (yB-yD) (yC-yD);(zA-zD) (zB-zD) (zC-zD)];
                        
                        
                        at = At(1,1);
                        bt = At(1,2);
                        ct = At(1,3);
                        dt = At(2,1);
                        et = At(2,2);
                        ft = At(2,3);
                        gt = At(3,1);
                        ht = At(3,2);
                        it = At(3,3);
                        
                        At_det = (At(1,1)*At(2,2)*At(3,3))+ (At(1,2)*At(2,3)*At(3,1)) + (At(1,3)*At(2,1)*At(3,2))...
                            - (At(3,1)*At(2,2)*At(1,3)) -(At(2,1)*At(1,2)*At(3,3)) -(At(1,1)*At(3,2)*At(2,3));
                        At_m =[(et*it - ft*ht) -(bt*it - ct*ht) (bt*ft - ct*et); -(dt*it - ft*gt) (at*it - ct*gt) -(at*ft - ct*dt);...
                            (dt*ht - et*gt) -(at*ht - bt*gt) (at*et - bt*dt)];
                        At_inv = (1/At_det)*At_m;
                        
                        lambda(1:3,inds_of_parts_found) = At_inv*rt(:,inds_of_parts_found);
                        
                        lambda(4,inds_of_parts_found) = 1-lambda(1,inds_of_parts_found)-lambda(2,inds_of_parts_found)-lambda(3,inds_of_parts_found);
                        
                        condition_of_truth_tetra = (lambda(1,inds_of_parts_found)>=0 & lambda(1,inds_of_parts_found)<=1) & (lambda(2,inds_of_parts_found)>=0 & lambda(2,inds_of_parts_found)<=1) & (lambda(3,inds_of_parts_found)>=0 & lambda(3,inds_of_parts_found)<=1) & (lambda(4,inds_of_parts_found)>=-50*eps & lambda(4,inds_of_parts_found)<=(1+50*eps));
                        %
                        
                        if sum(condition_of_truth_tetra)>0
                            
                            ind_tetra_cond_truth_where = find(condition_of_truth_tetra==1);
                            
                            starting_ele_part(inds_of_parts_found(ind_tetra_cond_truth_where)) = j;
                            found_in_ele(inds_of_parts_found(ind_tetra_cond_truth_where)) = true;
                            sum_condition_of_truth_tetra = sum_condition_of_truth_tetra + sum(condition_of_truth_tetra);
                            %% Move particle within the element, if on node (not the most elegant way, but the most simple one,and the introduced inaccuracy is
                            % far below real-life measurement accuracy. Moreover, the movement is saved. This work-around is, therefore, made transparent.)
                            for ll = inds_of_parts_found
                                lambda_pos = find(lambda(:,ll)>0);
                                lambda_zero = setdiff([1,2,3,4],lambda_pos);
                                
                                if ((xloc_part(ll)==xA && yloc_part(ll)==yA) || (xloc_part(ll)==xB && yloc_part(ll)==yB) || (xloc_part(ll)==xC && yloc_part(ll)==yC) || (xloc_part(ll)==xD && yloc_part(ll)==yD)) || ~isempty(lambda_zero)
                                    if length(lambda_pos)==2
                                        lambda(lambda_pos) = lambda(lambda_pos) - 5000*eps;
                                        lambda(lambda_zero) = 5000*eps;
                                        xloc_part(ll) = [xA,xB,xC,xD]*lambda';
                                        yloc_part(ll) = [yA,yB,yC,yD]*lambda';
                                        zloc_part_glob(ll) = [zA,zB,zC,zD]*lambda';
                                    elseif length(lambda_pos)==1
                                        lambda(lambda_pos) = lambda(lambda_pos) - 6000*eps;
                                        lambda(lambda_zero) = 2000*eps;
                                        xloc_part(ll) = [xA,xB,xC,xD]*lambda';
                                        yloc_part(ll) = [yA,yB,yC,yD]*lambda';
                                        zloc_part_glob(ll) = [zA,zB,zC,zD]*lambda';
                                    end
                                end
                            end
                            
                            
                            
                        end
                        
                        
                        if sum_condition_of_truth_tetra==sum(condition_of_truth)
                            break
                        end
                        
                    end % for tetra
                    
                    if sum_condition_of_truth_tetra==sum(condition_of_truth)
                        break % vertical elements
                    end
                    
                end % for list of elements
                
                if sum(found_in_ele)==no_particles
                    break % triangles
                end
            end % if in triangle
            
            while_iter = while_iter + 1;
            if while_iter==elem_per_layer % emergency break
                break
            end
            
        end % for in triangles
        
    end
else
    
    starting_ele_part = elem_ID_parts;
    
end

%% Transform the global z-coordinate in local coordinates

z_glob = zloc_part_glob;

node_ids1 = node_ids(starting_ele_part',1);
node_ids2 = node_ids(starting_ele_part',2);
node_ids3 = node_ids(starting_ele_part',3);
node_ids4 = node_ids(starting_ele_part',4);
node_ids5 = node_ids(starting_ele_part',5);
node_ids6 = node_ids(starting_ele_part',6);

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

l_pos = [xloc_part';yloc_part';z_glob'];
l_dir_norm_l = [0;0;-1];
l_dir_norm_u = [0;0;1];

% lower intersection point
p_pos_x_l = (xA + xB + xC)/3;
p_pos_y_l = (yA + yB + yC)/3;
p_pos_z_l = (zA + zB + zC)/3;

p_pos_l = [p_pos_x_l';p_pos_y_l';p_pos_z_l'];

p_dir1_l = [(xB - xA)';(yB - yA)';(zB - zA)'];
p_dir2_l = [(xC - xA)';(yC - yA)';(zC - zA)'];

n_face_l = [(p_dir1_l(2,:).*p_dir2_l(3,:)) - (p_dir1_l(3,:).*p_dir2_l(2,:)); (p_dir1_l(3,:).*p_dir2_l(1,:)) - (p_dir1_l(1,:).*p_dir2_l(3,:)); ...
    (p_dir1_l(1,:).*p_dir2_l(2,:)) - (p_dir1_l(2,:).*p_dir2_l(1,:))];

n_face_abs_l_sq = n_face_l.^2;
n_face_abs_l = sqrt(sum(n_face_abs_l_sq,1));
n_face_l = n_face_l./n_face_abs_l;

dist_position_vec_l = p_pos_l - l_pos;

d_line_l = (dist_position_vec_l(1,:).*n_face_l(1,:) + dist_position_vec_l(2,:).*n_face_l(2,:) + dist_position_vec_l(3,:).*n_face_l(3,:))...
    ./(l_dir_norm_l(1,:).*n_face_l(1,:) + l_dir_norm_l(2,:).*n_face_l(2,:) + l_dir_norm_l(3,:).*n_face_l(3,:));

P_intsec_l = l_pos + ((d_line_l')*(l_dir_norm_l'))';
P_intsec_l_z = P_intsec_l(3,:);

% upper intersection point
p_pos_x_u = (xD + xE + xF)/3;
p_pos_y_u = (yD + yE + yF)/3;
p_pos_z_u = (zD + zE + zF)/3;

p_pos_u = [p_pos_x_u';p_pos_y_u';p_pos_z_u'];

p_dir1_u = [(xE - xD)';(yE - yD)';(zE - zD)'];
p_dir2_u = [(xF - xD)';(yF - yA)';(zF - zD)'];

n_face_u = [(p_dir1_u(2,:).*p_dir2_u(3,:)) - (p_dir1_u(3,:).*p_dir2_u(2,:)); (p_dir1_u(3,:).*p_dir2_u(1,:)) - (p_dir1_u(1,:).*p_dir2_u(3,:)); ...
    (p_dir1_u(1,:).*p_dir2_u(2,:)) - (p_dir1_u(2,:).*p_dir2_u(1,:))];

n_face_abs_l_sq = n_face_u.^2;
n_face_abs_u = sqrt(sum(n_face_abs_l_sq,1));
n_face_u = n_face_u./n_face_abs_u;

% Compute actual intersection point:
dist_position_vec_u = p_pos_u - l_pos; % vector between position vectors
d_line_u = (dist_position_vec_u(1,:).*n_face_u(1,:) + dist_position_vec_u(2,:).*n_face_u(2,:) + dist_position_vec_u(3,:).*n_face_u(3,:))...
    ./(l_dir_norm_u(1,:).*n_face_u(1,:) + l_dir_norm_u(2,:).*n_face_u(2,:) + l_dir_norm_u(3,:).*n_face_u(3,:));

% upper intersection point:
P_intsec_u = l_pos + ((d_line_u')*(l_dir_norm_u'))';
P_intsec_u_z = P_intsec_u(3,:);

zloc_part = (z_glob - P_intsec_l_z')./(P_intsec_u_z' - P_intsec_l_z');

%%

xP_start = xloc_part;
yP_start = yloc_part;
zP_start = zloc_part;

if on_gpu
    xP_start = gpuArray(xP_start);
    yP_start = gpuArray(yP_start);
    zP_start = gpuArray(zP_start);
end

%% Now the actual tracking is initialized

current_element_of_part_i = starting_ele_part;

matrix_of_all_exit_points_x = NaN(preallocation_no_iterations,no_particles);
matrix_of_all_exit_points_y = NaN(preallocation_no_iterations,no_particles);
matrix_of_all_exit_points_z = NaN(preallocation_no_iterations,no_particles);
matrix_of_all_flight_times = NaN(preallocation_no_iterations,no_particles);
matrix_of_all_element_IDs = NaN(preallocation_no_iterations,no_particles);
matrix_of_all_zone_IDs = NaN(preallocation_no_iterations,no_particles);

matrix_of_all_exit_points_x(1,:) = xloc_part;
matrix_of_all_exit_points_y(1,:) = yloc_part;
matrix_of_all_exit_points_z(1,:) = zloc_part_glob;

no_poss_neighbs = 5;

cell_PT_x = {1,no_particles};
cell_PT_y = {1,no_particles};
cell_PT_z = {1,no_particles};
cell_PT_t = {1,no_particles};
cell_PT_elemIDs = {1,no_particles};
cell_PT_zoneIDs = {1,no_particles};

if on_gpu
    
    if ~gpu_slim
        current_element_of_part_i = gpuArray(current_element_of_part_i);
        no_poss_neighbs = gpuArray(no_poss_neighbs);
    end
    
    
    matrix_of_all_exit_points_x = gpuArray(matrix_of_all_exit_points_x);
    matrix_of_all_exit_points_y = gpuArray(matrix_of_all_exit_points_y);
    matrix_of_all_exit_points_z = gpuArray(matrix_of_all_exit_points_z);
    matrix_of_all_flight_times = gpuArray(matrix_of_all_flight_times);
    matrix_of_all_element_IDs = gpuArray(matrix_of_all_element_IDs);
    matrix_of_all_zone_IDs = gpuArray(matrix_of_all_zone_IDs);
    
    xloc_part = gpuArray(xloc_part);
    yloc_part = gpuArray(yloc_part);
    zloc_part_glob = gpuArray(zloc_part_glob);
    zloc_part = gpuArray(zloc_part);
    
end


tracking_time_step = ones(1,no_particles);

tracking_finished = false(1,no_particles);

while_iter = 1;
P_entry_x = xP_start';
P_entry_y = yP_start';
P_entry_z = zP_start';


cumulative_time_jj = zeros(no_particles,1);

last_Lf_ending_nodes = zeros(no_particles,4);

exit_horizontal = false(1,no_particles);
exit_on_top = false(1,no_particles);
exit_on_bottom = false(1,no_particles);
tracking_finished = false(1,no_particles);

time_of_flight_so_far = zeros(1,no_particles);

transient_indices = ones(1,no_particles);
transient_indices_all = transient_indices;

time_is_exceeded = false(1,no_particles);
cond_for_max_simulation_time_all = false(1,no_particles);

if on_gpu
    
    cumulative_time_jj = gpuArray(cumulative_time_jj);
    time_of_flight_so_far = gpuArray(time_of_flight_so_far);
    if ~gpu_slim
        tracking_time_step = gpuArray(tracking_time_step);
        while_iter = gpuArray(while_iter);
        
        last_Lf_ending_nodes = gpuArray(last_Lf_ending_nodes);
        exit_horizontal = gpuArray(exit_horizontal);
        exit_on_top = gpuArray(exit_on_top);
        exit_on_bottom = gpuArray(exit_on_bottom);
        tracking_finished = gpuArray(tracking_finished);
        
        transient_indices = gpuArray(transient_indices);
        transient_indices_all = gpuArray(transient_indices_all);
        time_is_exceeded = gpuArray(time_is_exceeded);
        cond_for_max_simulation_time_all = gpuArray(cond_for_max_simulation_time_all);
    end
    
    
end

if transient
    next_time_update = ones(1,no_particles)* next_time_update;
    
end

if on_gpu
    next_time_update = gpuArray(next_time_update);
end

%% Legend for rewriting it in matrix-form for transient simulations

% a_x(element_index) -> a_xmat(sub2ind(size_a_b_mat,element_index,transient_indices))'
% a_y(element_index) -> a_ymat(sub2ind(size_a_b_mat,element_index,transient_indices))'
% a_z(element_index) -> a_zmat(sub2ind(size_a_b_mat,element_index,transient_indices))'
% b_div_horiz(element_index) -> b_div_horiz_mat(sub2ind(size_a_b_mat,element_index,transient_indices))'
% b_div_vert(element_index) -> b_div_vert_mat(sub2ind(size_a_b_mat,element_index,transient_indices))'
% q_locn(element_index,4) -> q_locn_mat(sub2ind(size_qlocn_mat,element_index,[4 4 4],transient_indices))'
% q_locn(element_index,5) -> q_locn_mat(sub2ind(size_qlocn_mat,element_index,[5 5 5],transient_indices))'

rel_indices_v = 1:1:no_particles;

size_a_b_mat = [no_elem no_of_tracking_time_steps];
size_qlocn_mat = [no_elem 5 no_of_tracking_time_steps];

if ~transient
    a_xmat = a_xmat';
    a_ymat = a_ymat';
    a_zmat = a_zmat';
    b_div_horiz_mat = b_div_horiz_mat';
    b_div_vert_mat = b_div_vert_mat';
end

rel_indices = rel_indices_v;
terminatedWell_and_failed_relIDs = [];

if on_gpu && ~gpu_slim
    rel_indices_v = gpuArray(rel_indices_v);
    size_a_b_mat = gpuArray(size_a_b_mat);
    size_qlocn_mat = gpuArray(size_qlocn_mat);
    rel_indices = gpuArray(rel_indices);
    terminatedWell_and_failed_relIDs = gpuArray(terminatedWell_and_failed_relIDs);
    no_particles = gpuArray(no_particles);
end


%% Actual tracking loop:
% Some code is not embedded in a function to increase the speed of the tracking code
% a tiny bit more, also I don't like to have lot's of individual files..

tic

while ((length(tracking_finished)>1 || sum(tracking_finished)~=1)) && ~isempty(tracking_finished)
    
    time_is_exceeded = time_is_exceeded(~tracking_finished);
    transient_indices = transient_indices(~tracking_finished);
    no_active_particles = sum(~tracking_finished);
    element_index = current_element_of_part_i(~tracking_finished);
    
    time_of_flights = ones(1,no_poss_neighbs,no_active_particles)*-1;
    time_of_flight = ones(1,no_active_particles)*-1;
    
    if on_gpu
        time_of_flights = gpuArray(time_of_flights);
        time_of_flight = gpuArray(time_of_flight);
    end
    
    % Divergence flags
    divflow = false(1,no_active_particles);
    divflow_v = false(1,no_active_particles);
    
    cond_divflow_h = abs(b_div_horiz_mat(sub2ind(size_a_b_mat,element_index,transient_indices))')>crit_div_horiz;
    cond_divflow_v = abs(b_div_vert_mat(sub2ind(size_a_b_mat,element_index,transient_indices))')>crit_div_vert;
    
    divflow(cond_divflow_h) = true;
    divflow_v(cond_divflow_v) = true;
    
    divflow = divflow';
    divflow_v = divflow_v';
    
    tracking_finished_for_parts = tracking_finished;
    tracking_finished_for_parts(terminatedWell_and_failed_relIDs) = [];
    
    P_entry_x = P_entry_x(~tracking_finished_for_parts);
    P_entry_y = P_entry_y(~tracking_finished_for_parts);
    P_entry_z = P_entry_z(~tracking_finished_for_parts);
    
    if sum(tracking_finished)>=1
        where_tracking_finished = find(1==tracking_finished);
        
        last_Lf_ending_nodes(where_tracking_finished,:) = [];
        exit_on_top(where_tracking_finished) = [];
        exit_on_bottom(where_tracking_finished) = [];
        
        rel_indices(where_tracking_finished) = [];
        
    end
    
    % REMINDER P_entry_x and P_entry_y only for active particles
    v_x_current = a_xmat(sub2ind(size_a_b_mat,element_index,transient_indices))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index,transient_indices))'.*P_entry_x'; % IS this correct????
    v_y_current = a_ymat(sub2ind(size_a_b_mat,element_index,transient_indices))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index,transient_indices))'.*P_entry_y';
    
    v_xy_norm = sqrt((v_x_current.^2) + (v_y_current.^2));
    v_x_norm = v_x_current./v_xy_norm;
    v_y_norm = v_y_current./v_xy_norm;
    
    v_xy_norm_vec = [v_x_norm';v_y_norm'];
    s = 42.0; % just because
    
    P_entry_xy = [P_entry_x;P_entry_y];
    
    P2_v_xy = P_entry_xy + s*v_xy_norm_vec;
    P2_v_x = P2_v_xy(1,:);
    P2_v_y = P2_v_xy(2,:);
    
    
    node1 = node_ids(element_index,1);
    node2 = node_ids(element_index,2);
    node3 = node_ids(element_index,3);
    node4 = node_ids(element_index,4);
    node5 = node_ids(element_index,5);
    node6 = node_ids(element_index,6);
    
    xA = coordinates(node1,1);
    yA = coordinates(node1,2);
    zA = coordinates(node1,3);
    
    xB = coordinates(node2,1);
    yB = coordinates(node2,2);
    zB = coordinates(node2,3);
    
    xC = coordinates(node3,1);
    yC = coordinates(node3,2);
    zC = coordinates(node3,3);
    
    xD = coordinates(node4,1);
    yD = coordinates(node4,2);
    zD = coordinates(node4,3);
    
    xE = coordinates(node5,1);
    yE = coordinates(node5,2);
    zE = coordinates(node5,3);
    
    xF = coordinates(node6,1);
    yF = coordinates(node6,2);
    zF = coordinates(node6,3);
    
    Xpri = zeros(6,3,no_active_particles);
    
    % Assign Xpri in a vectorized manner:
    Xpri(1,1,:) = xA;
    Xpri(1,2,:) = yA;
    Xpri(1,3,:) = zA;
    
    Xpri(2,1,:) = xB;
    Xpri(2,2,:) = yB;
    Xpri(2,3,:) = zB;
    
    Xpri(3,1,:) = xC;
    Xpri(3,2,:) = yC;
    Xpri(3,3,:) = zC;
    
    Xpri(4,1,:) = xD;
    Xpri(4,2,:) = yD;
    Xpri(4,3,:) = zD;
    
    Xpri(5,1,:) = xE;
    Xpri(5,2,:) = yE;
    Xpri(5,3,:) = zE;
    
    Xpri(6,1,:) = xF;
    Xpri(6,2,:) = yF;
    Xpri(6,3,:) = zF;
    
    
    A_cur = A_horiz(element_index);
    
    xloc_centro = [(xA + xB + xC)./3,(yA + yB + yC)./3]';
    
    index_of_possible_exit_faces = zeros(1,5,no_active_particles);
    
    x1 = P_entry_x';
    y1 = P_entry_y';
    x2 = P2_v_x';
    y2 = P2_v_y';
    
    
    possible_exit_points = zeros(3,5,no_active_particles);
    
    if on_gpu
        possible_exit_points = gpuArray(possible_exit_points);
    end
    
    last_Lf_ending_nodes1_cmp = last_Lf_ending_nodes(:,1);
    last_Lf_ending_nodes2_cmp = last_Lf_ending_nodes(:,2);
    last_Lf_ending_nodes3_cmp = last_Lf_ending_nodes(:,3);
    last_Lf_ending_nodes4_cmp = last_Lf_ending_nodes(:,4);
    
    cond_neighb1 = (~((node1==last_Lf_ending_nodes1_cmp & node2==last_Lf_ending_nodes2_cmp) | (node2==last_Lf_ending_nodes1_cmp & node1==last_Lf_ending_nodes2_cmp)) | last_Lf_ending_nodes4_cmp==0 );
    
    if sum(cond_neighb1)>0
        % first
        cond_neighb1_nodiv = (cond_neighb1 & (~divflow));
        no_entries = sum(cond_neighb1_nodiv);
        
        if no_entries>0
            
            inds_cond_neighb1_nodiv = find(1==cond_neighb1_nodiv);
            
            x1_1div = x1(cond_neighb1_nodiv);
            y1_1div = y1(cond_neighb1_nodiv);
            x2_1div = x2(cond_neighb1_nodiv);
            y2_1div = y2(cond_neighb1_nodiv);
            
            x3 = xA(cond_neighb1_nodiv);
            y3 = yA(cond_neighb1_nodiv);
            x4 = xB(cond_neighb1_nodiv);
            y4 = yB(cond_neighb1_nodiv);
            
            
            possible_exit_points(1,1,inds_cond_neighb1_nodiv) = (((x4-x3).*((x2_1div.*y1_1div) - (x1_1div.*y2_1div))) - ((x2_1div-x1_1div).*((x4.*y3) - (x3.*y4))))./(((y4-y3).*(x2_1div-x1_1div)) - ((y2_1div-y1_1div).*(x4-x3)));
            possible_exit_points(2,1,inds_cond_neighb1_nodiv) = (((y1_1div-y2_1div).*((x4.*y3) - (x3.*y4))) - ((y3-y4).*((x2_1div.*y1_1div) - (x1_1div.*y2_1div))))./(((y4-y3).*(x2_1div-x1_1div)) - ((y2_1div-y1_1div).*(x4-x3)));
            
            exit_point_vec = [reshape(possible_exit_points(1,1,inds_cond_neighb1_nodiv),[1,no_entries]);reshape(possible_exit_points(2,1,inds_cond_neighb1_nodiv),[1,no_entries])];
            entry_point_vec = [x1_1div';y1_1div'];
            
            sign_of_vec_pointer = exit_point_vec - entry_point_vec;
            norm_vec_nodes_vv = coordinates(node1(cond_neighb1_nodiv),1:2)' - coordinates(node2(cond_neighb1_nodiv),1:2)';
            norm_vec_nodes_v = abs(norm_vec_nodes_vv);
            
            norm_vec_nodes = norm_vec_nodes_v./(sqrt((norm_vec_nodes_v(1,:)).^2 + (norm_vec_nodes_v(2,:)).^2));
            
            
            sign_of_velocities = [v_x_current(cond_neighb1_nodiv)';v_y_current(cond_neighb1_nodiv)'];
            
            sign_of_vec_pointer(sign_of_vec_pointer>-1e-12 & sign_of_vec_pointer<1e-12) = 0;
            sign_of_vec_pointer(sign_of_vec_pointer>0)=1;
            sign_of_vec_pointer(sign_of_vec_pointer<0)=-1;
            
            sign_of_velocities(sign_of_velocities>0)=1;
            sign_of_velocities(sign_of_velocities<0)=-1;
            
            distance_vec = [P_entry_x(cond_neighb1_nodiv);P_entry_y(cond_neighb1_nodiv)] -  exit_point_vec;
            euclid_distance = sqrt(sum(distance_vec.^2));
            
            v_xy_abs = (sqrt((a_xmat(sub2ind(size_a_b_mat,element_index(cond_neighb1_nodiv),transient_indices(cond_neighb1_nodiv)))'.^2) + (a_ymat(sub2ind(size_a_b_mat,element_index(cond_neighb1_nodiv),transient_indices(cond_neighb1_nodiv)))'.^2)));
            time_of_flight = euclid_distance'./v_xy_abs;
            
            cond_of_truth_assign_times = ((sign_of_vec_pointer==sign_of_velocities) | sign_of_velocities==0 ) & norm_vec_nodes~=v_xy_norm_vec(:,inds_cond_neighb1_nodiv);
            sum_cond_of_truth_assign_times = sum(cond_of_truth_assign_times);
            inds_cond_of_truth_assign_times = find(2==sum_cond_of_truth_assign_times);
            
            if sum(sum_cond_of_truth_assign_times)>0
                index_of_possible_exit_faces(1,1,inds_cond_of_truth_assign_times)=1;
                time_of_flights(1,1,inds_cond_of_truth_assign_times) = time_of_flight(inds_cond_of_truth_assign_times);
            end
            
            
        end
        % elseif divflow
        cond_neighb1_div = (cond_neighb1 & divflow);
        no_entries = sum(cond_neighb1_div);
        
        if no_entries>0
            
            inds_cond_neighb1_div = find(1==cond_neighb1_div);
            
            x3 = xA(cond_neighb1_div);
            y3 = yA(cond_neighb1_div);
            x4 = xB(cond_neighb1_div);
            y4 = yB(cond_neighb1_div);
            
            rv = [x4-x3,y4-y3]';
            r_norm = sqrt(sum(rv.^2));
            norm_r = rv./r_norm;
            
            s_param = (((a_ymat(sub2ind(size_a_b_mat,element_index(cond_neighb1_div),transient_indices(cond_neighb1_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb1_div),transient_indices(cond_neighb1_div)))'.*y3)./(a_ymat(sub2ind(size_a_b_mat,element_index(cond_neighb1_div),transient_indices(cond_neighb1_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb1_div),transient_indices(cond_neighb1_div)))'.*P_entry_y(inds_cond_neighb1_div)')) - ((a_xmat(sub2ind(size_a_b_mat,element_index(cond_neighb1_div),transient_indices(cond_neighb1_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb1_div),transient_indices(cond_neighb1_div)))'.*x3)./(a_xmat(sub2ind(size_a_b_mat,element_index(cond_neighb1_div),transient_indices(cond_neighb1_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb1_div),transient_indices(cond_neighb1_div)))'.*P_entry_x(inds_cond_neighb1_div)')))./...
                ((b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb1_div),transient_indices(cond_neighb1_div)))'.*norm_r(1,:)'./(a_xmat(sub2ind(size_a_b_mat,element_index(cond_neighb1_div),transient_indices(cond_neighb1_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb1_div),transient_indices(cond_neighb1_div)))'.*P_entry_x(inds_cond_neighb1_div)')) - (b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb1_div),transient_indices(cond_neighb1_div)))'.*norm_r(2,:)'./(a_ymat(sub2ind(size_a_b_mat,element_index(cond_neighb1_div),transient_indices(cond_neighb1_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb1_div),transient_indices(cond_neighb1_div)))'.*P_entry_y(inds_cond_neighb1_div)')));
            % HERE
            s_param_rep = repmat(s_param,[1,2])';
            
            P_exit_with_div = [x3,y3]' + s_param_rep.*norm_r;
            
            possible_exit_points(1,1,inds_cond_neighb1_div) = P_exit_with_div(1,:);
            possible_exit_points(2,1,inds_cond_neighb1_div) = P_exit_with_div(2,:);
            
            if on_gpu
                time_of_flight = (1./b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb1_div),transient_indices(cond_neighb1_div)))') .* log(complex((a_xmat(sub2ind(size_a_b_mat,element_index(cond_neighb1_div),transient_indices(cond_neighb1_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb1_div),transient_indices(cond_neighb1_div)))'.*P_exit_with_div(1,:)')...
                    ./(a_xmat(sub2ind(size_a_b_mat,element_index(cond_neighb1_div),transient_indices(cond_neighb1_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb1_div),transient_indices(cond_neighb1_div)))'.*P_entry_x((inds_cond_neighb1_div))')));
            else
                time_of_flight = (1./b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb1_div),transient_indices(cond_neighb1_div)))') .* log((a_xmat(sub2ind(size_a_b_mat,element_index(cond_neighb1_div),transient_indices(cond_neighb1_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb1_div),transient_indices(cond_neighb1_div)))'.*P_exit_with_div(1,:)')...
                    ./(a_xmat(sub2ind(size_a_b_mat,element_index(cond_neighb1_div),transient_indices(cond_neighb1_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb1_div),transient_indices(cond_neighb1_div)))'.*P_entry_x((inds_cond_neighb1_div))'));
            end
            
            if ~isreal(time_of_flight)
                imag_time_of_flight = imag(time_of_flight);
                logic_imag_time_of_flight = imag_time_of_flight == 0;
            else
                logic_imag_time_of_flight = true(size(time_of_flight));
            end
            
            cond_of_truth_assign_times = time_of_flight>0 & logic_imag_time_of_flight;
            inds_cond_of_truth_assign_times = find(1==cond_of_truth_assign_times);
            
            if sum(cond_of_truth_assign_times)>0
                time_of_flights(1,1,inds_cond_neighb1_div(inds_cond_of_truth_assign_times)) = time_of_flight(inds_cond_of_truth_assign_times);
                index_of_possible_exit_faces(1,1,inds_cond_neighb1_div(inds_cond_of_truth_assign_times)) = 1;
                
            end
            
        end
        
        if on_gpu
            cond_of_truth_forVert = imag(time_of_flights(1,1,:))==0 & time_of_flights(1,1,:)>0 & index_of_possible_exit_faces(1,1,:)>0;
        else
            cond_of_truth_forVert = isreal(time_of_flights(1,1,:)) & time_of_flights(1,1,:)>0 & index_of_possible_exit_faces(1,1,:)>0;
        end
        
        reshape_cond_of_truth_forVert = reshape(cond_of_truth_forVert,[no_active_particles,1]);
        
        if sum(cond_of_truth_forVert)>0
            
            true_time_divflow_v = reshape_cond_of_truth_forVert & divflow_v;
            
            if sum(true_time_divflow_v)>0
                inds_div_flow_vert = find(1==true_time_divflow_v);
                
                a_z_eff = a_zmat(sub2ind(size_a_b_mat,element_index(inds_div_flow_vert),transient_indices(inds_div_flow_vert)))';%a_z(element_index(inds_div_flow_vert));%/detJ_eff;
                b_div_vert_eff = b_div_vert_mat(sub2ind(size_a_b_mat,element_index(inds_div_flow_vert),transient_indices(inds_div_flow_vert)))';%b_div_vert(element_index(inds_div_flow_vert));%/detJ_eff;
                
                time_of_flights_rel_reshape = reshape(time_of_flights(1,1,inds_div_flow_vert),[length(inds_div_flow_vert),1]);
                
                possible_exit_points(3,1,inds_div_flow_vert) = ((a_z_eff + b_div_vert_eff.*P_entry_z(inds_div_flow_vert)').*exp(b_div_vert_eff.*time_of_flights_rel_reshape) - a_z_eff)./b_div_vert_eff;
            end
            
            true_time_nodivflow_v = reshape_cond_of_truth_forVert & ~divflow_v;
            
            if sum(true_time_nodivflow_v)
                
                inds_nodiv_flow_vert = find(1==true_time_nodivflow_v);
                
                a_z_eff = a_zmat(sub2ind(size_a_b_mat,element_index(inds_nodiv_flow_vert),transient_indices(inds_nodiv_flow_vert)))';%a_z(element_index(inds_nodiv_flow_vert));%/detJ_eff;
                time_of_flights_rel_reshape = reshape(time_of_flights(1,1,inds_nodiv_flow_vert),[length(inds_nodiv_flow_vert),1]);
                
                possible_exit_points(3,1,inds_nodiv_flow_vert) = a_z_eff.*time_of_flights_rel_reshape + P_entry_z(inds_nodiv_flow_vert)';
            end
            
        end
        
    end
    
    cond_neighb2 = (~((node1==last_Lf_ending_nodes1_cmp & node3==last_Lf_ending_nodes2_cmp) | (node3==last_Lf_ending_nodes1_cmp & node1==last_Lf_ending_nodes2_cmp)) | last_Lf_ending_nodes4_cmp==0 );
    
    if sum(cond_neighb2)>0
        % second
        cond_neighb2_nodiv = (cond_neighb2 & (~divflow));
        no_entries = sum(cond_neighb2_nodiv);
        
        if no_entries>0
            
            inds_cond_neighb2_nodiv = find(1==cond_neighb2_nodiv);
            
            x1_1div = x1(cond_neighb2_nodiv);
            y1_1div = y1(cond_neighb2_nodiv);
            x2_1div = x2(cond_neighb2_nodiv);
            y2_1div = y2(cond_neighb2_nodiv);
            
            x3 = xA(cond_neighb2_nodiv);
            y3 = yA(cond_neighb2_nodiv);
            x4 = xC(cond_neighb2_nodiv);
            y4 = yC(cond_neighb2_nodiv);
            
            possible_exit_points(1,2,inds_cond_neighb2_nodiv) = (((x4-x3).*((x2_1div.*y1_1div) - (x1_1div.*y2_1div))) - ((x2_1div-x1_1div).*((x4.*y3) - (x3.*y4))))./(((y4-y3).*(x2_1div-x1_1div)) - ((y2_1div-y1_1div).*(x4-x3)));
            possible_exit_points(2,2,inds_cond_neighb2_nodiv) = (((y1_1div-y2_1div).*((x4.*y3) - (x3.*y4))) - ((y3-y4).*((x2_1div.*y1_1div) - (x1_1div.*y2_1div))))./(((y4-y3).*(x2_1div-x1_1div)) - ((y2_1div-y1_1div).*(x4-x3)));
            
            exit_point_vec = [reshape(possible_exit_points(1,2,inds_cond_neighb2_nodiv),[1,no_entries]);reshape(possible_exit_points(2,2,inds_cond_neighb2_nodiv),[1,no_entries])];
            entry_point_vec = [x1_1div';y1_1div'];
            
            sign_of_vec_pointer = exit_point_vec - entry_point_vec;
            norm_vec_nodes_vv = coordinates(node1(cond_neighb2_nodiv),1:2)' - coordinates(node3(cond_neighb2_nodiv),1:2)';
            norm_vec_nodes_v = abs(norm_vec_nodes_vv);
            
            norm_vec_nodes = norm_vec_nodes_v./(sqrt((norm_vec_nodes_v(1,:)).^2 + (norm_vec_nodes_v(2,:)).^2));
            
            sign_of_velocities = [v_x_current(cond_neighb2_nodiv)';v_y_current(cond_neighb2_nodiv)'];
            
            sign_of_vec_pointer(sign_of_vec_pointer>-1e-12 & sign_of_vec_pointer<1e-12) = 0;
            sign_of_vec_pointer(sign_of_vec_pointer>0)=1;
            sign_of_vec_pointer(sign_of_vec_pointer<0)=-1;
            
            sign_of_velocities(sign_of_velocities>0)=1;
            sign_of_velocities(sign_of_velocities<0)=-1;
            
            distance_vec = [P_entry_x(cond_neighb2_nodiv);P_entry_y(cond_neighb2_nodiv)] -  exit_point_vec;
            euclid_distance = sqrt(sum(distance_vec.^2));
            
            v_xy_abs = (sqrt((a_xmat(sub2ind(size_a_b_mat,element_index(cond_neighb2_nodiv),transient_indices(cond_neighb2_nodiv)))'.^2) + (a_ymat(sub2ind(size_a_b_mat,element_index(cond_neighb2_nodiv),transient_indices(cond_neighb2_nodiv)))'.^2)));
            
            time_of_flight = euclid_distance'./v_xy_abs;
            
            cond_of_truth_assign_times = ((sign_of_vec_pointer==sign_of_velocities) | sign_of_velocities==0 ) & norm_vec_nodes~=v_xy_norm_vec(:,inds_cond_neighb2_nodiv);
            sum_cond_of_truth_assign_times = sum(cond_of_truth_assign_times);
            inds_cond_of_truth_assign_times = find(2==sum_cond_of_truth_assign_times);
            
            if sum(sum_cond_of_truth_assign_times)>0
                index_of_possible_exit_faces(1,2,inds_cond_of_truth_assign_times)=2;
                time_of_flights(1,2,inds_cond_of_truth_assign_times) = time_of_flight(inds_cond_of_truth_assign_times);
            end
            
        end
        
        % elseif divflow
        cond_neighb2_div = (cond_neighb2 & divflow);
        no_entries = sum(cond_neighb2_div);
        
        if no_entries>0
            
            inds_cond_neighb2_div = find(1==cond_neighb2_div);
            
            x3 = xA(cond_neighb2_div);
            y3 = yA(cond_neighb2_div);
            x4 = xC(cond_neighb2_div);
            y4 = yC(cond_neighb2_div);
            
            rv = [x4-x3,y4-y3]';
            r_norm = sqrt(sum(rv.^2));
            norm_r = rv./r_norm;
            
            s_param = (((a_ymat(sub2ind(size_a_b_mat,element_index(cond_neighb2_div),transient_indices(cond_neighb2_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb2_div),transient_indices(cond_neighb2_div)))'.*y3)./(a_ymat(sub2ind(size_a_b_mat,element_index(cond_neighb2_div),transient_indices(cond_neighb2_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb2_div),transient_indices(cond_neighb2_div)))'.*P_entry_y(inds_cond_neighb2_div)')) - ((a_xmat(sub2ind(size_a_b_mat,element_index(cond_neighb2_div),transient_indices(cond_neighb2_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb2_div),transient_indices(cond_neighb2_div)))'.*x3)./(a_xmat(sub2ind(size_a_b_mat,element_index(cond_neighb2_div),transient_indices(cond_neighb2_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb2_div),transient_indices(cond_neighb2_div)))'.*P_entry_x(inds_cond_neighb2_div)')))./...
                ((b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb2_div),transient_indices(cond_neighb2_div)))'.*norm_r(1,:)'./(a_xmat(sub2ind(size_a_b_mat,element_index(cond_neighb2_div),transient_indices(cond_neighb2_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb2_div),transient_indices(cond_neighb2_div)))'.*P_entry_x(inds_cond_neighb2_div)')) - (b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb2_div),transient_indices(cond_neighb2_div)))'.*norm_r(2,:)'./(a_ymat(sub2ind(size_a_b_mat,element_index(cond_neighb2_div),transient_indices(cond_neighb2_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb2_div),transient_indices(cond_neighb2_div)))'.*P_entry_y(inds_cond_neighb2_div)')));
            
            s_param_rep = repmat(s_param,[1,2])';
            
            P_exit_with_div = [x3,y3]' + s_param_rep.*norm_r;
            
            possible_exit_points(1,2,inds_cond_neighb2_div) = P_exit_with_div(1,:);
            possible_exit_points(2,2,inds_cond_neighb2_div) = P_exit_with_div(2,:);
            
            
            if on_gpu
                time_of_flight = (1./b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb2_div),transient_indices(cond_neighb2_div)))') .* log(complex((a_xmat(sub2ind(size_a_b_mat,element_index(cond_neighb2_div),transient_indices(cond_neighb2_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb2_div),transient_indices(cond_neighb2_div)))'.*P_exit_with_div(1,:)')...
                    ./(a_xmat(sub2ind(size_a_b_mat,element_index(cond_neighb2_div),transient_indices(cond_neighb2_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb2_div),transient_indices(cond_neighb2_div)))'.*P_entry_x((inds_cond_neighb2_div))')));
            else
                time_of_flight = (1./b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb2_div),transient_indices(cond_neighb2_div)))') .* log((a_xmat(sub2ind(size_a_b_mat,element_index(cond_neighb2_div),transient_indices(cond_neighb2_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb2_div),transient_indices(cond_neighb2_div)))'.*P_exit_with_div(1,:)')...
                    ./(a_xmat(sub2ind(size_a_b_mat,element_index(cond_neighb2_div),transient_indices(cond_neighb2_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb2_div),transient_indices(cond_neighb2_div)))'.*P_entry_x((inds_cond_neighb2_div))'));
            end
            
            
            if ~isreal(time_of_flight)
                imag_time_of_flight = imag(time_of_flight);
                logic_imag_time_of_flight = imag_time_of_flight == 0;
            else
                logic_imag_time_of_flight = true(size(time_of_flight));
            end
            
            cond_of_truth_assign_times = time_of_flight>0 & logic_imag_time_of_flight;
            inds_cond_of_truth_assign_times = find(1==cond_of_truth_assign_times);
            
            if sum(cond_of_truth_assign_times)>0
                time_of_flights(1,2,inds_cond_neighb2_div(inds_cond_of_truth_assign_times)) = time_of_flight(inds_cond_of_truth_assign_times);
                index_of_possible_exit_faces(1,2,inds_cond_neighb2_div(inds_cond_of_truth_assign_times)) = 2;
            end
            
        end
        
        if on_gpu
            cond_of_truth_forVert = imag(time_of_flights(1,2,:))==0 & time_of_flights(1,2,:)>0 & index_of_possible_exit_faces(1,2,:)>0;
        else
            cond_of_truth_forVert = isreal(time_of_flights(1,2,:)) & time_of_flights(1,2,:)>0 & index_of_possible_exit_faces(1,2,:)>0;
        end
        
        reshape_cond_of_truth_forVert = reshape(cond_of_truth_forVert,[no_active_particles,1]);
        
        if sum(cond_of_truth_forVert)>0
            
            true_time_divflow_v = reshape_cond_of_truth_forVert & divflow_v;
            
            if sum(true_time_divflow_v)>0
                inds_div_flow_vert = find(1==true_time_divflow_v);
                
                a_z_eff = a_zmat(sub2ind(size_a_b_mat,element_index(inds_div_flow_vert),transient_indices(inds_div_flow_vert)))';%a_z(element_index(inds_div_flow_vert));%/detJ_eff;
                b_div_vert_eff = b_div_vert_mat(sub2ind(size_a_b_mat,element_index(inds_div_flow_vert),transient_indices(inds_div_flow_vert)))';
                
                time_of_flights_rel_reshape = reshape(time_of_flights(1,2,inds_div_flow_vert),[length(inds_div_flow_vert),1]);
                
                possible_exit_points(3,2,inds_div_flow_vert) = ((a_z_eff + b_div_vert_eff.*P_entry_z(inds_div_flow_vert)').*exp(b_div_vert_eff.*time_of_flights_rel_reshape) - a_z_eff)./b_div_vert_eff;
            end
            
            true_time_nodivflow_v = reshape_cond_of_truth_forVert & ~divflow_v;
            
            if sum(true_time_nodivflow_v)
                
                inds_nodiv_flow_vert = find(1==true_time_nodivflow_v);
                
                a_z_eff = a_zmat(sub2ind(size_a_b_mat,element_index(inds_nodiv_flow_vert),transient_indices(inds_nodiv_flow_vert)))';
                
                time_of_flights_rel_reshape = reshape(time_of_flights(1,2,inds_nodiv_flow_vert),[length(inds_nodiv_flow_vert),1]);
                
                possible_exit_points(3,2,inds_nodiv_flow_vert) = a_z_eff.*time_of_flights_rel_reshape + P_entry_z(inds_nodiv_flow_vert)';
            end
            
        end
        
        
    end
    
    
    % HERE IS THE END OF THE 2ND NEIGHBOUR
    
    cond_neighb3 = (~((node2==last_Lf_ending_nodes1_cmp & node3==last_Lf_ending_nodes2_cmp) | (node3==last_Lf_ending_nodes1_cmp & node2==last_Lf_ending_nodes2_cmp)) | last_Lf_ending_nodes4_cmp==0 );
    
    if sum(cond_neighb3)>0
        % third
        cond_neighb3_nodiv = (cond_neighb3 & (~divflow));
        no_entries = sum(cond_neighb3_nodiv);
        
        if no_entries>0
            
            inds_cond_neighb3_nodiv = find(1==cond_neighb3_nodiv);
            
            x1_1div = x1(cond_neighb3_nodiv);
            y1_1div = y1(cond_neighb3_nodiv);
            x2_1div = x2(cond_neighb3_nodiv);
            y2_1div = y2(cond_neighb3_nodiv);
            
            x3 = xB(cond_neighb3_nodiv);
            y3 = yB(cond_neighb3_nodiv);
            x4 = xC(cond_neighb3_nodiv);
            y4 = yC(cond_neighb3_nodiv);
            
            
            possible_exit_points(1,3,inds_cond_neighb3_nodiv) = (((x4-x3).*((x2_1div.*y1_1div) - (x1_1div.*y2_1div))) - ((x2_1div-x1_1div).*((x4.*y3) - (x3.*y4))))./(((y4-y3).*(x2_1div-x1_1div)) - ((y2_1div-y1_1div).*(x4-x3)));
            possible_exit_points(2,3,inds_cond_neighb3_nodiv) = (((y1_1div-y2_1div).*((x4.*y3) - (x3.*y4))) - ((y3-y4).*((x2_1div.*y1_1div) - (x1_1div.*y2_1div))))./(((y4-y3).*(x2_1div-x1_1div)) - ((y2_1div-y1_1div).*(x4-x3)));
            
            exit_point_vec = [reshape(possible_exit_points(1,3,inds_cond_neighb3_nodiv),[1,no_entries]);reshape(possible_exit_points(2,3,inds_cond_neighb3_nodiv),[1,no_entries])];
            entry_point_vec = [x1_1div';y1_1div'];
            
            sign_of_vec_pointer = exit_point_vec - entry_point_vec;
            norm_vec_nodes_vv = coordinates(node2(cond_neighb3_nodiv),1:2)' - coordinates(node3(cond_neighb3_nodiv),1:2)';
            norm_vec_nodes_v = abs(norm_vec_nodes_vv);
            
            norm_vec_nodes = norm_vec_nodes_v./(sqrt((norm_vec_nodes_v(1,:)).^2 + (norm_vec_nodes_v(2,:)).^2));
            
            sign_of_velocities = [v_x_current(cond_neighb3_nodiv)';v_y_current(cond_neighb3_nodiv)'];
            
            sign_of_vec_pointer(sign_of_vec_pointer>-1e-12 & sign_of_vec_pointer<1e-12) = 0;
            sign_of_vec_pointer(sign_of_vec_pointer>0)=1;
            sign_of_vec_pointer(sign_of_vec_pointer<0)=-1;
            
            sign_of_velocities(sign_of_velocities>0)=1;
            sign_of_velocities(sign_of_velocities<0)=-1;
            
            distance_vec = [P_entry_x(cond_neighb3_nodiv);P_entry_y(cond_neighb3_nodiv)] -  exit_point_vec;
            euclid_distance = sqrt(sum(distance_vec.^2));
            
            v_xy_abs = (sqrt((a_xmat(sub2ind(size_a_b_mat,element_index(cond_neighb3_nodiv),transient_indices(cond_neighb3_nodiv)))'.^2) + (a_ymat(sub2ind(size_a_b_mat,element_index(cond_neighb3_nodiv),transient_indices(cond_neighb3_nodiv)))'.^2)));
            
            time_of_flight = euclid_distance'./v_xy_abs;
            
            cond_of_truth_assign_times = ((sign_of_vec_pointer==sign_of_velocities) | sign_of_velocities==0 ) & norm_vec_nodes~=v_xy_norm_vec(:,inds_cond_neighb3_nodiv);
            sum_cond_of_truth_assign_times = sum(cond_of_truth_assign_times);
            inds_cond_of_truth_assign_times = find(2==sum_cond_of_truth_assign_times);
            
            if sum(sum_cond_of_truth_assign_times)>0
                index_of_possible_exit_faces(1,3,inds_cond_of_truth_assign_times)=3;
                time_of_flights(1,3,inds_cond_of_truth_assign_times) = time_of_flight(inds_cond_of_truth_assign_times);
                
            end
            
            
        end
        % elseif divflow
        cond_neighb3_div = (cond_neighb3 & divflow);
        no_entries = sum(cond_neighb3_div);
        
        if no_entries>0
            
            inds_cond_neighb3_div = find(1==cond_neighb3_div);
            
            x3 = xB(cond_neighb3_div);
            y3 = yB(cond_neighb3_div);
            x4 = xC(cond_neighb3_div);
            y4 = yC(cond_neighb3_div);
            
            rv = [x4-x3,y4-y3]';
            r_norm = sqrt(sum(rv.^2));
            norm_r = rv./r_norm;
            
            s_param = (((a_ymat(sub2ind(size_a_b_mat,element_index(cond_neighb3_div),transient_indices(cond_neighb3_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb3_div),transient_indices(cond_neighb3_div)))'.*y3)./(a_ymat(sub2ind(size_a_b_mat,element_index(cond_neighb3_div),transient_indices(cond_neighb3_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb3_div),transient_indices(cond_neighb3_div)))'.*P_entry_y(inds_cond_neighb3_div)')) - ((a_xmat(sub2ind(size_a_b_mat,element_index(cond_neighb3_div),transient_indices(cond_neighb3_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb3_div),transient_indices(cond_neighb3_div)))'.*x3)./(a_xmat(sub2ind(size_a_b_mat,element_index(cond_neighb3_div),transient_indices(cond_neighb3_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb3_div),transient_indices(cond_neighb3_div)))'.*P_entry_x(inds_cond_neighb3_div)')))./...
                ((b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb3_div),transient_indices(cond_neighb3_div)))'.*norm_r(1,:)'./(a_xmat(sub2ind(size_a_b_mat,element_index(cond_neighb3_div),transient_indices(cond_neighb3_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb3_div),transient_indices(cond_neighb3_div)))'.*P_entry_x(inds_cond_neighb3_div)')) - (b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb3_div),transient_indices(cond_neighb3_div)))'.*norm_r(2,:)'./(a_ymat(sub2ind(size_a_b_mat,element_index(cond_neighb3_div),transient_indices(cond_neighb3_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb3_div),transient_indices(cond_neighb3_div)))'.*P_entry_y(inds_cond_neighb3_div)')));
            
            s_param_rep = repmat(s_param,[1,2])';
            
            P_exit_with_div = [x3,y3]' + s_param_rep.*norm_r;
            
            possible_exit_points(1,3,inds_cond_neighb3_div) = P_exit_with_div(1,:);
            possible_exit_points(2,3,inds_cond_neighb3_div) = P_exit_with_div(2,:);
            
            if on_gpu
                time_of_flight = (1./b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb3_div),transient_indices(cond_neighb3_div)))') .* log(complex((a_xmat(sub2ind(size_a_b_mat,element_index(cond_neighb3_div),transient_indices(cond_neighb3_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb3_div),transient_indices(cond_neighb3_div)))'.*P_exit_with_div(1,:)')...
                    ./(a_xmat(sub2ind(size_a_b_mat,element_index(cond_neighb3_div),transient_indices(cond_neighb3_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb3_div),transient_indices(cond_neighb3_div)))'.*P_entry_x((inds_cond_neighb3_div))')));
            else
                time_of_flight = (1./b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb3_div),transient_indices(cond_neighb3_div)))') .* log((a_xmat(sub2ind(size_a_b_mat,element_index(cond_neighb3_div),transient_indices(cond_neighb3_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb3_div),transient_indices(cond_neighb3_div)))'.*P_exit_with_div(1,:)')...
                    ./(a_xmat(sub2ind(size_a_b_mat,element_index(cond_neighb3_div),transient_indices(cond_neighb3_div)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(cond_neighb3_div),transient_indices(cond_neighb3_div)))'.*P_entry_x((inds_cond_neighb3_div))'));
            end
            
            
            if ~isreal(time_of_flight)
                imag_time_of_flight = imag(time_of_flight);
                logic_imag_time_of_flight = imag_time_of_flight == 0;
            else
                logic_imag_time_of_flight = true(size(time_of_flight));
            end
            
            cond_of_truth_assign_times = time_of_flight>0 & logic_imag_time_of_flight;
            inds_cond_of_truth_assign_times = find(1==cond_of_truth_assign_times);
            
            if sum(cond_of_truth_assign_times)>0
                
                time_of_flights(1,3,inds_cond_neighb3_div(inds_cond_of_truth_assign_times)) = time_of_flight(inds_cond_of_truth_assign_times);
                index_of_possible_exit_faces(1,3,inds_cond_neighb3_div(inds_cond_of_truth_assign_times)) = 3;
                
            end
            
        end
        
        if on_gpu
            cond_of_truth_forVert =imag(time_of_flights(1,3,:))==0 & time_of_flights(1,3,:)>0 & index_of_possible_exit_faces(1,3,:)>0;
        else
            cond_of_truth_forVert = isreal(time_of_flights(1,3,:)) & time_of_flights(1,3,:)>0 & index_of_possible_exit_faces(1,3,:)>0;
        end
        
        reshape_cond_of_truth_forVert = reshape(cond_of_truth_forVert,[no_active_particles,1]);
        
        if sum(cond_of_truth_forVert)>0
            
            true_time_divflow_v = reshape_cond_of_truth_forVert & divflow_v;
            
            if sum(true_time_divflow_v)>0
                inds_div_flow_vert = find(1==true_time_divflow_v);
                % divflow_v
                
                a_z_eff = a_zmat(sub2ind(size_a_b_mat,element_index(inds_div_flow_vert),transient_indices(inds_div_flow_vert)))';
                b_div_vert_eff = b_div_vert_mat(sub2ind(size_a_b_mat,element_index(inds_div_flow_vert),transient_indices(inds_div_flow_vert)))';
                
                time_of_flights_rel_reshape = reshape(time_of_flights(1,3,inds_div_flow_vert),[length(inds_div_flow_vert),1]);
                
                possible_exit_points(3,3,inds_div_flow_vert) = ((a_z_eff + b_div_vert_eff.*P_entry_z(inds_div_flow_vert)').*exp(b_div_vert_eff.*time_of_flights_rel_reshape) - a_z_eff)./b_div_vert_eff;
            end
            
            true_time_nodivflow_v = reshape_cond_of_truth_forVert & ~divflow_v;
            
            if sum(true_time_nodivflow_v)
                
                inds_nodiv_flow_vert = find(1==true_time_nodivflow_v);
                
                a_z_eff = a_zmat(sub2ind(size_a_b_mat,element_index(inds_nodiv_flow_vert),transient_indices(inds_nodiv_flow_vert)))';
                
                time_of_flights_rel_reshape = reshape(time_of_flights(1,3,inds_nodiv_flow_vert),[length(inds_nodiv_flow_vert),1]);
                
                possible_exit_points(3,3,inds_nodiv_flow_vert) = a_z_eff.*time_of_flights_rel_reshape + P_entry_z(inds_nodiv_flow_vert)';
            end
            
        end
        
        
    end
    
    % HERE IS THE END OF THE 3RD NEIGHBOUR
    
    v_z_top = a_zmat(sub2ind(size_a_b_mat,element_index,transient_indices))' + b_div_vert_mat(sub2ind(size_a_b_mat,element_index,transient_indices))'*1;
    v_z_bottom =  a_zmat(sub2ind(size_a_b_mat,element_index,transient_indices))';
    
    cond_neighb4 = ~exit_on_top' & (sign(q_locn_mat(sub2ind(size_qlocn_mat,element_index,ones(1,no_active_particles)*4,transient_indices))')==1 | sign(q_locn_mat(sub2ind(size_qlocn_mat,element_index,ones(1,no_active_particles)*4,transient_indices))')==0 );
    
    if sum(cond_neighb4)>0
        
        inds_poss_exit_bot = find(1==cond_neighb4);
        possible_exit_points(3,4,inds_poss_exit_bot) = 1;
        
        if sum(~divflow_v)>0
            inds_nodiv_vert = find(0==divflow_v);
            time_of_flight = -1*((P_entry_z(inds_nodiv_vert))'./(a_zmat(sub2ind(size_a_b_mat,element_index(inds_nodiv_vert),transient_indices(inds_nodiv_vert)))'));
            
            if ~isreal(time_of_flight)
                imag_time_of_flight = imag(time_of_flight);
                logic_imag_time_of_flight = imag_time_of_flight == 0;
            else
                logic_imag_time_of_flight = true(size(time_of_flight));
            end
            
            cond_intmed_nodiv_bot = time_of_flight>0 & logic_imag_time_of_flight;
            
            if sum(cond_intmed_nodiv_bot)>0
                true_inds_pos_bot_nodiv = find(1==cond_intmed_nodiv_bot);
                index_of_possible_exit_faces(1,4,inds_nodiv_vert(true_inds_pos_bot_nodiv)) = 4;
                time_of_flights(1,4,inds_nodiv_vert(true_inds_pos_bot_nodiv)) = time_of_flight(true_inds_pos_bot_nodiv);
                
            end
            
        end
        
        if sum(divflow_v)>0
            inds_div_vert = find(1==divflow_v);
            
            if on_gpu
                time_of_flight = (1./b_div_vert_mat(sub2ind(size_a_b_mat,element_index(inds_div_vert),transient_indices(inds_div_vert)))').* log(complex((a_zmat(sub2ind(size_a_b_mat,element_index(inds_div_vert),transient_indices(inds_div_vert)))' + b_div_vert_mat(sub2ind(size_a_b_mat,element_index(inds_div_vert),transient_indices(inds_div_vert)))'*0)...
                    ./(a_zmat(sub2ind(size_a_b_mat,element_index(inds_div_vert),transient_indices(inds_div_vert)))' + b_div_vert_mat(sub2ind(size_a_b_mat,element_index(inds_div_vert),transient_indices(inds_div_vert)))'.*P_entry_z(inds_div_vert)')));
            else
                time_of_flight = (1./b_div_vert_mat(sub2ind(size_a_b_mat,element_index(inds_div_vert),transient_indices(inds_div_vert)))').* log((a_zmat(sub2ind(size_a_b_mat,element_index(inds_div_vert),transient_indices(inds_div_vert)))' + b_div_vert_mat(sub2ind(size_a_b_mat,element_index(inds_div_vert),transient_indices(inds_div_vert)))'*0)...
                    ./(a_zmat(sub2ind(size_a_b_mat,element_index(inds_div_vert),transient_indices(inds_div_vert)))' + b_div_vert_mat(sub2ind(size_a_b_mat,element_index(inds_div_vert),transient_indices(inds_div_vert)))'.*P_entry_z(inds_div_vert)'));
            end
            
            if ~isreal(time_of_flight)
                imag_time_of_flight = imag(time_of_flight);
                logic_imag_time_of_flight = imag_time_of_flight == 0;
            else
                logic_imag_time_of_flight = true(size(time_of_flight));
            end
            
            cond_intmed_div_bot = time_of_flight>0 &  logic_imag_time_of_flight;
            
            if sum(cond_intmed_div_bot)>0
                true_inds_pos_bot_div = find(1==cond_intmed_div_bot);
                index_of_possible_exit_faces(1,4,inds_div_vert(true_inds_pos_bot_div)) = 4;
                time_of_flights(1,4,inds_div_vert(true_inds_pos_bot_div)) = time_of_flight(true_inds_pos_bot_div);
                
            end
            
        end
        
        reshape_time_of_flights_bot = reshape(time_of_flights(1,4,:),[no_active_particles,1]);
        
        cond_nodiv_horiz_exbot = ~divflow & reshape_time_of_flights_bot>0;
        
        if sum(cond_nodiv_horiz_exbot)>0
            inds_cond_nodiv_horiz_exbot = find(1==cond_nodiv_horiz_exbot);
            time_of_flights_rel_reshape = reshape(time_of_flights(1,4,inds_cond_nodiv_horiz_exbot),[length(inds_cond_nodiv_horiz_exbot),1]);
            
            possible_exit_points(1,4,inds_cond_nodiv_horiz_exbot) = a_xmat(sub2ind(size_a_b_mat,element_index(inds_cond_nodiv_horiz_exbot),transient_indices(inds_cond_nodiv_horiz_exbot)))'.*time_of_flights_rel_reshape + P_entry_x(inds_cond_nodiv_horiz_exbot)';
            possible_exit_points(2,4,inds_cond_nodiv_horiz_exbot) = a_ymat(sub2ind(size_a_b_mat,element_index(inds_cond_nodiv_horiz_exbot),transient_indices(inds_cond_nodiv_horiz_exbot)))'.*time_of_flights_rel_reshape + P_entry_y(inds_cond_nodiv_horiz_exbot)';
            
        end
        
        cond_div_horiz_exbot = divflow & reshape_time_of_flights_bot>0;
        
        if sum(cond_div_horiz_exbot)>0
            inds_cond_div_horiz_exbot = find(1==cond_div_horiz_exbot);
            time_of_flights_rel_reshape = reshape(time_of_flights(1,4,inds_cond_div_horiz_exbot),[length(inds_cond_div_horiz_exbot),1]);
            
            possible_exit_points(1,4,inds_cond_div_horiz_exbot)= ((a_xmat(sub2ind(size_a_b_mat,element_index(inds_cond_div_horiz_exbot),transient_indices(inds_cond_div_horiz_exbot)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(inds_cond_div_horiz_exbot),transient_indices(inds_cond_div_horiz_exbot)))'.*P_entry_x(inds_cond_div_horiz_exbot)').*exp(b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(inds_cond_div_horiz_exbot),transient_indices(inds_cond_div_horiz_exbot)))'.* time_of_flights_rel_reshape) - a_xmat(sub2ind(size_a_b_mat,element_index(inds_cond_div_horiz_exbot),transient_indices(inds_cond_div_horiz_exbot)))')./b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(inds_cond_div_horiz_exbot),transient_indices(inds_cond_div_horiz_exbot)))';
            possible_exit_points(2,4,inds_cond_div_horiz_exbot)= ((a_ymat(sub2ind(size_a_b_mat,element_index(inds_cond_div_horiz_exbot),transient_indices(inds_cond_div_horiz_exbot)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(inds_cond_div_horiz_exbot),transient_indices(inds_cond_div_horiz_exbot)))'.*P_entry_y(inds_cond_div_horiz_exbot)').*exp(b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(inds_cond_div_horiz_exbot),transient_indices(inds_cond_div_horiz_exbot)))'.* time_of_flights_rel_reshape) - a_ymat(sub2ind(size_a_b_mat,element_index(inds_cond_div_horiz_exbot),transient_indices(inds_cond_div_horiz_exbot)))')./b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(inds_cond_div_horiz_exbot),transient_indices(inds_cond_div_horiz_exbot)))';
            
        end
        
    end
    
    
    cond_neighb5 = ~exit_on_bottom' & (sign(q_locn_mat(sub2ind(size_qlocn_mat,element_index,ones(1,no_active_particles)*5,transient_indices))')==1 | sign(q_locn_mat(sub2ind(size_qlocn_mat,element_index,ones(1,no_active_particles)*5,transient_indices))')==0 );
    
    if sum(cond_neighb5)>0
        
        inds_poss_exit_bot = find(1==cond_neighb5);
        possible_exit_points(3,5,inds_poss_exit_bot) = 0;
        
        if sum(~divflow_v)>0
            
            inds_nodiv_vert = find(0==divflow_v);
            
            time_of_flight = ((1-P_entry_z(inds_nodiv_vert)')./(a_zmat(sub2ind(size_a_b_mat,element_index(inds_nodiv_vert),transient_indices(inds_nodiv_vert)))'));
            
            if ~isreal(time_of_flight)
                imag_time_of_flight = imag(time_of_flight);
                logic_imag_time_of_flight = imag_time_of_flight == 0;
            else
                logic_imag_time_of_flight = true(size(time_of_flight));
            end
            
            cond_intmed_nodiv_bot = time_of_flight>0 & logic_imag_time_of_flight;
            
            if sum(cond_intmed_nodiv_bot)>0
                true_inds_pos_bot_nodiv = find(1==cond_intmed_nodiv_bot);
                index_of_possible_exit_faces(1,5,inds_nodiv_vert(true_inds_pos_bot_nodiv)) = 5;
                time_of_flights(1,5,inds_nodiv_vert(true_inds_pos_bot_nodiv)) = time_of_flight(true_inds_pos_bot_nodiv);
            end
            
        end
        
        
        if sum(divflow_v)>0
            
            inds_div_vert = find(1==divflow_v);
            
            if on_gpu
                time_of_flight = (1./b_div_vert_mat(sub2ind(size_a_b_mat,element_index(inds_div_vert),transient_indices(inds_div_vert)))').* log(complex((a_zmat(sub2ind(size_a_b_mat,element_index(inds_div_vert),transient_indices(inds_div_vert)))' + b_div_vert_mat(sub2ind(size_a_b_mat,element_index(inds_div_vert),transient_indices(inds_div_vert)))'*1)...
                    ./(a_zmat(sub2ind(size_a_b_mat,element_index(inds_div_vert),transient_indices(inds_div_vert)))' + b_div_vert_mat(sub2ind(size_a_b_mat,element_index(inds_div_vert),transient_indices(inds_div_vert)))'.*P_entry_z(inds_div_vert)')));
            else
                time_of_flight = (1./b_div_vert_mat(sub2ind(size_a_b_mat,element_index(inds_div_vert),transient_indices(inds_div_vert)))').* log((a_zmat(sub2ind(size_a_b_mat,element_index(inds_div_vert),transient_indices(inds_div_vert)))' + b_div_vert_mat(sub2ind(size_a_b_mat,element_index(inds_div_vert),transient_indices(inds_div_vert)))'*1)...
                    ./(a_zmat(sub2ind(size_a_b_mat,element_index(inds_div_vert),transient_indices(inds_div_vert)))' + b_div_vert_mat(sub2ind(size_a_b_mat,element_index(inds_div_vert),transient_indices(inds_div_vert)))'.*P_entry_z(inds_div_vert)'));
            end
            
            if ~isreal(time_of_flight)
                imag_time_of_flight = imag(time_of_flight);
                logic_imag_time_of_flight = imag_time_of_flight == 0;
            else
                logic_imag_time_of_flight = true(size(time_of_flight));
            end
            
            cond_intmed_div_bot = time_of_flight>0 & logic_imag_time_of_flight;
            
            if sum(cond_intmed_div_bot)>0
                
                true_inds_pos_bot_div = find(1==cond_intmed_div_bot);
                index_of_possible_exit_faces(1,5,inds_div_vert(true_inds_pos_bot_div)) = 5;
                time_of_flights(1,5,inds_div_vert(true_inds_pos_bot_div)) = time_of_flight(true_inds_pos_bot_div);
            end
            
        end
        
        reshape_time_of_flights_top = reshape(time_of_flights(1,5,:),[no_active_particles,1]);
        cond_nodiv_horiz_extop = ~divflow & reshape_time_of_flights_top>0;
        
        if sum(cond_nodiv_horiz_extop)>0
            inds_cond_nodiv_horiz_extop = find(1==cond_nodiv_horiz_extop);
            time_of_flights_rel_reshape = reshape(time_of_flights(1,5,inds_cond_nodiv_horiz_extop),[length(inds_cond_nodiv_horiz_extop),1]);
            
            possible_exit_points(1,5,inds_cond_nodiv_horiz_extop) = a_xmat(sub2ind(size_a_b_mat,element_index(inds_cond_nodiv_horiz_extop),transient_indices(inds_cond_nodiv_horiz_extop)))'.*time_of_flights_rel_reshape + P_entry_x(inds_cond_nodiv_horiz_extop)';
            possible_exit_points(2,5,inds_cond_nodiv_horiz_extop) = a_ymat(sub2ind(size_a_b_mat,element_index(inds_cond_nodiv_horiz_extop),transient_indices(inds_cond_nodiv_horiz_extop)))'.*time_of_flights_rel_reshape + P_entry_y(inds_cond_nodiv_horiz_extop)';
            
        end
        
        cond_div_horiz_extop = divflow & reshape_time_of_flights_top>0;
        
        if sum(cond_div_horiz_extop)>0
            inds_cond_div_horiz_extop = find(1==cond_div_horiz_extop);
            time_of_flights_rel_reshape = reshape(time_of_flights(1,5,inds_cond_div_horiz_extop),[length(inds_cond_div_horiz_extop),1]);
            
            possible_exit_points(1,5,inds_cond_div_horiz_extop)= ((a_xmat(sub2ind(size_a_b_mat,element_index(inds_cond_div_horiz_extop),transient_indices(inds_cond_div_horiz_extop)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(inds_cond_div_horiz_extop),transient_indices(inds_cond_div_horiz_extop)))'.*P_entry_x(inds_cond_div_horiz_extop)').*exp(b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(inds_cond_div_horiz_extop),transient_indices(inds_cond_div_horiz_extop)))'.* time_of_flights_rel_reshape) - a_xmat(sub2ind(size_a_b_mat,element_index(inds_cond_div_horiz_extop),transient_indices(inds_cond_div_horiz_extop)))')./b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(inds_cond_div_horiz_extop),transient_indices(inds_cond_div_horiz_extop)))';
            possible_exit_points(2,5,inds_cond_div_horiz_extop)= ((a_ymat(sub2ind(size_a_b_mat,element_index(inds_cond_div_horiz_extop),transient_indices(inds_cond_div_horiz_extop)))' + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(inds_cond_div_horiz_extop),transient_indices(inds_cond_div_horiz_extop)))'.*P_entry_y(inds_cond_div_horiz_extop)').*exp(b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(inds_cond_div_horiz_extop),transient_indices(inds_cond_div_horiz_extop)))'.* time_of_flights_rel_reshape) - a_ymat(sub2ind(size_a_b_mat,element_index(inds_cond_div_horiz_extop),transient_indices(inds_cond_div_horiz_extop)))')./b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(inds_cond_div_horiz_extop),transient_indices(inds_cond_div_horiz_extop)))';
            
        end
        
        
    end
    
    
    if on_gpu
        time_of_flights = real(time_of_flights);
    end
    
    
    rel_particle_IDs_in_well = [];
    particleElemIDs_in_well = [];
    rel_particle_IDs_of_failure = [];
    
    possibly_failed_indices = ones(no_active_particles,1);
    
    P_entry_vor_x = P_entry_x';
    P_entry_vor_y = P_entry_y';
    P_entry_vor_z = P_entry_z';
    
    %% Time to fly...
    
    sum_index_of_possible_exit_faces = sum(index_of_possible_exit_faces,2);
    sum_index_of_possible_exit_faces_reshape = reshape(sum_index_of_possible_exit_faces,[no_active_particles,1]);
    
    % Wells
    cond_exist_well_truth = sum_index_of_possible_exit_faces_reshape==0 & wells_exist;
    sum_cond_exist_well_truth = sum(cond_exist_well_truth);
    
    if sum_cond_exist_well_truth>0
        ind_cond_exist_well_truth = find(1==cond_exist_well_truth);
        
        particleElemIDs_in_well = zeros(sum_cond_exist_well_truth,1);
        rel_particle_IDs_in_well = zeros(sum_cond_exist_well_truth,1);
        
        find_part_in_well = element_index(ind_cond_exist_well_truth)==ele_list_wells_all;
        sum_find_part_in_well = sum(find_part_in_well,1);
        sum_sum_find_part_in_well = sum(sum_find_part_in_well);
        
        if sum_sum_find_part_in_well>0
            find_sum_find_part_in_well = find(1==sum_find_part_in_well)';
            particleElemIDs_in_well = element_index(ind_cond_exist_well_truth(find_sum_find_part_in_well))';
            rel_particle_IDs_in_well = ind_cond_exist_well_truth(find_sum_find_part_in_well);
            possibly_failed_indices(ind_cond_exist_well_truth(find_sum_find_part_in_well)) = false;
        end
        
    end
    
    if sum_cond_exist_well_truth>0
        if ~isempty(rel_particle_IDs_in_well)
            time_of_flights(:,:,rel_particle_IDs_in_well) = 0;
        end
    end
    
    % possible failure of trajectories (should not happen, but better account for it..)
    condition_of_failure = sum_index_of_possible_exit_faces_reshape==0 & possibly_failed_indices;
    inds_of_failure = [];
    if sum(condition_of_failure)>0
        rel_particle_IDs_of_failure = find(1==condition_of_failure);
        inds_of_failure = rel_indices(rel_particle_IDs_of_failure);
        
        
        matrix_of_all_exit_points_x(while_iter+1,inds_of_failure) = -Inf;
        matrix_of_all_exit_points_y(while_iter+1,inds_of_failure) = -Inf;
        matrix_of_all_exit_points_z(while_iter+1,inds_of_failure) = -Inf;
        matrix_of_all_flight_times(while_iter+1,inds_of_failure) = -Inf;
        fprintf('Something went wrong here => Ask yourself: does your model fulfill the requirements/constraints of this code? \n If you are sure that the answer is "yes", ask Philipp to fix the issue!')
    end
    
    % correct for wells and failed particles
    
    terminatedWell_and_failed_relIDs = [rel_particle_IDs_in_well;rel_particle_IDs_of_failure];
    
    successfull_relIDs = [1:1:no_active_particles]';
    successfull_relIDs(terminatedWell_and_failed_relIDs) = [];
    
    if on_gpu
        terminatedWell_and_failed_relIDs = gpuArray(terminatedWell_and_failed_relIDs);
        successfull_relIDs = gpuArray(successfull_relIDs);
    end
    
    % failed: -Inf
    % well: 0
    
    time_of_flights_inf = time_of_flights(:,:,successfull_relIDs);
    time_of_flights_inf(time_of_flights_inf<0) = Inf;
    [time_of_flights_min,right_exit_indices] = min(time_of_flights_inf,[],2);
    
    current_time_of_flight_vv = reshape(time_of_flights_min,[length(successfull_relIDs),1]);
    right_index_vv = reshape(right_exit_indices,[length(successfull_relIDs),1]);
    
    if on_gpu
        current_time_of_flight = NaN(no_active_particles,1,'gpuArray');
        right_index = NaN(no_active_particles,1,'gpuArray');
    else
        current_time_of_flight = NaN(no_active_particles,1);
        right_index = NaN(no_active_particles,1);
    end
    
    current_time_of_flight(successfull_relIDs) = current_time_of_flight_vv;
    right_index(successfull_relIDs) = right_index_vv;
    
    current_time_of_flight(rel_particle_IDs_in_well) = 0;
    right_index(rel_particle_IDs_in_well) = 0;
    
    current_time_of_flight(rel_particle_IDs_of_failure) = -Inf;
    right_index(rel_particle_IDs_of_failure) = -1;
    
    old_cumulative_time_jj = cumulative_time_jj;
    
    cumulative_time_jj(rel_indices(successfull_relIDs)) = old_cumulative_time_jj(rel_indices(successfull_relIDs)) + (current_time_of_flight(successfull_relIDs));
    
    exit_on_lower_layer = false(no_active_particles,1);
    exit_on_upper_layer = false(no_active_particles,1);
    
    right_index_bot = find(4==right_index);
    right_index_top = find(5==right_index);
    
    if ~isempty(right_index_bot)
        exit_on_lower_layer(right_index_bot) = true;
    end
    
    if ~isempty(right_index_top)
        exit_on_upper_layer(right_index_top) = true;
    end
    
    
    possible_exit_points_longmat = reshape(possible_exit_points,[3,no_active_particles*5]);
    right_matIDs_for_exit_point = ((successfull_relIDs-1)*5 + right_index(successfull_relIDs))';
    exit_point_v =  possible_exit_points_longmat(:,right_matIDs_for_exit_point);
    
    length_failed_inds = length(terminatedWell_and_failed_relIDs);
    
    if on_gpu
        no_nan_to_fill_up = nan(3,length_failed_inds,'gpuArray');
        exit_point = nan(3,no_active_particles,'gpuArray');
    else
        no_nan_to_fill_up = nan(3,length_failed_inds);
        exit_point = nan(3,no_active_particles);
    end
    
    exit_point(:,successfull_relIDs) = exit_point_v;
    
    ending_in_element = false(1,no_active_particles);
    
    logical_successfull_IDs = ones(1,no_active_particles);
    logical_successfull_IDs(~successfull_relIDs) = false;
    
    if transient
        
        cond_for_transient_velocity_update = cumulative_time_jj(rel_indices(successfull_relIDs))> next_time_update(rel_indices(successfull_relIDs))' & next_time_update(rel_indices(successfull_relIDs))'<steady_state_time;
        
        if sum(cond_for_transient_velocity_update)>0
            
            part_inds_toupdate = find(1==cond_for_transient_velocity_update)';
            
            next_time_update_current =  next_time_update(rel_indices(part_inds_toupdate));
            
            tracking_time_step(rel_indices(part_inds_toupdate)) = tracking_time_step(rel_indices(part_inds_toupdate)) + 1;
            next_time_update(rel_indices(part_inds_toupdate)) = target_times(FVM_corr_target_time_IDs(tracking_time_step(rel_indices(part_inds_toupdate))));
            
            current_time_of_flight(part_inds_toupdate) = next_time_update_current - time_of_flight_so_far(rel_indices(part_inds_toupdate));
            
            cond_divflow_horiz = divflow(cond_for_transient_velocity_update)';
            
            
            if sum(cond_divflow_horiz)>0
                exit_point(1,part_inds_toupdate(cond_divflow_horiz)) = ((a_xmat(sub2ind(size_a_b_mat,element_index(part_inds_toupdate(cond_divflow_horiz)),transient_indices(part_inds_toupdate(cond_divflow_horiz)))) + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(part_inds_toupdate(cond_divflow_horiz)),transient_indices(part_inds_toupdate(cond_divflow_horiz)))).*P_entry_x(part_inds_toupdate(cond_divflow_horiz))).*exp(b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(part_inds_toupdate(cond_divflow_horiz)),transient_indices(part_inds_toupdate(cond_divflow_horiz)))).*current_time_of_flight(part_inds_toupdate(cond_divflow_horiz))') - a_xmat(sub2ind(size_a_b_mat,element_index(part_inds_toupdate(cond_divflow_horiz)),transient_indices(part_inds_toupdate(cond_divflow_horiz)))))./b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(part_inds_toupdate(cond_divflow_horiz)),transient_indices(part_inds_toupdate(cond_divflow_horiz))));
                
                exit_point(2,part_inds_toupdate(cond_divflow_horiz)) = ((a_ymat(sub2ind(size_a_b_mat,element_index(part_inds_toupdate(cond_divflow_horiz)),transient_indices(part_inds_toupdate(cond_divflow_horiz)))) + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(part_inds_toupdate(cond_divflow_horiz)),transient_indices(part_inds_toupdate(cond_divflow_horiz)))).*P_entry_y(part_inds_toupdate(cond_divflow_horiz))).*exp(b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(part_inds_toupdate(cond_divflow_horiz)),transient_indices(part_inds_toupdate(cond_divflow_horiz)))).*current_time_of_flight(part_inds_toupdate(cond_divflow_horiz))') - a_ymat(sub2ind(size_a_b_mat,element_index(part_inds_toupdate(cond_divflow_horiz)),transient_indices(part_inds_toupdate(cond_divflow_horiz)))))./b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(part_inds_toupdate(cond_divflow_horiz)),transient_indices(part_inds_toupdate(cond_divflow_horiz))));
            end
            
            if sum(~cond_divflow_horiz)>0
                exit_point(1,part_inds_toupdate(~cond_divflow_horiz)) = a_xmat(sub2ind(size_a_b_mat,element_index(part_inds_toupdate(~cond_divflow_horiz)),transient_indices(part_inds_toupdate(~cond_divflow_horiz)))).*current_time_of_flight(part_inds_toupdate(~cond_divflow_horiz))' + P_entry_x(part_inds_toupdate(~cond_divflow_horiz));
                exit_point(2,part_inds_toupdate(~cond_divflow_horiz)) = a_ymat(sub2ind(size_a_b_mat,element_index(part_inds_toupdate(~cond_divflow_horiz)),transient_indices(part_inds_toupdate(~cond_divflow_horiz)))).*current_time_of_flight(part_inds_toupdate(~cond_divflow_horiz))' + P_entry_y(part_inds_toupdate(~cond_divflow_horiz));
            end
            
            cond_divflow_vert = divflow_v(cond_for_transient_velocity_update)';
            
            if sum(cond_divflow_vert)>0
                exit_point(3,part_inds_toupdate(cond_divflow_vert)) = ((a_zmat(sub2ind(size_a_b_mat,element_index(part_inds_toupdate(cond_divflow_vert)),transient_indices(part_inds_toupdate(cond_divflow_vert)))) + b_div_vert_mat(sub2ind(size_a_b_mat,element_index(part_inds_toupdate(cond_divflow_vert)),transient_indices(part_inds_toupdate(cond_divflow_vert)))).*P_entry_z(part_inds_toupdate(cond_divflow_vert))).*exp(b_div_vert_mat(sub2ind(size_a_b_mat,element_index(part_inds_toupdate(cond_divflow_vert)),transient_indices(part_inds_toupdate(cond_divflow_vert)))).*current_time_of_flight(part_inds_toupdate(cond_divflow_vert))') - a_zmat(sub2ind(size_a_b_mat,element_index(part_inds_toupdate(cond_divflow_vert)),transient_indices(part_inds_toupdate(cond_divflow_vert)))))./b_div_vert_mat(sub2ind(size_a_b_mat,element_index(part_inds_toupdate(cond_divflow_vert)),transient_indices(part_inds_toupdate(cond_divflow_vert))));
            end
            
            if sum(~cond_divflow_vert)>0
                exit_point(3,part_inds_toupdate(~cond_divflow_vert)) = a_zmat(sub2ind(size_a_b_mat,element_index(part_inds_toupdate(~cond_divflow_vert)),transient_indices(part_inds_toupdate(~cond_divflow_vert)))).*current_time_of_flight(part_inds_toupdate(~cond_divflow_vert))' + P_entry_z(part_inds_toupdate(~cond_divflow_vert));
            end
            
            ending_in_element(part_inds_toupdate) = true;
            cumulative_time_jj(rel_indices(successfull_relIDs)) = old_cumulative_time_jj(rel_indices(successfull_relIDs)) + (current_time_of_flight(successfull_relIDs));
            
        end
    end
    
    
    cond_for_max_simulation_time = cumulative_time_jj(rel_indices(successfull_relIDs))> max_simulation_time;
    
    if sum(cond_for_max_simulation_time)>0
        part_inds_max_sim_exceeded = find(1==cond_for_max_simulation_time)';
        
        current_time_of_flight(part_inds_max_sim_exceeded) = max_simulation_time - time_of_flight_so_far(rel_indices(part_inds_max_sim_exceeded));
        
        
        cond_divflow_horiz = divflow(cond_for_max_simulation_time)';
        
        if sum(cond_divflow_horiz)>0
            exit_point(1,part_inds_max_sim_exceeded(cond_divflow_horiz)) = ((a_xmat(sub2ind(size_a_b_mat,element_index(part_inds_max_sim_exceeded(cond_divflow_horiz)),transient_indices(part_inds_max_sim_exceeded(cond_divflow_horiz)))) + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(part_inds_max_sim_exceeded(cond_divflow_horiz)),transient_indices(part_inds_max_sim_exceeded(cond_divflow_horiz)))).*P_entry_x(part_inds_max_sim_exceeded(cond_divflow_horiz))).*exp(b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(part_inds_max_sim_exceeded(cond_divflow_horiz)),transient_indices(part_inds_max_sim_exceeded(cond_divflow_horiz)))).*current_time_of_flight(part_inds_max_sim_exceeded(cond_divflow_horiz))') - a_xmat(sub2ind(size_a_b_mat,element_index(part_inds_max_sim_exceeded(cond_divflow_horiz)),transient_indices(part_inds_max_sim_exceeded(cond_divflow_horiz)))))./b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(part_inds_max_sim_exceeded(cond_divflow_horiz)),transient_indices(part_inds_max_sim_exceeded(cond_divflow_horiz))));
            
            exit_point(2,part_inds_max_sim_exceeded(cond_divflow_horiz)) = ((a_ymat(sub2ind(size_a_b_mat,element_index(part_inds_max_sim_exceeded(cond_divflow_horiz)),transient_indices(part_inds_max_sim_exceeded(cond_divflow_horiz)))) + b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(part_inds_max_sim_exceeded(cond_divflow_horiz)),transient_indices(part_inds_max_sim_exceeded(cond_divflow_horiz)))).*P_entry_y(part_inds_max_sim_exceeded(cond_divflow_horiz))).*exp(b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(part_inds_max_sim_exceeded(cond_divflow_horiz)),transient_indices(part_inds_max_sim_exceeded(cond_divflow_horiz)))).*current_time_of_flight(part_inds_max_sim_exceeded(cond_divflow_horiz))') - a_ymat(sub2ind(size_a_b_mat,element_index(part_inds_max_sim_exceeded(cond_divflow_horiz)),transient_indices(part_inds_max_sim_exceeded(cond_divflow_horiz)))))./b_div_horiz_mat(sub2ind(size_a_b_mat,element_index(part_inds_max_sim_exceeded(cond_divflow_horiz)),transient_indices(part_inds_max_sim_exceeded(cond_divflow_horiz))));
        end
        
        if sum(~cond_divflow_horiz)>0
            exit_point(1,part_inds_max_sim_exceeded(~cond_divflow_horiz)) = a_xmat(sub2ind(size_a_b_mat,element_index(part_inds_max_sim_exceeded(~cond_divflow_horiz)),transient_indices(part_inds_max_sim_exceeded(~cond_divflow_horiz)))).*current_time_of_flight(part_inds_max_sim_exceeded(~cond_divflow_horiz))' + P_entry_x(part_inds_max_sim_exceeded(~cond_divflow_horiz));
            exit_point(2,part_inds_max_sim_exceeded(~cond_divflow_horiz)) = a_ymat(sub2ind(size_a_b_mat,element_index(part_inds_max_sim_exceeded(~cond_divflow_horiz)),transient_indices(part_inds_max_sim_exceeded(~cond_divflow_horiz)))).*current_time_of_flight(part_inds_max_sim_exceeded(~cond_divflow_horiz))' + P_entry_y(part_inds_max_sim_exceeded(~cond_divflow_horiz));
        end
        
        
        cond_divflow_vert = divflow_v(cond_for_max_simulation_time)';
        
        if sum(cond_divflow_vert)>0
            exit_point(3,part_inds_max_sim_exceeded(cond_divflow_vert)) = ((a_zmat(sub2ind(size_a_b_mat,element_index(part_inds_max_sim_exceeded(cond_divflow_vert)),transient_indices(part_inds_max_sim_exceeded(cond_divflow_vert)))) + b_div_vert_mat(sub2ind(size_a_b_mat,element_index(part_inds_max_sim_exceeded(cond_divflow_vert)),transient_indices(part_inds_max_sim_exceeded(cond_divflow_vert)))).*P_entry_z(part_inds_max_sim_exceeded(cond_divflow_vert))).*exp(b_div_vert_mat(sub2ind(size_a_b_mat,element_index(part_inds_max_sim_exceeded(cond_divflow_vert)),transient_indices(part_inds_max_sim_exceeded(cond_divflow_vert)))).*current_time_of_flight(part_inds_max_sim_exceeded(cond_divflow_vert))') - a_zmat(sub2ind(size_a_b_mat,element_index(part_inds_max_sim_exceeded(cond_divflow_vert)),transient_indices(part_inds_max_sim_exceeded(cond_divflow_vert)))))./b_div_vert_mat(sub2ind(size_a_b_mat,element_index(part_inds_max_sim_exceeded(cond_divflow_vert)),transient_indices(part_inds_max_sim_exceeded(cond_divflow_vert))));
        end
        
        if sum(~cond_divflow_vert)>0
            exit_point(3,part_inds_max_sim_exceeded(~cond_divflow_vert)) = a_zmat(sub2ind(size_a_b_mat,element_index(part_inds_max_sim_exceeded(~cond_divflow_vert)),transient_indices(part_inds_max_sim_exceeded(~cond_divflow_vert)))).*current_time_of_flight(part_inds_max_sim_exceeded(~cond_divflow_vert))' + P_entry_z(part_inds_max_sim_exceeded(~cond_divflow_vert));
        end
        
        
        ending_in_element(part_inds_max_sim_exceeded) = true;
        cumulative_time_jj(rel_indices(successfull_relIDs)) = old_cumulative_time_jj(rel_indices(successfull_relIDs)) + (current_time_of_flight(successfull_relIDs));
        
        time_is_exceeded(part_inds_max_sim_exceeded) = true;
        
    end
    
    
    if transient
        
        if sum(cond_for_transient_velocity_update)>0
            transient_indices(part_inds_toupdate) = transient_indices(part_inds_toupdate) + 1;
        end
        
    end
    
    
    time_of_flight_so_far(rel_indices(successfull_relIDs)) = time_of_flight_so_far(rel_indices(successfull_relIDs)) + current_time_of_flight(successfull_relIDs)';
    
    P_entry_x = exit_point(1,successfull_relIDs);
    P_entry_y = exit_point(2,successfull_relIDs);
    P_entry_z = exit_point(3,successfull_relIDs);
    
    %% save horizontal exit location/ time:
    
    if on_gpu
        exit_pint = real(exit_point);
    end
    
    matrix_of_all_exit_points_x(while_iter+1,rel_indices(successfull_relIDs)) = exit_point(1,successfull_relIDs); % HINT FOR REWRITING: ONLY ASSIGN SUCCESSFULL IDs
    matrix_of_all_exit_points_y(while_iter+1,rel_indices(successfull_relIDs)) = exit_point(2,successfull_relIDs);
    
    % for saving: local to global coordinates
    
    if sum(exit_on_lower_layer)>0
        where_exit_on_lower_layer = find(1==exit_on_lower_layer)';
        exit_point(3,where_exit_on_lower_layer) = 0;
    end
    
    if sum(exit_on_upper_layer)>0
        where_exit_on_upper_layer = find(1==exit_on_upper_layer)';
        exit_point(3,where_exit_on_upper_layer) = 1;
    end
    
    exit_point_glob3 = loc_to_glob_coord_mixC_mat(Xpri,exit_point,A_horiz(element_index),no_active_particles,on_gpu);
    
    exit_point_glob = exit_point_glob3(:,3);
    matrix_of_all_exit_points_z(while_iter+1,rel_indices(successfull_relIDs)) = exit_point_glob(successfull_relIDs);
    
    
    matrix_of_all_flight_times(while_iter,rel_indices(successfull_relIDs)) = current_time_of_flight(successfull_relIDs)';
    matrix_of_all_element_IDs(while_iter,rel_indices(successfull_relIDs)) = element_index(successfull_relIDs)';
    matrix_of_all_zone_IDs(while_iter,rel_indices(successfull_relIDs)) = zone_for_elements(element_index(successfull_relIDs))';
    
    while_iter = while_iter + 1;
    
    
    %% Determine next element and exclude current or then (past) entry/exit face:
    tracking_finished = true(1,no_active_particles);
    
    if sum(~ending_in_element)>0
        
        % Exclude old exit-face as new exit-face:
        
        where_right_index_1 = find(1==right_index)';
        where_right_index_2 = find(2==right_index)';
        where_right_index_3 = find(3==right_index)';
        where_right_index_4 = find(4==right_index)';
        where_right_index_5 = find(5==right_index)';
        
        len_index_4 = length(where_right_index_4);
        len_index_5 = length(where_right_index_5);
        
        if ~isempty(where_right_index_1)
            last_Lf_ending_nodes(where_right_index_1,:) = [node1(where_right_index_1),node2(where_right_index_1),node4(where_right_index_1),node5(where_right_index_1)];
        end
        if ~isempty(where_right_index_2)
            last_Lf_ending_nodes(where_right_index_2,:) = [node1(where_right_index_2),node3(where_right_index_2),node4(where_right_index_2),node6(where_right_index_2)];
        end
        if ~isempty(where_right_index_3)
            last_Lf_ending_nodes(where_right_index_3,:) = [node2(where_right_index_3),node3(where_right_index_3),node5(where_right_index_3),node6(where_right_index_3)];
        end
        if ~isempty(where_right_index_4)
            last_Lf_ending_nodes(where_right_index_4,:) = [node1(where_right_index_4),node2(where_right_index_4),node3(where_right_index_4),zeros(len_index_4,1)];
        end
        if ~isempty(where_right_index_5)
            last_Lf_ending_nodes(where_right_index_5,:) = [node4(where_right_index_5),node5(where_right_index_5),node6(where_right_index_5),zeros(len_index_5,1)];
        end
        
        % Determine next element,
        current_element_of_part_vor=current_element_of_part_i;
        current_tracking_neighb = all_neighb_indices(element_index,:);
        
        cond_neighb1_exit = current_tracking_neighb(:,1)~=0;
        if sum(cond_neighb1_exit)>0
            where_neighb1_exist = find(1==cond_neighb1_exit);
            
            logical_comp_nodes1 = last_Lf_ending_nodes(where_neighb1_exist,1) == node_ids(current_tracking_neighb(where_neighb1_exist,1),1:3); % Change add indices
            logical_comp_nodes2 = last_Lf_ending_nodes(where_neighb1_exist,2) == node_ids(current_tracking_neighb(where_neighb1_exist,1),1:3);
            logical_comp_nodes3 = last_Lf_ending_nodes(where_neighb1_exist,3) == node_ids(current_tracking_neighb(where_neighb1_exist,1),4:6);
            logical_comp_nodes4 = last_Lf_ending_nodes(where_neighb1_exist,4) == node_ids(current_tracking_neighb(where_neighb1_exist,1),4:6);
            
            next_elem_is_neighb1_sum = sum(logical_comp_nodes1,2) + sum(logical_comp_nodes2,2) + sum(logical_comp_nodes3,2) + sum(logical_comp_nodes4,2);
            next_elem_is_neighb1_logic = next_elem_is_neighb1_sum;
            next_elem_is_neighb1_logic(next_elem_is_neighb1_logic~=4) = 0;
            next_elem_is_neighb1_logic(next_elem_is_neighb1_logic==4) = 1;
            
            if sum(next_elem_is_neighb1_logic)>0
                which_elems_enter_first_neighb = find(1==next_elem_is_neighb1_logic);
                current_element_of_part_i(where_neighb1_exist(which_elems_enter_first_neighb)) = current_tracking_neighb(where_neighb1_exist(which_elems_enter_first_neighb),1);
                exit_horizontal(where_neighb1_exist(which_elems_enter_first_neighb)) = true;
                exit_on_top(where_neighb1_exist(which_elems_enter_first_neighb)) = false;
                exit_on_bottom(where_neighb1_exist(which_elems_enter_first_neighb)) = false;
                tracking_finished(where_neighb1_exist(which_elems_enter_first_neighb))  = false;
            end
        end
        
        
        
        cond_neighb2_exit = current_tracking_neighb(:,2)~=0;
        if sum(cond_neighb2_exit)>0
            where_neighb2_exist = find(1==cond_neighb2_exit);
            
            logical_comp_nodes1 = last_Lf_ending_nodes(where_neighb2_exist,1) == node_ids(current_tracking_neighb(where_neighb2_exist,2),1:3);
            logical_comp_nodes2 = last_Lf_ending_nodes(where_neighb2_exist,2) == node_ids(current_tracking_neighb(where_neighb2_exist,2),1:3);
            logical_comp_nodes3 = last_Lf_ending_nodes(where_neighb2_exist,3) == node_ids(current_tracking_neighb(where_neighb2_exist,2),4:6);
            logical_comp_nodes4 = last_Lf_ending_nodes(where_neighb2_exist,4) == node_ids(current_tracking_neighb(where_neighb2_exist,2),4:6);
            
            next_elem_is_neighb2_sum = sum(logical_comp_nodes1,2) + sum(logical_comp_nodes2,2) + sum(logical_comp_nodes3,2) + sum(logical_comp_nodes4,2);
            next_elem_is_neighb2_logic = next_elem_is_neighb2_sum;
            next_elem_is_neighb2_logic(next_elem_is_neighb2_logic~=4) = 0;
            next_elem_is_neighb2_logic(next_elem_is_neighb2_logic==4) = 1;
            
            if sum(next_elem_is_neighb2_logic)>0
                which_elems_enter_second_neighb = find(1==next_elem_is_neighb2_logic);
                
                current_element_of_part_i(where_neighb2_exist(which_elems_enter_second_neighb)) = current_tracking_neighb(where_neighb2_exist(which_elems_enter_second_neighb),2);
                exit_horizontal(where_neighb2_exist(which_elems_enter_second_neighb)) = true;
                exit_on_top(where_neighb2_exist(which_elems_enter_second_neighb)) = false;
                exit_on_bottom(where_neighb2_exist(which_elems_enter_second_neighb)) = false;
                tracking_finished(where_neighb2_exist(which_elems_enter_second_neighb))  = false;
            end
        end
        
        
        cond_neighb3_exit = current_tracking_neighb(:,3)~=0;
        if sum(cond_neighb3_exit)>0
            where_neighb3_exist = find(1==cond_neighb3_exit);
            
            logical_comp_nodes1 = last_Lf_ending_nodes(where_neighb3_exist,1) == node_ids(current_tracking_neighb(where_neighb3_exist,3),1:3);
            logical_comp_nodes2 = last_Lf_ending_nodes(where_neighb3_exist,2) == node_ids(current_tracking_neighb(where_neighb3_exist,3),1:3);
            logical_comp_nodes3 = last_Lf_ending_nodes(where_neighb3_exist,3) == node_ids(current_tracking_neighb(where_neighb3_exist,3),4:6);
            logical_comp_nodes4 = last_Lf_ending_nodes(where_neighb3_exist,4) == node_ids(current_tracking_neighb(where_neighb3_exist,3),4:6);
            
            next_elem_is_neighb3_sum = sum(logical_comp_nodes1,2) + sum(logical_comp_nodes2,2) + sum(logical_comp_nodes3,2) + sum(logical_comp_nodes4,2);
            next_elem_is_neighb3_logic = next_elem_is_neighb3_sum;
            next_elem_is_neighb3_logic(next_elem_is_neighb3_logic~=4) = 0;
            next_elem_is_neighb3_logic(next_elem_is_neighb3_logic==4) = 1;
            
            if sum(next_elem_is_neighb3_logic)>0
                which_elems_enter_third_neighb = find(1==next_elem_is_neighb3_logic);
                current_element_of_part_i(where_neighb3_exist(which_elems_enter_third_neighb)) = current_tracking_neighb(where_neighb3_exist(which_elems_enter_third_neighb),3);
                exit_horizontal(where_neighb3_exist(which_elems_enter_third_neighb)) = true;
                exit_on_top(where_neighb3_exist(which_elems_enter_third_neighb)) = false;
                exit_on_bottom(where_neighb3_exist(which_elems_enter_third_neighb)) = false;
                tracking_finished(where_neighb3_exist(which_elems_enter_third_neighb))  = false;
            end
        end
        
        
        
        cond_neighb4_exit = current_tracking_neighb(:,4)~=0;
        if sum(cond_neighb4_exit)>0
            where_neighb4_exist = find(1==cond_neighb4_exit);
            
            logical_comp_nodes1 = last_Lf_ending_nodes(where_neighb4_exist,1) == node_ids(current_tracking_neighb(where_neighb4_exist,4),4:6);
            logical_comp_nodes2 = last_Lf_ending_nodes(where_neighb4_exist,2) == node_ids(current_tracking_neighb(where_neighb4_exist,4),4:6);
            logical_comp_nodes3 = last_Lf_ending_nodes(where_neighb4_exist,3) == node_ids(current_tracking_neighb(where_neighb4_exist,4),4:6);
            
            next_elem_is_neighb4_sum = sum(logical_comp_nodes1,2) + sum(logical_comp_nodes2,2) + sum(logical_comp_nodes3,2);
            next_elem_is_neighb4_logic = next_elem_is_neighb4_sum;
            next_elem_is_neighb4_logic(next_elem_is_neighb4_logic~=3) = 0;
            next_elem_is_neighb4_logic(next_elem_is_neighb4_logic==3) = 1;
            
            if sum(next_elem_is_neighb4_logic)>0
                
                which_elems_enter_fourth_neighb = find(1==next_elem_is_neighb4_logic);
                
                current_element_of_part_i(where_neighb4_exist(which_elems_enter_fourth_neighb)) = current_tracking_neighb(where_neighb4_exist(which_elems_enter_fourth_neighb),4);
                exit_horizontal(where_neighb4_exist(which_elems_enter_fourth_neighb)) = true;
                exit_on_top(where_neighb4_exist(which_elems_enter_fourth_neighb)) = false;
                exit_on_bottom(where_neighb4_exist(which_elems_enter_fourth_neighb)) = false;
                tracking_finished(where_neighb4_exist(which_elems_enter_fourth_neighb))  = false;
                
            end
            
        end
        
        
        cond_neighb5_exit = current_tracking_neighb(:,5)~=0;
        
        if sum(cond_neighb5_exit)>0
            where_neighb5_exist = find(1==cond_neighb5_exit);
            
            logical_comp_nodes1 = last_Lf_ending_nodes(where_neighb5_exist,1) == node_ids(current_tracking_neighb(where_neighb5_exist,5),1:3);
            logical_comp_nodes2 = last_Lf_ending_nodes(where_neighb5_exist,2) == node_ids(current_tracking_neighb(where_neighb5_exist,5),1:3);
            logical_comp_nodes3 = last_Lf_ending_nodes(where_neighb5_exist,3) == node_ids(current_tracking_neighb(where_neighb5_exist,5),1:3);
            
            next_elem_is_neighb5_sum = sum(logical_comp_nodes1,2) + sum(logical_comp_nodes2,2) + sum(logical_comp_nodes3,2);
            next_elem_is_neighb5_logic = next_elem_is_neighb5_sum;
            next_elem_is_neighb5_logic(next_elem_is_neighb5_logic~=3) = 0;
            next_elem_is_neighb5_logic(next_elem_is_neighb5_logic==3) = 1;
            
            if sum(next_elem_is_neighb5_logic)>0
                
                which_elems_enter_fifth_neighb = find(1==next_elem_is_neighb5_logic);
                
                current_element_of_part_i(where_neighb5_exist(which_elems_enter_fifth_neighb)) = current_tracking_neighb(where_neighb5_exist(which_elems_enter_fifth_neighb),5);
                exit_horizontal(where_neighb5_exist(which_elems_enter_fifth_neighb)) = true;
                exit_on_top(where_neighb5_exist(which_elems_enter_fifth_neighb)) = false;
                exit_on_bottom(where_neighb5_exist(which_elems_enter_fifth_neighb)) = false;
                tracking_finished(where_neighb5_exist(which_elems_enter_fifth_neighb))  = false;
                
            end
            
        end
        
        %% TEST if exit point is on a node in the horizontal:
        
        cond_for_exit_on_node = ((exit_point(1,:)'==xA | exit_point(1,:)'==xB | exit_point(1,:)'==xC) & (exit_point(2,:)'==yA | exit_point(2,:)'==yB | exit_point(2,:)'==yC)) & ~(ending_in_element');
        sum_cond_for_exit_on_node = sum(cond_for_exit_on_node);
        
        
        if sum_cond_for_exit_on_node>0
            where_exit_on_node = find(1==cond_for_exit_on_node);
            
            next_ele = current_element_of_part_i(where_exit_on_node)';
            
            node1_n = node_ids(next_ele,1);
            node2_n = node_ids(next_ele,2);
            node3_n = node_ids(next_ele,3);
            node4_n = node_ids(next_ele,4);
            node5_n = node_ids(next_ele,5);
            node6_n = node_ids(next_ele,6);
            
            xA_n = coordinates(node1_n,1);
            yA_n = coordinates(node1_n,2);
            
            xB_n = coordinates(node2_n,1);
            yB_n = coordinates(node2_n,2);
            
            xC_n = coordinates(node3_n,1);
            yC_n = coordinates(node3_n,2);
            
            z_mean_bot_n = (coordinates(node1_n,3) + coordinates(node2_n,3) + coordinates(node3_n,3))/3;
            z_mean_top_n = (coordinates(node4_n,3) + coordinates(node5_n,3) + coordinates(node6_n,3))/3;
            
            cond_x = exit_point(1,where_exit_on_node)'== [xA_n,xB_n,xC_n];
            cond_y = exit_point(2,where_exit_on_node)'== [yA_n,yB_n,yC_n];
            cond_xy = cond_x & cond_y;
            cond_xy = double(cond_xy);
            cond_xy(:,1) = cond_xy(:,1)*1;
            cond_xy(:,2) = cond_xy(:,2)*2;
            cond_xy(:,3) = cond_xy(:,3)*3;
            sum_cond_xy = sum(cond_xy,2);
            
            on_node_i = sum_cond_xy;
            
            intmed_help = 6 - on_node_i;
            not_on_nodes_jj = zeros(sum_cond_for_exit_on_node,2);
            
            find_5 = find(5==intmed_help);
            find_4 = find(4==intmed_help);
            find_3 = find(3==intmed_help);
            
            if ~isempty(find_5)
                not_on_nodes_jj(find_5,:) = [2,3];
            end
            if ~isempty(find_4)
                not_on_nodes_jj(find_4,:) = [1,3];
            end
            if ~isempty(find_3)
                not_on_nodes_jj(find_3,:) = [1,2];
            end
            
            
            lambda_pos = 1 - 2e6*eps;
            lambda_neg = 1e6*eps;
            
            coord_x = [xA_n;xB_n;xC_n];
            coord_y = [yA_n;yB_n;yC_n];
            
            lin_mulitplier_inds = 0:3:(length(on_node_i)-1)*3;
            lin_on_node_i = on_node_i + lin_mulitplier_inds';
            
            not_on_nodes_jj_first = not_on_nodes_jj(:,1);
            not_on_nodes_jj_second = not_on_nodes_jj(:,2);
            
            lin_not_on_nodes_jj_first = not_on_nodes_jj_first + lin_mulitplier_inds';
            lin_not_on_nodes_jj_second =  not_on_nodes_jj_second + lin_mulitplier_inds';
            
            exit_point(1,where_exit_on_node) = lambda_pos*coord_x(lin_on_node_i)' + lambda_neg*coord_x(lin_not_on_nodes_jj_first)' + lambda_neg*coord_x(lin_not_on_nodes_jj_second)';
            exit_point(2,where_exit_on_node) = lambda_pos*coord_y(lin_on_node_i)' + lambda_neg*coord_y(lin_not_on_nodes_jj_first)' + lambda_neg*coord_y(lin_not_on_nodes_jj_second)';
            
            assoc_exit_points_vert = exit_point(3,where_exit_on_node);
            exit_down = assoc_exit_points_vert<10^-(10);
            exit_top  = assoc_exit_points_vert> 1 - 10^-(10);
            
            if sum(exit_down)>0
                where_exactly = find(1==exit_down);
                exit_point(3,where_exit_on_node(where_exactly)) = exit_point(3,where_exit_on_node(where_exactly)) + 10^-8*(z_mean_bot_n(where_exactly) - exit_point(3,where_exit_on_node(where_exactly)));
            end
            if sum(exit_top)>0
                where_exactly = find(1==exit_top);
                exit_point(3,where_exit_on_node(where_exactly)) = exit_point(3,where_exit_on_node(where_exactly)) + 10^-8*(z_mean_top_n(where_exactly) - exit_point(3,where_exit_on_node(where_exactly)));
            end
            
            matrix_of_all_exit_points_x(while_iter,rel_indices(where_exit_on_node)) = exit_point(1,where_exit_on_node);
            matrix_of_all_exit_points_y(while_iter,rel_indices(where_exit_on_node)) = exit_point(2,where_exit_on_node);
            
            P_entry_x(where_exit_on_node) = exit_point(1,where_exit_on_node);
            P_entry_y(where_exit_on_node) = exit_point(2,where_exit_on_node);
        end
        
    end
    
    if  sum(ending_in_element)>0
        
        
        where_ending_in_element = find(1==ending_in_element);
        
        if sum(cond_for_max_simulation_time)>0
            inds_max_sim_exceeded = part_inds_max_sim_exceeded;
        else
            inds_max_sim_exceeded = [];
        end
        
        last_Lf_ending_nodes(where_ending_in_element,:) = 0;
        
        exit_horizontal(where_ending_in_element) = false;
        exit_on_top(where_ending_in_element) = false;
        exit_on_bottom(where_ending_in_element) = false;
        
        current_element_of_part_i(where_ending_in_element) = element_index(where_ending_in_element);
        
        tracking_finished(where_ending_in_element)  = false;
        tracking_finished(inds_max_sim_exceeded) = true;
        
    end
    
    
    
    if ~isempty(terminatedWell_and_failed_relIDs)
        tracking_finished(terminatedWell_and_failed_relIDs) = true;
    end
    
    
end

%% Save it all after the while loop

for jj = 1:no_particles
    
    if ~on_gpu
        x_save_fin = matrix_of_all_exit_points_x(:,jj);
        y_save_fin = matrix_of_all_exit_points_y(:,jj);
        z_save_fin = matrix_of_all_exit_points_z(:,jj);
        
        times_of_flight_fin = matrix_of_all_flight_times(:,jj);
        all_elem_ids_of_flight = matrix_of_all_element_IDs(:,jj);
        all_zone_ids_of_flight = matrix_of_all_zone_IDs(:,jj);
    else
        x_save_fin = gather(matrix_of_all_exit_points_x(:,jj));
        y_save_fin = gather(matrix_of_all_exit_points_y(:,jj));
        z_save_fin = gather(matrix_of_all_exit_points_z(:,jj));
        
        times_of_flight_fin = gather(matrix_of_all_flight_times(:,jj));
        all_elem_ids_of_flight = gather(matrix_of_all_element_IDs(:,jj));
        all_zone_ids_of_flight = gather(matrix_of_all_zone_IDs(:,jj));
        
    end
    
    x_save_fin(isnan(x_save_fin)) = [];
    y_save_fin(isnan(y_save_fin)) = [];
    z_save_fin(isnan(z_save_fin)) = [];
    times_of_flight_fin(isnan(times_of_flight_fin)) = [];
    all_elem_ids_of_flight(isnan(all_elem_ids_of_flight)) = [];
    all_zone_ids_of_flight(isnan(all_zone_ids_of_flight)) = [];
    
    cell_PT_x{jj} = x_save_fin;
    cell_PT_y{jj} = y_save_fin;
    cell_PT_z{jj} = z_save_fin;
    cell_PT_t{jj} = times_of_flight_fin;
    cell_PT_elemIDs{jj} = all_elem_ids_of_flight;
    cell_PT_zoneIDs{jj} = all_zone_ids_of_flight;
    
    
end
toc
%%

x_coord_part_name = strcat('x_coord_PT_',prefix,'.mat');
save(x_coord_part_name,'cell_PT_x')

y_coord_part_name = strcat('y_coord_PT_',prefix,'.mat');
save(y_coord_part_name,'cell_PT_y')

z_coord_part_name = strcat('z_coord_PT_',prefix,'.mat');
save(z_coord_part_name,'cell_PT_z')

travel_times_part_name = strcat('travel_times_PT_',prefix,'.mat');
save(travel_times_part_name,'cell_PT_t')

element_IDs_name = strcat('element_IDs_PT_',prefix,'.mat');
save(element_IDs_name,'cell_PT_elemIDs')

zone_IDs_name = strcat('zone_IDs_PT_',prefix,'.mat');
save(zone_IDs_name,'cell_PT_zoneIDs')


if plotting
    figure(2)
    
    hold on
    
    for i = 1:elem_per_layer
        
        node_ids_sh1 = node_ids(i,1);
        node_ids_sh2 = node_ids(i,2);
        node_ids_sh3 = node_ids(i,3);
        
        xA = coordinates(node_ids_sh1,1);
        yA = coordinates(node_ids_sh1,2);
        zA = coordinates(node_ids_sh1,3);
        xB = coordinates(node_ids_sh2,1);
        yB = coordinates(node_ids_sh2,2);
        zB = coordinates(node_ids_sh2,3);
        xC = coordinates(node_ids_sh3,1);
        yC = coordinates(node_ids_sh3,2);
        zC = coordinates(node_ids_sh3,3);
        
        line([xA xB xC xA],[yA yB yC yA],[zA zB zC zA],'Color',[0.6 0.6 0.6])
        
    end
    
    xlabel('x [m]');ylabel('y [m]');
    
    frame_plot_x = (max(coordinates(:,1)) - min(coordinates(:,1)))/20;
    frame_plot_y = (max(coordinates(:,2)) - min(coordinates(:,2)))/20;
    
    axis([min(coordinates(:,1))-frame_plot_x max(coordinates(:,1))+frame_plot_x min(coordinates(:,2))-frame_plot_y max(coordinates(:,2))+frame_plot_y]);
    
    daspect([1 1 1]);
    
    set(gca,'FontSize',12)
    
    for nn = 1:no_particles
        
        x_plot = cell_PT_x{nn};
        y_plot = cell_PT_y{nn};
        z_plot = cell_PT_z{nn};
        
        plot(x_plot,y_plot,'r','LineWidth',1.15);
        
    end
    
    hold off
end


disp('ready!')
