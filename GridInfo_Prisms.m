%% Metainformation:
%
% Author:           Philipp Selzer
% Version:          2.0
% Date:             20.11.2019
% Last update:      24.03.2022
% Purpose:          This code computes all necessary geometrical informations
%                   for the FVM-flux-correction and the particle tracking
%                   scheme on prism. It has to be run once per grid.

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


%% ===================================================================
% Read all relevant HGS output (and input) files
%% ===================================================================

% read the binary node ids of the elements
n_dim = 3;
filename1 = strcat(prefix,'o.elements_pm');

no_nodes_per_element = 6;
no_elem =138565; % TO ADAPT
no_layers = 35; % TO ADAPT
elem_per_layer = no_elem/no_layers;
elem_top_layer = elem_per_layer; %

[node_ids] = read_HGS_binary_elements(no_nodes_per_element,no_elem,filename1);

% read the binary coordinates
filename2 = strcat(prefix,'o.coordinates_pm');

no_nodes = 74412; % TO ADAPT LINE
[coordinates] = read_HGS_binary_coordinates(n_dim,no_nodes,filename2);


%% Dirichlet boundaries: TO ADAPT SECTION

% % extract all nodes which are part of the Dirichlet boundary
filename_lbc1 = strcat(prefix,'o.node_set.bc_ammer_in');
filename_lbc2 = strcat(prefix,'o.node_set.bc_ammer_out');
filename_lbc3 = strcat(prefix,'o.node_set.bc_neckar');

inds_lbc1 = importdata(filename_lbc1);
inds_lbc1 = inds_lbc1(2:end);

inds_lbc2 = importdata(filename_lbc2);
inds_lbc2 = inds_lbc2(2:end);
%
inds_lbc3 = importdata(filename_lbc3);
inds_lbc3 = inds_lbc3(2:end);

diri_nodes = [inds_lbc1;inds_lbc2;inds_lbc3];
diri_nodes = unique(diri_nodes);


%% Neumann unlike zero nodes (i.e., flow boundary condition)
neumann_unlike0_exists = false; % TO ADAPT LINE

neumann_nodes = [];

if neumann_unlike0_exists
    filename_lbcN1 = strcat(prefix,'o.node_set.Neumannlateral_bc2'); % TO ADAPT LINE
    
    inds_lbcN2 = importdata(filename_lbcN1);
    inds_lbcN2 = inds_lbcN2(2:end);
    
    neumann_nodes = [inds_lbcN2];
    neumann_nodes = unique(neumann_nodes);
    
end

diri_neumann_nodes = [diri_nodes,neumann_nodes];

%% End read all relevant HGS output (and input) files

% ===================================================================
% Preallocations

no_poss_neighbs = 5; % number of possible neighbouring elements

vector_of_d_bounds_mult = zeros(no_elem,1,'uint32');
vector_of_n_bounds_mult = zeros(no_elem,1,'uint32');
face_index_d_mult = zeros(no_elem,1,'uint32');
face_index_n_mult = zeros(no_elem,1,'uint32');
d_mult_counter = 1;
n_mult_counter = 1;

face_index_mat = zeros(no_elem,no_poss_neighbs,'uint32');
face_flag_D_N = zeros(no_elem,no_poss_neighbs,'uint32');

centro_of_elem = zeros(no_elem,3);

tetrahedral_nodes = zeros(no_elem*3,4);
counter = 1;

Vol_prisms = zeros(no_elem,1);

for i = 1:no_elem
    
    node_ids1 = node_ids(i,1);
    node_ids2 = node_ids(i,2);
    node_ids3 = node_ids(i,3);
    node_ids4 = node_ids(i,4);
    node_ids5 = node_ids(i,5);
    node_ids6 = node_ids(i,6);
    
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
    
    d_nodes = intersect(node_ids(i,:),diri_nodes);
    
    n_nodes = intersect(node_ids(i,:),neumann_nodes);
    
    % 1st neighbour: 2,3,5,6
    % 2nd neighbour: 1,3,4,6
    % 3rd neighbour: 1,2,4,5
    % 4th neighbour: 1,2,3
    % 5th neighbour: 4,5,6
    counter_more_than4 = 1;
    if (length(d_nodes)==3 || length(d_nodes)==4)
        intsect_1st = intersect([node_ids2;node_ids3;node_ids5;node_ids6],d_nodes);
        intsect_2nd = intersect([node_ids1;node_ids3;node_ids4;node_ids6],d_nodes);
        intsect_3rd = intersect([node_ids1;node_ids2;node_ids4;node_ids5],d_nodes);
        intsect_4th = intersect([node_ids1;node_ids2;node_ids3],d_nodes);
        intsect_5th = intersect([node_ids4;node_ids5;node_ids6],d_nodes);
        
        if length(intsect_1st)==4
            face_index_mat(i,1) = 1;
            face_flag_D_N(i,1) = 1;
            face_index_d_mult(d_mult_counter) = 1;
            vector_of_d_bounds_mult(d_mult_counter) = i;
            d_mult_counter = d_mult_counter + 1;
        elseif length(intsect_2nd)==4
            face_index_mat(i,1) = 2;
            face_flag_D_N(i,1) = 1;
            face_index_d_mult(d_mult_counter) = 2;
            vector_of_d_bounds_mult(d_mult_counter) = i;
            d_mult_counter = d_mult_counter + 1;
        elseif length(intsect_3rd)==4
            face_index_mat(i,1) = 3;
            face_flag_D_N(i,1) = 1;
            face_index_d_mult(d_mult_counter) = 3;
            vector_of_d_bounds_mult(d_mult_counter) = i;
            d_mult_counter = d_mult_counter + 1;
        elseif length(intsect_4th)==3
            face_index_mat(i,1) = 4;
            face_flag_D_N(i,1) = 1;
            face_index_d_mult(d_mult_counter) = 4;
            vector_of_d_bounds_mult(d_mult_counter) = i;
            d_mult_counter = d_mult_counter + 1;
        elseif length(intsect_5th)==3
            face_index_mat(i,1) = 5;
            face_flag_D_N(i,1) = 1;
            face_index_d_mult(d_mult_counter) = 5;
            vector_of_d_bounds_mult(d_mult_counter) = i;
            d_mult_counter = d_mult_counter + 1;
        end
        
        
    elseif length(d_nodes)>4
        intsect_1st = intersect([node_ids2;node_ids3;node_ids5;node_ids6],d_nodes);
        intsect_2nd = intersect([node_ids1;node_ids3;node_ids4;node_ids6],d_nodes);
        intsect_3rd = intersect([node_ids1;node_ids2;node_ids4;node_ids5],d_nodes);
        intsect_4th = intersect([node_ids1;node_ids2;node_ids3],d_nodes);
        intsect_5th = intersect([node_ids4;node_ids5;node_ids6],d_nodes);
        
        if length(intsect_1st)==4
            face_index_mat(i,counter_more_than4) = 1;
            face_flag_D_N(i,counter_more_than4) = 1;
            counter_more_than4 = counter_more_than4 + 1 ;
            
            vector_of_d_bounds_mult(d_mult_counter) = i;
            face_index_d_mult(d_mult_counter) = 1;
            d_mult_counter = d_mult_counter + 1;
        end
        
        if length(intsect_2nd)==4
            face_index_mat(i,counter_more_than4) = 2;
            face_flag_D_N(i,counter_more_than4) = 1;
            counter_more_than4 = counter_more_than4 + 1 ;
            
            vector_of_d_bounds_mult(d_mult_counter) = i;
            face_index_d_mult(d_mult_counter) = 2;
            d_mult_counter = d_mult_counter + 1;
        end
        
        if length(intsect_3rd)==4
            face_index_mat(i,counter_more_than4) = 3;
            face_flag_D_N(i,counter_more_than4) = 1;
            counter_more_than4 = counter_more_than4 + 1 ;
            
            vector_of_d_bounds_mult(d_mult_counter) = i;
            face_index_d_mult(d_mult_counter) = 3;
            d_mult_counter = d_mult_counter + 1;
        end
        
        if length(intsect_4th)==3
            face_index_mat(i,counter_more_than4) = 4;
            face_flag_D_N(i,counter_more_than4) = 1;
            counter_more_than4 = counter_more_than4 + 1 ;
            
            vector_of_d_bounds_mult(d_mult_counter) = i;
            face_index_d_mult(d_mult_counter) = 4;
            d_mult_counter = d_mult_counter + 1;
        end
        
        if length(intsect_5th)==3
            face_index_mat(i,counter_more_than4) = 5;
            face_flag_D_N(i,counter_more_than4) = 1;
            counter_more_than4 = counter_more_than4 + 1 ;
            
            vector_of_d_bounds_mult(d_mult_counter) = i;
            face_index_d_mult(d_mult_counter) = 5;
            d_mult_counter = d_mult_counter + 1;
        end
        
    end
    
    if (length(n_nodes)==3 || length(n_nodes)==4)
        
        intsect_1st = intersect([node_ids2;node_ids3;node_ids5;node_ids6],n_nodes);
        intsect_2nd = intersect([node_ids1;node_ids3;node_ids4;node_ids6],n_nodes);
        intsect_3rd = intersect([node_ids1;node_ids2;node_ids4;node_ids5],n_nodes);
        intsect_4th = intersect([node_ids1;node_ids2;node_ids3],n_nodes);
        intsect_5th = intersect([node_ids4;node_ids5;node_ids6],n_nodes);
        
        if length(intsect_1st)==4
            face_index_mat(i,1) = 1;
            face_flag_D_N(i,1) = 2;
            face_index_n_mult(n_mult_counter) = 1;
            vector_of_n_bounds_mult(n_mult_counter) = i;
            n_mult_counter = n_mult_counter + 1;
        elseif length(intsect_2nd)==4
            face_index_mat(i,1) = 2;
            face_flag_D_N(i,1) = 2;
            face_index_n_mult(n_mult_counter) = 2;
            vector_of_n_bounds_mult(n_mult_counter) = i;
            n_mult_counter = n_mult_counter + 1;
        elseif length(intsect_3rd)==4
            face_index_mat(i,1) = 3;
            face_flag_D_N(i,1) = 2;
            face_index_n_mult(n_mult_counter) = 3;
            vector_of_n_bounds_mult(n_mult_counter) = i;
            n_mult_counter = n_mult_counter + 1;
        elseif length(intsect_4th)==3
            face_index_mat(i,1) = 4;
            face_flag_D_N(i,1) = 2;
            face_index_n_mult(n_mult_counter) = 4;
            vector_of_n_bounds_mult(n_mult_counter) = i;
            n_mult_counter = n_mult_counter + 1;
        elseif length(intsect_5th)==3
            face_index_mat(i,1) = 5;
            face_flag_D_N(i,1) = 2;
            face_index_n_mult(n_mult_counter) = 5;
            vector_of_n_bounds_mult(n_mult_counter) = i;
            n_mult_counter = n_mult_counter + 1;
        end
        
        
    elseif length(n_nodes)>4
        
        intsect_1st = intersect([node_ids2;node_ids3;node_ids5;node_ids6],d_nodes);
        intsect_2nd = intersect([node_ids1;node_ids3;node_ids4;node_ids6],d_nodes);
        intsect_3rd = intersect([node_ids1;node_ids2;node_ids4;node_ids5],d_nodes);
        intsect_4th = intersect([node_ids1;node_ids2;node_ids3],d_nodes);
        intsect_5th = intersect([node_ids4;node_ids5;node_ids6],d_nodes);
        
        if length(intsect_1st)==4
            face_index_mat(i,counter_more_than4) = 1;
            face_flag_D_N(i,counter_more_than4) = 2;
            counter_more_than4 = counter_more_than4 + 1 ;
            
            vector_of_n_bounds_mult(n_mult_counter) = i;
            face_index_n_mult(n_mult_counter) = 1;
            n_mult_counter = n_mult_counter + 1;
        end
        
        if length(intsect_2nd)==4
            face_index_mat(i,counter_more_than4) = 2;
            face_flag_D_N(i,counter_more_than4) = 2;
            counter_more_than4 = counter_more_than4 + 1 ;
            
            vector_of_n_bounds_mult(n_mult_counter) = i;
            face_index_n_mult(n_mult_counter) = 2;
            n_mult_counter = n_mult_counter + 1;
        end
        
        if length(intsect_3rd)==4
            face_index_mat(i,counter_more_than4) = 3;
            face_flag_D_N(i,counter_more_than4) = 2;
            counter_more_than4 = counter_more_than4 + 1 ;
            
            vector_of_n_bounds_mult(n_mult_counter) = i;
            face_index_n_mult(n_mult_counter) = 3;
            n_mult_counter = n_mult_counter + 1;
        end
        
        if length(intsect_4th)==3
            face_index_mat(i,counter_more_than4) = 4;
            face_flag_D_N(i,counter_more_than4) = 2;
            counter_more_than4 = counter_more_than4 + 1 ;
            
            vector_of_n_bounds_mult(n_mult_counter) = i;
            face_index_n_mult(n_mult_counter) = 4;
            n_mult_counter = n_mult_counter + 1;
        end
        
        if length(intsect_5th)==3
            face_index_mat(i,counter_more_than4) = 5;
            face_flag_D_N(i,counter_more_than4) = 2;
            counter_more_than4 = counter_more_than4 + 1 ;
            
            vector_of_n_bounds_mult(n_mult_counter) = i;
            face_index_n_mult(n_mult_counter) = 5;
            n_mult_counter = n_mult_counter + 1;
        end
        
        
    end
    
    a=[(xA + xB + xC)/3,(yA + yB + yC)/3,(zA+zB+zC+zD+zE+zF)/6];
    centro_of_elem(i,:)=a;
    
    
    %% Split prisms in tetrahedra for search for particle initial elements later on:
    
    pos_VI1 = find(min(node_ids(i,1:3))==node_ids(i,1:3));
    pos_VI2_3 = setdiff(1:1:3,pos_VI1);
    pos_VI4 = pos_VI1 + 3;
    pos_VI5_6 = pos_VI2_3 + 3;
    
    VI2_VI6 = [node_ids(i,pos_VI2_3(1)),node_ids(i,pos_VI5_6(2))];
    VI3_VI5 = [node_ids(i,pos_VI2_3(2)),node_ids(i,pos_VI5_6(1))];
    
    if min(VI2_VI6)<min(VI3_VI5)
        tetra_1 = [node_ids(i,pos_VI1),node_ids(i,pos_VI2_3),node_ids(i,pos_VI5_6(2))]; % VI_1_2_3_6
        tetra_2 = [node_ids(i,pos_VI1),node_ids(i,pos_VI2_3(1)),node_ids(i,pos_VI5_6(2)),node_ids(i,pos_VI5_6(1))]; % VI_1_2_6_5
        tetra_3 = [node_ids(i,pos_VI1),node_ids(i,pos_VI5_6),node_ids(i,pos_VI4)]; % VI_1_5_6_4
        
    elseif min(VI3_VI5)<min(VI2_VI6)
        tetra_1 = [node_ids(i,pos_VI1),node_ids(i,pos_VI2_3),node_ids(i,pos_VI5_6(1))]; % VI_1_2_3_5
        tetra_2 = [node_ids(i,pos_VI1),node_ids(i,pos_VI5_6(1)),node_ids(i,pos_VI2_3(2)),node_ids(i,pos_VI5_6(2))]; % VI_1_5_3_6
        tetra_3 = [node_ids(i,pos_VI1),node_ids(i,pos_VI5_6),node_ids(i,pos_VI4)]; % VI_1_5_6_4
    end
    
    tetrahedral_nodes(counter,:) = tetra_1;
    tetrahedral_nodes(counter+1,:) = tetra_2;
    tetrahedral_nodes(counter+2,:) = tetra_3;
    
    
    counter = counter+3;
    
    
end

tetrahedral_nodes_name = strcat('tetrahedral_nodes_',prefix,'.mat');
save(tetrahedral_nodes_name,'tetrahedral_nodes')

vector_of_d_bounds_mult(vector_of_d_bounds_mult==0) = [];
face_index_d_mult(face_index_d_mult==0) = [];
elem_diri_bound = unique(vector_of_d_bounds_mult);

vector_of_n_bounds_mult(vector_of_n_bounds_mult==0) = [];
face_index_n_mult(face_index_n_mult==0) = [];
elem_neumann_bound = unique(vector_of_n_bounds_mult);


%% Initializations

neighbours_of_all_elements = zeros(no_elem,5);

all_neighb_indices = zeros(no_elem,5);
A_f_all = zeros(no_elem,5);
A_quad = zeros(no_elem,3);
A_tri = zeros(no_elem,2);

l_elem_all = zeros(no_elem,no_poss_neighbs);
l_neighb_all = zeros(no_elem,no_poss_neighbs);

A_horiz = zeros(no_elem,1);
A_top = zeros(no_elem,1);
hQ_all = zeros(no_elem,n_dim);

Aniso_coeff = zeros(no_elem,no_poss_neighbs,3);

corr_fact_orth_save = zeros(no_elem,no_poss_neighbs);
diri_counter = 1;

active_DiriNeumann_bounds = [elem_diri_bound;elem_neumann_bound];

Vol_prisms2 = zeros(no_elem,1);

for j= 1:no_elem % element loop
    % ===================================================================
    % Element at Dirichlet-boundary:
    find_jDN = find(j==active_DiriNeumann_bounds);
    
    if ~isempty(find_jDN)
        flag_at_DiriNeumann = true;
    else
        flag_at_DiriNeumann = false;
    end
    
    % ===================================================================
    % Determination of all neighbouring elements:
    
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
    
    neighb_ele = zeros(5,1);
    ind_neighb_neu = zeros(1,no_poss_neighbs);
    corresp_nodes = zeros(n_dim +1,no_poss_neighbs);
    
    current_layer_add = floor((j-1)/elem_per_layer);
    
    if j<=elem_per_layer
        
        for i = 1:elem_per_layer
            
            % First possible lateral neighbour BC:
            if ((isequal(node_ids(i + (current_layer_add*elem_per_layer),1),node_ids2) || isequal(node_ids(i + (current_layer_add*elem_per_layer),2),node_ids2) || isequal(node_ids(i + (current_layer_add*elem_per_layer),3),node_ids2)) ...
                    && (isequal(node_ids(i + (current_layer_add*elem_per_layer),1),node_ids3) || isequal(node_ids(i + (current_layer_add*elem_per_layer),2),node_ids3) || isequal(node_ids(i + (current_layer_add*elem_per_layer),3),node_ids3))...
                    && (isequal(node_ids(i + (current_layer_add*elem_per_layer),4),node_ids5) || isequal(node_ids(i + (current_layer_add*elem_per_layer),5),node_ids5) || isequal(node_ids(i + (current_layer_add*elem_per_layer),6),node_ids5)) ...
                    && (isequal(node_ids(i + (current_layer_add*elem_per_layer),4),node_ids6) || isequal(node_ids(i + (current_layer_add*elem_per_layer),5),node_ids6) || isequal(node_ids(i + (current_layer_add*elem_per_layer),6),node_ids6)) ...
                    && i + (current_layer_add*elem_per_layer) ~= j)
                neighb_ele(1) = 1;
                corresp_nodes(1,1) = node_ids2;
                corresp_nodes(2,1) = node_ids3;
                corresp_nodes(3,1) = node_ids5;
                corresp_nodes(4,1) = node_ids6;
                ind_neighb_neu(1) = i + (current_layer_add*elem_per_layer);
                
            end
            
            % Second possible lateral neighbour AC:
            if ((isequal(node_ids(i + (current_layer_add*elem_per_layer),1),node_ids1) || isequal(node_ids(i + (current_layer_add*elem_per_layer),2),node_ids1) || isequal(node_ids(i + (current_layer_add*elem_per_layer),3),node_ids1)) ...
                    && (isequal(node_ids(i + (current_layer_add*elem_per_layer),1),node_ids3) || isequal(node_ids(i + (current_layer_add*elem_per_layer),2),node_ids3) || isequal(node_ids(i + (current_layer_add*elem_per_layer),3),node_ids3))...
                    && (isequal(node_ids(i + (current_layer_add*elem_per_layer),4),node_ids4) || isequal(node_ids(i + (current_layer_add*elem_per_layer),5),node_ids4) || isequal(node_ids(i + (current_layer_add*elem_per_layer),6),node_ids4))...
                    && (isequal(node_ids(i + (current_layer_add*elem_per_layer),4),node_ids6) || isequal(node_ids(i + (current_layer_add*elem_per_layer),5),node_ids6) || isequal(node_ids(i + (current_layer_add*elem_per_layer),6),node_ids6))...
                    && i + (current_layer_add*elem_per_layer) ~= j)
                neighb_ele(2) = 2;
                corresp_nodes(1,2) = node_ids1;
                corresp_nodes(2,2) = node_ids3;
                corresp_nodes(3,2) = node_ids4;
                corresp_nodes(4,2) = node_ids6;
                ind_neighb_neu(2) = i + (current_layer_add*elem_per_layer);
            end
            
            % Third possible lateral neighbour AB:
            if ((isequal(node_ids(i + (current_layer_add*elem_per_layer),1),node_ids1) || isequal(node_ids(i + (current_layer_add*elem_per_layer),2),node_ids1) || isequal(node_ids(i + (current_layer_add*elem_per_layer),3),node_ids1)) ...
                    && (isequal(node_ids(i + (current_layer_add*elem_per_layer),1),node_ids2) || isequal(node_ids(i + (current_layer_add*elem_per_layer),2),node_ids2) || isequal(node_ids(i + (current_layer_add*elem_per_layer),3),node_ids2))...
                    && (isequal(node_ids(i + (current_layer_add*elem_per_layer),4),node_ids4) || isequal(node_ids(i + (current_layer_add*elem_per_layer),5),node_ids4) || isequal(node_ids(i + (current_layer_add*elem_per_layer),6),node_ids4))...
                    && (isequal(node_ids(i + (current_layer_add*elem_per_layer),4),node_ids5) || isequal(node_ids(i + (current_layer_add*elem_per_layer),5),node_ids5) || isequal(node_ids(i + (current_layer_add*elem_per_layer),6),node_ids5))...
                    && i + (current_layer_add*elem_per_layer) ~= j)
                neighb_ele(3) = 3;
                corresp_nodes(1,3) = node_ids1;
                corresp_nodes(2,3) = node_ids2;
                corresp_nodes(3,3) = node_ids4;
                corresp_nodes(4,3) = node_ids5;
                ind_neighb_neu(3) = i + (current_layer_add*elem_per_layer);
            end
            
            % Lower neighbour:
            if current_layer_add>0
                if ((isequal(node_ids(i + ((current_layer_add-1)*elem_per_layer),4),node_ids1) || isequal(node_ids(i + ((current_layer_add-1)*elem_per_layer),5),node_ids1) || isequal(node_ids(i + ((current_layer_add-1)*elem_per_layer),6),node_ids1)) ...
                        && (isequal(node_ids(i + ((current_layer_add-1)*elem_per_layer),4),node_ids2) || isequal(node_ids(i + ((current_layer_add-1)*elem_per_layer),5),node_ids2) || isequal(node_ids(i + ((current_layer_add-1)*elem_per_layer),6),node_ids2))...
                        && (isequal(node_ids(i + ((current_layer_add-1)*elem_per_layer),4),node_ids3) || isequal(node_ids(i + ((current_layer_add-1)*elem_per_layer),5),node_ids3) || isequal(node_ids(i + ((current_layer_add-1)*elem_per_layer),6),node_ids3))...
                        && i + ((current_layer_add-1)*elem_per_layer) ~= j)
                    neighb_ele(4) = 4;
                    corresp_nodes(1,4) = node_ids1;
                    corresp_nodes(2,4) = node_ids2;
                    corresp_nodes(3,4) = node_ids3;
                    ind_neighb_neu(4) = i + ((current_layer_add-1)*elem_per_layer);
                end
            end
            % Upper neighbour:
            if current_layer_add< no_layers-1
                if ((isequal(node_ids(i + ((current_layer_add+1)*elem_per_layer),1),node_ids4) || isequal(node_ids(i + ((current_layer_add+1)*elem_per_layer),2),node_ids4) || isequal(node_ids(i + ((current_layer_add+1)*elem_per_layer),3),node_ids4)) ...
                        && (isequal(node_ids(i + ((current_layer_add+1)*elem_per_layer),1),node_ids5) || isequal(node_ids(i + ((current_layer_add+1)*elem_per_layer),2),node_ids5) || isequal(node_ids(i + ((current_layer_add+1)*elem_per_layer),3),node_ids5))...
                        && (isequal(node_ids(i + ((current_layer_add+1)*elem_per_layer),1),node_ids6) || isequal(node_ids(i + ((current_layer_add+1)*elem_per_layer),2),node_ids6) || isequal(node_ids(i + ((current_layer_add+1)*elem_per_layer),3),node_ids6))...
                        && i + ((current_layer_add+1)*elem_per_layer) ~= j)
                    neighb_ele(5) = 5;
                    corresp_nodes(1,5) = node_ids4;
                    corresp_nodes(2,5) = node_ids5;
                    corresp_nodes(3,5) = node_ids6;
                    ind_neighb_neu(5) = i + ((current_layer_add+1)*elem_per_layer);
                end
            end
            
        end
        
    elseif j>elem_per_layer
        
        layer_beneath = j - elem_per_layer;
        neighbs_elem_beneath = all_neighb_indices(layer_beneath,:);
        
        
        if neighbs_elem_beneath(1)~=0 && neighbs_elem_beneath(1)<no_elem
            neighb_ele(1) = 1;
            corresp_nodes(1,1) = node_ids2;
            corresp_nodes(2,1) = node_ids3;
            corresp_nodes(3,1) = node_ids5;
            corresp_nodes(4,1) = node_ids6;
            ind_neighb_neu(1) = neighbs_elem_beneath(1) + elem_per_layer;
        end
        
        if neighbs_elem_beneath(2)~=0 && neighbs_elem_beneath(2)<no_elem
            neighb_ele(2) = 2;
            corresp_nodes(1,2) = node_ids1;
            corresp_nodes(2,2) = node_ids3;
            corresp_nodes(3,2) = node_ids4;
            corresp_nodes(4,2) = node_ids6;
            ind_neighb_neu(2) = neighbs_elem_beneath(2) + elem_per_layer;
        end
        
        if neighbs_elem_beneath(3)~=0 && neighbs_elem_beneath(3)<no_elem
            neighb_ele(3) = 3;
            corresp_nodes(1,3) = node_ids1;
            corresp_nodes(2,3) = node_ids2;
            corresp_nodes(3,3) = node_ids4;
            corresp_nodes(4,3) = node_ids5;
            ind_neighb_neu(3) = neighbs_elem_beneath(3) + elem_per_layer;
        end
        
        if current_layer_add>0
            neighb_ele(4) = 4;
            corresp_nodes(1,4) = node_ids1;
            corresp_nodes(2,4) = node_ids2;
            corresp_nodes(3,4) = node_ids3;
            ind_neighb_neu(4) = j - elem_per_layer;
        end
        
        if current_layer_add< no_layers-1
            neighb_ele(5) = 5;
            corresp_nodes(1,5) = node_ids4;
            corresp_nodes(2,5) = node_ids5;
            corresp_nodes(3,5) = node_ids6;
            ind_neighb_neu(5) = j + elem_per_layer;
        end
        
    end
    
    % ===================================================================
    % Delete empty columns in case the number of neighbourss is smaller than
    % the possible maximum number of neighbours (i.e. 3) and determine the actual
    % elemental ids of the neighbouring elements
    
    corresp_nodes_prev = corresp_nodes;
    check_sum = sum(corresp_nodes_prev);
    
    corresp_nodes( :, all(~corresp_nodes,1)) = [];
    ind_neighb = [find(neighb_ele==1)' find(neighb_ele==2)' find(neighb_ele==3)' find(neighb_ele==4)' find(neighb_ele==5)'];
    ind_neighb = ind_neighb_neu(ind_neighb);
    
    ind_neighb_lat = [find(neighb_ele==1)' find(neighb_ele==2)' find(neighb_ele==3)'];
    ind_neighb_lat = ind_neighb_neu(ind_neighb_lat);
    
    no_neighb = length(ind_neighb); % number of neighbours
    no_neighb_lat = length(ind_neighb_lat); % number of lateral_neighb
    
    l_quest = find(neighb_ele==4);
    u_quest = find(neighb_ele==5);
    
    if current_layer_add>0
        ind_lower_neighb = j - elem_per_layer;
    else
        ind_lower_neighb = [];
    end
    
    if current_layer_add<no_layers-1
        ind_upper_neighb = j + elem_per_layer;
    else
        ind_upper_neighb = [];
    end
    
    lower_neighb = ~isempty(l_quest);
    upper_neighb = ~isempty(u_quest);
    
    
    %%
    % ===================================================================
    % Determination of the lengths of the edges, the intersection points
    % of the lines between the centroids and the edges Also, determination
    % of the lengths between the centroid of an element i to the
    % edge between i and neighbour i-1 (i.e. l_i) and the corresponding
    % length between the centroid of element i-1 and the same edge (i.e.
    % l_i-1). Also, compuation of first total flux estimates (not orthogonal
    % to the edges) :
    % ===================================================================
    % Preallcoations
    L_f = zeros(no_neighb_lat,1); % lengths of the edges
    centros_neighbs = zeros(2,no_neighb_lat);
    xy_coord_L_f_nodes = zeros(2,2,no_neighb_lat);
    
    P2 = zeros(2,no_neighb_lat); % centroids of neighbouring elements
    P3_4 = zeros(2,2,no_neighb_lat); % end points of the edges
    
    % coordinates of the centroid of the current element i:
    x_centro_of_i = centro_of_elem(j,1);
    y_centro_of_i = centro_of_elem(j,2);
    z_centro_of_i = centro_of_elem(j,3);
    
    % initializations
    l_total = zeros(1,no_neighb_lat); % lengths between centroid of element i and edges
    points_of_intersection = zeros(2,no_neighb_lat); % points of intersection of the line between
    % the centroids and the shared edge
    l_elem = zeros(1,no_neighb); % part of l_total lying in i
    l_neighb = zeros(1,no_neighb); % part of l_total lying in i-1
    
    P1 = [x_centro_of_i y_centro_of_i]; % centroid of triangular base of element i
    
    for i = 1:no_neighb_lat
        
        node_1 = corresp_nodes(1,i);
        node_2 = corresp_nodes(2,i);
        
        x_coord_node1 = coordinates(node_1,1);
        y_coord_node1 = coordinates(node_1,2);
        x_coord_node2 = coordinates(node_2,1);
        y_coord_node2 = coordinates(node_2,2);
        
        xy_coord_L_f_nodes(1,1,i) = x_coord_node1;
        xy_coord_L_f_nodes(2,1,i) = y_coord_node1;
        xy_coord_L_f_nodes(1,2,i) = x_coord_node2;
        xy_coord_L_f_nodes(2,2,i) = y_coord_node2;
        
        L_f(i) = sqrt(((x_coord_node2 - x_coord_node1)^2) + ((y_coord_node2 - y_coord_node1)^2));
        
        % Determination of the actual centroids of the neighbouring elements
        centros_neighbs(1,i)=centro_of_elem(ind_neighb(i),1);
        centros_neighbs(2,i)=centro_of_elem(ind_neighb(i),2);
        
        P2 = centros_neighbs; % centroids of neighbouring elements
        P3_4 = xy_coord_L_f_nodes; % end points of the edges
        
        l_total(i) = sqrt((x_centro_of_i - centros_neighbs(1,i))^2 + (y_centro_of_i - centros_neighbs(2,i))^2);
        points_of_intersection(1,i)=(((P3_4(1,2,i) - P3_4(1,1,i)) * ((P2(1,i)*P1(2))-(P1(1)*P2(2,i)))) - ((P2(1,i)-P1(1))*((P3_4(1,2,i)*P3_4(2,1,i))-(P3_4(1,1,i)*P3_4(2,2,i)))))/ ...
            (((P3_4(2,2,i)-P3_4(2,1,i))*(P2(1,i)-P1(1)))-((P2(2,i)-P1(2))*(P3_4(1,2,i)-P3_4(1,1,i))));%x values of intersection point
        points_of_intersection(2,i)=(((P1(2)-P2(2,i))*((P3_4(1,2,i)*P3_4(2,1,i))-(P3_4(1,1,i)*P3_4(2,2,i))))-((P3_4(2,1,i)-P3_4(2,2,i))*((P2(1,i)*P1(2))-(P1(1)*P2(2,i)))))/ ...
            (((P3_4(2,2,i)-P3_4(2,1,i))*(P2(1,i)-P1(1)))-((P2(2,i)-P1(2))*(P3_4(1,2,i)-P3_4(1,1,i))));%y values of intersection point
        
        l_elem(i)= sqrt((P1(1)-points_of_intersection(1,i))^2 + (P1(2)-points_of_intersection(2,i))^2); % with P1 and points_of_intersection
        l_neighb(i)= sqrt((centros_neighbs(1,i) - points_of_intersection(1,i))^2 + (centros_neighbs(2,i)-points_of_intersection(2,i))^2);
        
    end
    
    if no_neighb_lat==0
        i = 0;
    end
    
    
    %% top down lengths
    
    if lower_neighb
        i = i+1;
        z_bound = (zA + zB + zC)/3;
        z_centro_elem = centro_of_elem(j,3);
        z_centro_neighb = centro_of_elem(ind_lower_neighb,3);
        l_elem(i) = z_centro_elem - z_bound;
        l_neighb(i) = z_bound - z_centro_neighb;
    end
    
    if upper_neighb
        i = i+1;
        z_bound = (zD + zE + zF)/3;
        z_centro_elem = centro_of_elem(j,3);
        z_centro_neighb = centro_of_elem(ind_upper_neighb,3);
        l_elem(i) =z_bound - z_centro_elem;
        l_neighb(i) = z_centro_neighb - z_bound;
    end
    
    
    
    %% ALl the base areas
    
    A_for_q_to_Q = zeros(1,no_neighb);
    h_for_q_to_Q = zeros(1,no_neighb_lat);
    
    
    %% A horizontal
    a_horiz = sqrt((xA - xB)^2 + (yA - yB)^2);
    b_horiz = sqrt((xA - xC)^2 + (yA - yC)^2);
    c_horiz = sqrt((xB - xC)^2 + (yB - yC)^2);
    
    s_horiz = (a_horiz+b_horiz+c_horiz)/2;
    A_horiz(j) = sqrt(s_horiz*(s_horiz-a_horiz)*(s_horiz-b_horiz)*(s_horiz-c_horiz));
    
    Vol_prisms(j) = A_horiz(j)*(abs((zD+zE+zF)/3) - ((zA+zB+zC)/3)); % Volumina prisms
    %% A_top
    a_top = sqrt((xD - xE)^2 + (yD - yE)^2 + (zD - zE)^2);
    b_top = sqrt((xD - xF)^2 + (yD - yF)^2 + (zD - zF)^2);
    c_top = sqrt((xE - xF)^2 + (yE - yF)^2 + (zE - zF)^2);
    
    s_top = (a_top+b_top+c_top)/2;
    A_top(j) = sqrt(s_top*(s_top-a_top)*(s_top-b_top)*(s_top-c_top));
    
    %%
    
    h_counter = 1;
    for i = 1:no_neighb
        
        node_1 = corresp_nodes(1,i);
        node_2 = corresp_nodes(2,i);
        node_3 = corresp_nodes(3,i);
        node_4 = corresp_nodes(4,i);
        
        x_coord_node1 = coordinates(node_1,1);
        y_coord_node1 = coordinates(node_1,2);
        z_coord_node1 = coordinates(node_1,3);
        x_coord_node2 = coordinates(node_2,1);
        y_coord_node2 = coordinates(node_2,2);
        z_coord_node2 = coordinates(node_2,3);
        x_coord_node3 = coordinates(node_3,1);
        y_coord_node3 = coordinates(node_3,2);
        z_coord_node3 = coordinates(node_3,3);
        
        
        if node_4==0
            a_Q = sqrt((x_coord_node1 - x_coord_node2)^2 + (y_coord_node1 - y_coord_node2)^2 + (z_coord_node1 - z_coord_node2)^2);
            b_Q = sqrt((x_coord_node1 - x_coord_node3)^2 + (y_coord_node1 - y_coord_node3)^2 + (z_coord_node1 - z_coord_node3)^2);
            c_Q = sqrt((x_coord_node2 - x_coord_node3)^2 + (y_coord_node2 - y_coord_node3)^2 + (z_coord_node2 - z_coord_node3)^2);
            
            s_Q = (a_Q + b_Q + c_Q)/2;
            A_for_q_to_Q(i) = sqrt(s_Q*(s_Q-a_Q)*(s_Q-b_Q)*(s_Q-c_Q)); % area of lower face
        else
            x_coord_node4 = coordinates(node_4,1);
            y_coord_node4 = coordinates(node_4,2);
            z_coord_node4 = coordinates(node_4,3);
            LQ = sqrt((x_coord_node1 - x_coord_node2)^2 + (y_coord_node1 - y_coord_node2)^2);
            hQ_d = (z_coord_node1 + z_coord_node2)/2;
            hQ_u = (z_coord_node3 + z_coord_node4)/2;
            hQ = sqrt((hQ_d - hQ_u)^2);
            A_for_q_to_Q(i) = LQ*hQ;
            
            h_for_q_to_Q(h_counter) = hQ;
            h_counter = h_counter + 1;
            
        end
    end
    
    
    % ===================================================================
    % First check if there is a flux which equals zero
    % and create a new vector which potentially include these fluxes
    % and which accounts for the boundaries as no flow boundaries as no
    % flux can be calculated when there is no neighb element.
    
    % introduce zeros again as placeholder, if no neighbour is there at the
    % repective place
    
    if no_neighb_lat<3
        missing_ind = find(check_sum(1:3)==0);
        existing_ind = find(check_sum(1:3)~=0);
        h_f_neu = zeros(1,3);
        
        h_f_neu(missing_ind) = 0;
        h_f_neu(existing_ind) = h_for_q_to_Q;
    else
        h_f_neu = h_for_q_to_Q;
        
    end
    
    
    is_at_bound = find(ind_neighb_neu==0);
    
    if no_neighb<no_poss_neighbs
        
        missing_ind = find(ind_neighb_neu==0);
        existing_ind = find(ind_neighb_neu~=0);
        
        A_f_neu = zeros(1,no_poss_neighbs);
        l_elem_neu = zeros(1,no_poss_neighbs);
        l_neighb_neu = zeros(1,no_poss_neighbs);
        
        A_f_neu(missing_ind) = 0;
        A_f_neu(existing_ind) = A_for_q_to_Q;
        
        l_elem_neu(missing_ind) = 0;
        l_elem_neu(existing_ind) = l_elem;
        
        l_neighb_neu(missing_ind) = 0;
        l_neighb_neu(existing_ind) = l_neighb;
        
    else
        A_f_neu = A_for_q_to_Q;
        l_elem_neu = l_elem;
        l_neighb_neu = l_neighb;
    end
    
    %% orthogonal projection.
    
    ind_Q_not_zero = find(corresp_nodes_prev(1,:)~=0);
    
    corr_fact_orth = NaN(no_neighb,1);
    
    
    for i = 1:no_neighb
        
        ind_neighb_act = ind_neighb(i);
        vec_between_centroids =  centro_of_elem(ind_neighb_act,:) - centro_of_elem(j,:);
        abs_vec_between_centroids = sqrt(vec_between_centroids(1)^2 + vec_between_centroids(2)^2 + vec_between_centroids(3)^2);
        
        node_1 = corresp_nodes(1,i);
        node_2 = corresp_nodes(2,i);
        node_3 = corresp_nodes(3,i);
        node_4 = corresp_nodes(4,i);
        
        x_coord_node1 = coordinates(node_1,1);
        y_coord_node1 = coordinates(node_1,2);
        z_coord_node1 = coordinates(node_1,3);
        x_coord_node2 = coordinates(node_2,1);
        y_coord_node2 = coordinates(node_2,2);
        z_coord_node2 = coordinates(node_2,3);
        x_coord_node3 = coordinates(node_3,1);
        y_coord_node3 = coordinates(node_3,2);
        z_coord_node3 = coordinates(node_3,3);
        
        
        A_dir1 = [x_coord_node1-x_coord_node2;y_coord_node1-y_coord_node2;z_coord_node1-z_coord_node2]; % one directional vector on corresponding face
        A_dir2 = [x_coord_node3-x_coord_node2;y_coord_node3-y_coord_node2;z_coord_node3-z_coord_node2]; % another directional vector on corresp. face
        Face_vec_normal = [(A_dir1(2)*A_dir2(3))-(A_dir1(3)*A_dir2(2));(A_dir1(3)*A_dir2(1))-(A_dir1(1)*A_dir2(3));...
            (A_dir1(1)*A_dir2(2)) - (A_dir1(2)*A_dir2(1))]; % a normal vector on face
        
        orth_part_of_vec_centro = ((vec_between_centroids(1)*Face_vec_normal(1) + vec_between_centroids(2)*Face_vec_normal(2) + vec_between_centroids(3)*Face_vec_normal(3))...
            /((Face_vec_normal(1)^2) + (Face_vec_normal(2))^2 + (Face_vec_normal(3))^2 ))*Face_vec_normal;
        abs_orth_part_of_vec_centro = sqrt((orth_part_of_vec_centro(1)^2) + (orth_part_of_vec_centro(2)^2) + (orth_part_of_vec_centro(3))^2);
        
        corr_fact_orth(i) = abs_orth_part_of_vec_centro/abs_vec_between_centroids;
        
        Face_vec_normal_abs = abs(Face_vec_normal);
        
    end
    
    corr_fact_orth_neu = zeros(1,no_poss_neighbs);
    
    if no_neighb<no_poss_neighbs
        corr_fact_orth_neu(existing_ind) = corr_fact_orth;
    else
        corr_fact_orth_neu = corr_fact_orth;
    end
    
    %% Insert IDs and Values for Nodes on Dirichlet-boundary:
    
    if flag_at_DiriNeumann
        
        no_DiriNeum_faces_on_ele = face_index_mat(j,:);
        no_DiriNeum_faces_on_ele(no_DiriNeum_faces_on_ele==0) = [];
        
        for ff = 1:length(no_DiriNeum_faces_on_ele)
            
            ind_neighb_neu(face_index_mat(j,ff)) = no_elem + diri_counter;
            
            %%
            if face_index_mat(j,ff)==1 % 1st neighbour: 2,3,5,6
                % B,C,E,F
                xDiri1 = xB; yDiri1 = yB; zDiri1 = zB;
                xDiri2 = xC; yDiri2 = yC; zDiri2 = zC;
                xDiri3 = xE; yDiri3 = yE; zDiri3 = zE;
                xDiri4 = xF; yDiri4 = yF; zDiri4 = zF;
            elseif face_index_mat(j,ff)==2 % 2nd neighbour: 1,3,4,6
                % A,C,D,F
                xDiri1 = xA; yDiri1 = yA; zDiri1 = zA;
                xDiri2 = xC; yDiri2 = yC; zDiri2 = zC;
                xDiri3 = xD; yDiri3 = yD; zDiri3 = zD;
                xDiri4 = xF; yDiri4 = yF; zDiri4 = zF;
            elseif face_index_mat(j,ff)==3 % 3rd neighbour: 1,2,4,5
                % A,B,D,E
                xDiri1 = xA; yDiri1 = yA; zDiri1 = zA;
                xDiri2 = xB; yDiri2 = yB; zDiri2 = zB;
                xDiri3 = xD; yDiri3 = yD; zDiri3 = zD;
                xDiri4 = xE; yDiri4 = yE; zDiri4 = zE;
            elseif face_index_mat(j,ff)==4 % 4th neighbour: 1,2,3
                % A,B,C
                xDiri1 = xA; yDiri1 = yA; zDiri1 = zA;
                xDiri2 = xB; yDiri2 = yB; zDiri2 = zB;
                xDiri3 = xC; yDiri3 = yC; zDiri3 = zC;
            elseif face_index_mat(j,ff)==5 % 5th neighbour: 4,5,6
                % D,E,F
                xDiri1 = xD; yDiri1 = yD; zDiri1 = zD;
                xDiri2 = xE; yDiri2 = yE; zDiri2 = zE;
                xDiri3 = xF; yDiri3 = yF; zDiri3 = zF;
            end
            
            if face_index_mat(j,ff)==1 || face_index_mat(j,ff)==2 || face_index_mat(j,ff)==3
                LQ = sqrt((xDiri1 - xDiri2)^2 + (yDiri1 - yDiri2)^2);
                hQ_d = (zDiri1 + zDiri2)/2;
                hQ_u = (zDiri3 + zDiri4)/2;
                hQ = sqrt((hQ_d - hQ_u)^2);
                
                A_f_neu(face_index_mat(j,ff)) = LQ*hQ;
                h_f_neu(face_index_mat(j,ff)) = hQ;
                
            elseif face_index_mat(j,ff)==4 || face_index_mat(j,ff)==5
                a_Q = sqrt((xDiri1 - xDiri2)^2 + (yDiri1 - yDiri2)^2 + (zDiri1 - zDiri2)^2);
                b_Q = sqrt((xDiri1 - xDiri3)^2 + (yDiri1 - yDiri3)^2 + (zDiri1 - zDiri3)^2);
                c_Q = sqrt((xDiri2 - xDiri3)^2 + (yDiri2 - yDiri3)^2 + (zDiri2 - zDiri3)^2);
                
                s_Q = (a_Q + b_Q + c_Q)/2;
                
                A_f_neu(face_index_mat(j,ff)) = sqrt(s_Q*(s_Q-a_Q)*(s_Q-b_Q)*(s_Q-c_Q)); % area of lower face
            end
            
            % lengths
            p1 = [xDiri1 - xDiri3; yDiri1 - yDiri3; zDiri1 - zDiri3];
            p2 = [xDiri2 - xDiri3; yDiri2 - yDiri3; zDiri2 - zDiri3];
            
            nE = [(p1(2)*p2(3)) - (p1(3)*p2(2)); (p1(3)*p2(1)) - (p1(1)*p2(3)); ...
                (p1(1)*p2(2)) - (p1(2)*p2(1))];
            nEnorm = nE/(sqrt(nE(1)^2 + nE(2)^2 +nE(3)^2));
            
            if face_index_mat(j,ff)==4 || face_index_mat(j,ff)==5
                r0 = [(xDiri1+xDiri2+xDiri3)/3;(yDiri1+yDiri2+yDiri3)/3;(zDiri1+zDiri2+zDiri3)/3];
            elseif face_index_mat(j,ff)==1 || face_index_mat(j,ff)==2 || face_index_mat(j,ff)==3
                r0 = [(xDiri1+xDiri2+xDiri3+xDiri4)/4;(yDiri1+yDiri2+yDiri3+yDiri4)/4;(zDiri1+zDiri2+zDiri3+zDiri4)/4];
            end
            dist = [x_centro_of_i-r0(1);y_centro_of_i-r0(2);z_centro_of_i-r0(3)];
            
            
            P_v = dist(1)*nEnorm(1) + dist(2)*nEnorm(2) + dist(3)*nEnorm(3);
            P_on_Diri = [x_centro_of_i;y_centro_of_i;z_centro_of_i] - (P_v*nEnorm);
            
            l_elem_neu(face_index_mat(j,ff)) = sqrt((x_centro_of_i - P_on_Diri(1))^2 + (y_centro_of_i - P_on_Diri(2))^2 + (z_centro_of_i - P_on_Diri(3))^2);
            
            nEnorm_abs = abs(nEnorm);
            diri_counter = diri_counter + 1;
        end
        
        
    end
    
    if A_f_neu(5)==0
        a_f = sqrt((xD - xE)^2 + (yD - yE)^2 + (zD - zE)^2);
        b_f = sqrt((xD - xF)^2 + (yD - yF)^2 + (zD - zF)^2);
        c_f = sqrt((xE - xF)^2 + (yE - yF)^2 + (zE - zF)^2);
        
        s_f = (a_f + b_f + c_f)/2;
        A_f_neu(5) = sqrt(s_f*(s_f-a_f)*(s_f-b_f)*(s_f-c_f));
        
    end
    
    all_neighb_indices(j,:) =ind_neighb_neu;
    A_f_all(j,:) = A_f_neu;
    l_elem_all(j,:) = l_elem_neu;
    l_neighb_all(j,:) = l_neighb_neu;
    corr_fact_orth_save(j,ind_Q_not_zero) = corr_fact_orth;
    hQ_all(j,:) = h_f_neu;
end

% Necessary only for the flux reconstruction

Vol_prisms_name = strcat('Vol_prisms_',prefix,'.mat');
save(Vol_prisms_name,'Vol_prisms')

corr_fact_orth_name = strcat('corr_fact_orth_save_',prefix,'.mat');
save(corr_fact_orth_name,'corr_fact_orth_save')

A_f_all_name = strcat('A_f_all_',prefix,'.mat');
save(A_f_all_name,'A_f_all')

l_elem_all_name = strcat('l_elem_all_',prefix,'.mat');
save(l_elem_all_name,'l_elem_all')

l_neighb_all_name = strcat('l_neighb_all_',prefix,'.mat');
save(l_neighb_all_name,'l_neighb_all')

centro_of_elem_name = strcat('centro_of_elem_',prefix,'.mat');
save(centro_of_elem_name,'centro_of_elem')

% Necessary for particle tracking and the flux reconstruction

A_horiz_name = strcat('A_horiz_',prefix,'.mat');
save(A_horiz_name,'A_horiz')

A_top_name = strcat('A_top_',prefix,'.mat');
save(A_top_name,'A_top')

all_neighb_indices_name = strcat('all_neighb_indices_',prefix,'.mat');
save(all_neighb_indices_name,'all_neighb_indices')

hQ_all_name = strcat('hQ_all_',prefix,'.mat');
save(hQ_all_name,'hQ_all')

face_index_mat_name = strcat('face_index_mat_',prefix,'.mat');
save(face_index_mat_name,'face_index_mat')

face_flag_D_N_name = strcat('face_flag_D_N_',prefix,'.mat');
save(face_flag_D_N_name,'face_flag_D_N')

Dirichlet_elements_mult_name = strcat('Dirichlet_elements_mult_',prefix,'.mat');
save(Dirichlet_elements_mult_name,'vector_of_d_bounds_mult')

face_index_mat_multD_name = strcat('face_index_mat_multD_',prefix,'.mat');
save(face_index_mat_multD_name,'face_index_d_mult')

Neumann_elements_mult_name = strcat('Neumann_elements_mult_',prefix,'.mat');
save(Neumann_elements_mult_name,'vector_of_n_bounds_mult')

face_index_mat_multN_name = strcat('face_index_mat_mult_',prefix,'.mat');
save(face_index_mat_multN_name,'face_index_n_mult')

disp('ready')
toc

