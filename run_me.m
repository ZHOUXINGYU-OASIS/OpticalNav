%% Main demo
clear all;
close all;
clc;
%% if you want to reproduce the results in Fig. 10
case_id = 1; % case I: 1; case II: 2
if_rand = 0; % if use randomly generated initial errors
if_load_LOS = 1;
[tk, err] = OD_LOS(case_id, if_rand, if_load_LOS);
%% if you want to reproduce the results in Fig. 11
% case_id = 1; % case I: 1; case II: 2
% if_rand = 0; % if use randomly generated initial errors
% if_load_LOS = 1;
% [tk, err] = OD_LOS_VBR(case_id, if_rand, if_load_LOS);
%% if you want to reproduce the MC results of VBR-1
% case_id = 1; % case I: 1; case II: 2
% [ErrData, TimeData, stdData, rmseData, time_cost, MDData] = MC_OD_LOS_VBR(case_id);
%% if you want to reproduce the MC results of VBR-5
% case_id = 1; % case I: 1; case II: 2
% [ErrData, TimeData, stdData, rmseData, time_cost, MDData] = MC_OD_LOS_VBR5(case_id);