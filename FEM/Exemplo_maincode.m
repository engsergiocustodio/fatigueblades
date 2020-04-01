clear all; close all; clc;

journal_name='Journa1_0.wbjn';
geom_name='Geometry.js';
model_name='Model1_0.js';

shape_name='shape235.txt';
root_name='elipses235.txt';
perfil_name='NACA653618.txt';
npontos_perfil=50;
xlongarina=0.0;
ylongarina=0.0;
dlongarina=140;
perfil_fim_longarina=7;

load_name='1_cargas0_5292.txt';
vel_rotation=0.5292; %rad/s

TURBINE_ANSYS(journal_name,geom_name,model_name,shape_name,root_name,perfil_name,npontos_perfil,xlongarina,ylongarina,dlongarina,perfil_fim_longarina,load_name,vel_rotation)