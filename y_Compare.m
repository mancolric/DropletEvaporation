clear;clc;close all;

load("yv_np1.mat")
y_Ander = yv_np1;
load("y_fixed.mat")
y_Fixed = y;

matriz = [y_Ander,y_Fixed];