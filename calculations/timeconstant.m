clear all; clc; close all;

FILE = './time_constant_data.xlsx';

data = xlsread(FILE);
time = data(:, 1);
temp = data(:, 2);

plot(time, temp);
