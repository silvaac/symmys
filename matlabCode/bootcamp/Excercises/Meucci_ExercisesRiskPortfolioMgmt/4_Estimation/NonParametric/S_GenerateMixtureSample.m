clear; clc; close all

a=.8;
m_Y=.1;
s_Y=.2;
m_Z=0;
s_Z=.15;

T=52;

P=rand(T,1);
Q=QuantileMixture(P,a,m_Y,s_Y,m_Z,s_Z);
