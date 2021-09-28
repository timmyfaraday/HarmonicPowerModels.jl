% Case to test space based matlab matrix
% And other hard to parse cases
% also test data without a generator cost model

function mpc = case_xfmr_Yy0
mpc.version = '2';
mpc.baseMVA =  10.00;

%% bus data
%	bus_id	type    Pd      Qd	    Gs	    Bs	    area	Vm	    Va	    baseKV  zone	Vmax	Vmin
mpc.bus = [
    1       3       0.000   0.000   0.00    0.00    1       1.00    0.00    150.0   1       1.10    0.90;
    2       1       50.00   15.88   0.00    0.00    1       1.00	0.00    36.0    1       1.10    0.90;
 %   3       1       0.000   0.000   0.00    0.00    1       1.00    0.00    36.0    1       1.10    0.90;
 %   4       1       12.60   4.001   0.00    0.00    1       1.00	0.00    10.0    1       1.10    0.90;
%   5       1       0.000   0.000   0.00    0.00    1       1.00    0.00    10.0    1       1.10    0.90;
%    6       1       1.420   0.000   0.00    0.00    1       1.00	0.00    0.69    1       1.10    0.90;
%    7       1       0.000   0.000   0.00    0.00    1       1.00    0.00    0.69    1       1.10    0.90;
%    8       1       0.050   0.000   0.00    0.00    1       1.00	0.00    0.40    1       1.10    0.90;
];

% note that it is realistic to add additional fundamental power consumption to
% nodes 2 and 4, at 40% of the preceeding transformer power and cos(phi) = 0.95
% today harmonic filters at 690 V
% dont forget to add capacitor banks again


%% bus harmonic data 
%column_names%  nh_1    nh_3    thdmax 
mpc.bus_harmonics = [
                1.000   0.0000   0.1; % 1
                1.000   0.0000  0.1; % 2
 %               1.000   0.0000  0.1; % 3
 %               1.000   0.0000   0.1; % 4
 %               1.000   0.0000   0.1; % 5
 %               1.000   0.0000   0.1; % 6
%                1.000   0.0000   0.1; % 7
 %               1.000   0.2000   0.1; % 8
]

%% Six pulse converter -- see document andreas 
%% TL Lights - https://www.dranetz.com/wp-content/uploads/2014/02/harmonics-understanding-thefacts-part3.pdf

%% generator data
%   bus     Pg      Qg      Qmax    Qmin    Vg      mBase   status  Pmax    Pmin
mpc.gen = [
    1       0.0     0.0     100.00  -100.00 1.05    100.0   1       100.00   -100;
 %   7       0.0     0.0     100.00    -100.00   1.05    100.0   1         100.00  -100;
];

mpc.gencost = [
	2	 0.0	 0.0	 3	   0.110000	   5.000000	   0.000000;
%	2	 0.0	 0.0	 3	   0.100000	   4.000000	   0.000000;
];

%% branch data
%   f_bus	t_bus	r	    x	    b	    rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [ 
 %   2       3       0.01    0.01       0       100    100    100           1   0       1       -60        60; 
 %   4       5       0.01    0.01       0       100    100    100           1   0       1       -60        60;
 %   6       7       0.01    0.01       0       100    100    100           1   0       1       -60        60;
];

%% transformer data
%column_names%  f_bus   t_bus   xsc     r1      r2      vg      gnd1    gnd2    re1     xe1     re2     xe2     rate_a  
mpc.xfmr = [
                1       2       0.0     0.1     0.1     'Yd11'  0       0       0.0     0.0     0.0     0.0     125;
  %              3       4       0.0     0.1     0.1     'Yy0'   0       0       0.0     0.0     0.0     0.0     31.5;
  %              5       6       0.0     0.1     0.1     'Dy11'  0       1       0.0     0.0     0.0     0.0     2.5;
   %             7       8       0.0     0.1     0.1     'Dz0'   0       1       0.0     0.0     0.0     0.0     0.25;
]
