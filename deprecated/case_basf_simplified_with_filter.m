% Case to test space based matlab matrix
% And other hard to parse cases
% also test data without a generator cost model

function mpc = case_xfmr_Yy0
mpc.version = '2';
mpc.baseMVA =  100.00;

%% bus data
%	bus_id	type    Pd      Qd	    Gs	    Bs	    area	Vm	    Va	    baseKV  zone	Vmax	Vmin
mpc.bus = [
    1       3       0.000   0.000   0.00    0.00    1       1.00    0.00    150.0   1       1.06    0.94;
    2       1       50.00   15.88   0.00    0.00    1       1.00	0.00    36.0    1       1.06    0.94;
    3       1       0.000   0.000   0.00    0.00    1       1.00    0.00    36.0    1       1.06    0.94;
    4       1       12.60   4.001   0.00    6.00    1       1.00	0.00    10.0    1       1.06    0.94;
    5       1       0.000   0.000   0.00    0.00    1       1.00    0.00    10.0    1       1.06    0.94;
    6       1       1.420   0.000   0.00    0.40    1       1.00	0.00    0.69    1       1.10    0.90;
    7       1       0.000   0.000   0.00    0.00    1       1.00    0.00    0.69    1       1.10    0.90;
    8       1       0.050   0.040   0.00    0.00    1       1.00	0.00    0.40    1       1.10    0.90;
];

% note that it is realistic to add additional fundamental power consumption to
% nodes 2 and 4, at 40% of the preceeding transformer power and cos(phi) = 0.95
% today harmonic filters at 690 V
% lighting -> cos phi of 0.6

%% bus harmonic data 
%column_names%  nh_1    nh_3     nh_5       nh_7    nh_9    nh_13     thdmax 
mpc.bus_harmonics = [
                1.000   0.0000  0.0000   0.0000  0.0000    0.0000  0.012; % 1
                1.000   0.0000  0.0000   0.0000  0.0000    0.0000  0.08; % 2
                1.000   0.0000  0.0000   0.0000  0.0000    0.0000 0.08; % 3
                1.000   0.0000  0.0000   0.0000  0.0000    0.0000 0.08; % 4
                1.000   0.0000  0.0000   0.0000  0.0000    0.0000 0.08; % 5
                1.000   0.0000  0.4500   0.2500  0.0000    0.0800 0.08; % 6
                1.000   0.0000  0.0000   0.0000  0.0000    0.0000 0.08; % 7
                1.000   0.2000  0.1000   0.0500  0.0600    0.0000 0.08; % 8
]

%% Six pulse converter -- see document andreas 
%% TL Lights - https://www.dranetz.com/wp-content/uploads/2014/02/harmonics-understanding-thefacts-part3.pdf

%% generator data
%   bus     Pg      Qg      Qmax    Qmin    Vg      mBase   status  Pmax    Pmin
mpc.gen = [
    1       0.0     0.0     100.00  -100.00 1.05    100.0   1       100.00   0;
    6       0.0     0.0     0.717    -0.717   1.05    100.0   1         0.717  -0.717;
];

mpc.gencost = [
	2	 0.0	 0.0	 3	   0.00000	   1.000000	   0.000000;
	2	 0.0	 0.0	 3	   0.0000	   0.00000	   0.000000;
];

%column_names%  isfilter
mpc.gen_extra = [
	0;
	1;
];

%% branch data
%   f_bus	t_bus	r	    x	    b	    rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [ 
    2       3       0.00060    0.00053       0.02513       125    100    100           1   0       1       -60        60; 
    4       5       0.00013    0.000085      0.002780      31.5    100    100           1   0       1       -60        60;
    6       7       0.000045 0.0000298       0.000973       2.5    100    100           1   0       1       -60        60;
];

%% transformer data
%column_names%  f_bus   t_bus   xsc     r1      r2      vg      gnd1    gnd2    re1     xe1     re2     xe2     rate_a  rsh
mpc.xfmr = [
                1       2       0.13     0.00211      0.00211     'Yd11'  0       0       0.0     0.0     0.0     0.0     125 1250; 
                3       4       0.0229    0.0003      0.0003     'Yy0'   0       0       0.0     0.0     0.0     0.0     31.5 6000;
                5       6       0.0107     0.0007     0.0007     'Dy11'  0       1       0.0     0.0     0.0000046     0.0000063     2.5 100000; % add PEN to re2 + jxe2
                7       8       0.0002     0.00005     0.00005     'Dz0'   0       1       0.0     0.0    0.000025     0.000013     0.25 100000; % ]
