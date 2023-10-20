% Industrial Network, see: `Harmonic Optimal Power Flow with Transformer
% Excitation` by F. Geth and T. Van Acker, pg. 8, § IV.B.

function mpc = industrial_network
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

%% bus harmonic data 
%column_names%  thdmax 
mpc.bus_harmonics = [
                0.016; % 1
                0.080; % 2
                0.080; % 3
                0.080; % 4
                0.080; % 5
                0.080; % 6
                0.080; % 7
                0.080; % 8
]

%% generator data
%   bus     Pg      Qg      Qmax    Qmin    Vg      mBase   status  Pmax    Pmin
mpc.gen = [
    1       0.0     0.0     100.00  -100.00 1.05    100.0   1       100.00  0;
];
mpc.gencost = [
	2	    0.0	    0.0	    3	    0.0     1.0     0.0;
];
%column_names%  isfilter
mpc.gen_extra = [
	0;
];

%% branch data
%   f_bus	t_bus	r	        x	        b	        rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [ 
    2       3       0.0         0.00053     0.02513     125     1000    1000    1       0       1       -60     60; 
    4       5       0.0         0.000085    0.002780    31.5    100     100     1       0       1       -60     60;
    6       7       0.0         0.0000298   0.000973    2.5     100     100     1       0       1       -60     60;
];

%% transformer data
%column_names%  f_bus   t_bus   xsc     rsh     r1      r2      vg      gnd1    gnd2    re1     xe1     re2         xe2         rate_a
mpc.xfmr = [
                1       2       0.13    99999.0 0.0     0.0     'Yd11'  0       0       0.0     0.0     0.0         0.0         125; 
                3       4       0.0229  99999.0 0.0     0.0     'Yy0'   0       0       0.0     0.0     0.0         0.0         31.5;
                5       6       0.0107  99999.0 0.0     0.0     'Dy11'  0       1       0.0     0.0     0.0         0.0000063   2.5;
                7       8       0.0002  99999.0 0.0     0.0     'Dz0'   0       1       0.0     0.0     0.0         0.000013    0.25; 
];