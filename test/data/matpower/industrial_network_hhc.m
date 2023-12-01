% Industrial Network, see: `Harmonic Optimal Power Flow with Transformer
% Excitation` by F. Geth and T. Van Acker, pg. 8, § IV.B.

function mpc = industrial_network
mpc.version = '2';
mpc.baseMVA =  100.00;

%% bus data
%	bus_id	type    Pd      Qd	    Gs	    Bs	    area	Vm	    Va	    baseKV  zone	Vmax	Vmin
mpc.bus = [
    1       3       0.000   0.000   0.00    0.00    1       1.00    0.00    150.0   1       1.06    0.94;
    2       1       50.00   15.88   0.00    0.00    1       1.00	0.00    35.4    1       1.06    0.94;
    3       1       0.000   0.000   0.00    0.00    1       1.00    0.00    35.4    1       1.06    0.94;
    4       1       12.60   4.001   0.00    6.00    1       1.00	0.00    10.0    1       1.06    0.94;
    5       1       0.000   0.000   0.00    0.00    1       1.00    0.00    10.0    1       1.06    0.94;
    6       1       1.420   0.000   0.00    0.40    1       1.00	0.00    0.69    1       1.10    0.90;
    7       1       0.000   0.000   0.00    0.00    1       1.00    0.00    0.69    1       1.10    0.90;
    8       1       0.050   0.040   0.00    0.00    1       1.00	0.00    0.40    1       1.10    0.90;
];

%% bus harmonic data 
%column_names%  standard 
mpc.bus_harmonics = [
                'IEC61000-3-6:2008';        % 1
                'IEC61000-2-4:2002, Cl. 2'; % 2
                'IEC61000-2-4:2002, Cl. 2'; % 3
                'IEC61000-2-4:2002, Cl. 2'; % 4
                'IEC61000-2-4:2002, Cl. 2'; % 5
                'IEC61000-2-4:2002, Cl. 2'; % 6
                'IEC61000-2-4:2002, Cl. 2'; % 7
                'IEC61000-2-4:2002, Cl. 2'; % 8
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
    2       3       0.00060     0.00053     0.02513     125     1000    1000    1       0       1       -60     60; 
    4       5       0.00013     0.000085    0.002780    31.5    100     100     1       0       1       -60     60;
    6       7       0.000045    0.0000298   0.000973    2.5     100     100     1       0       1       -60     60;
];

%% transformer data
%column_names%  f_bus   t_bus   xsc     rsh     r1      r2      vg      gnd1    gnd2    re1     xe1     re2         xe2         rate_a
mpc.xfmr = [
                1       2       0.13    12.90   0.00211 0.00211 'Yd11'  0       0       0.0     0.0     0.0         0.0         125.0; 
                3       4       0.0229  70.42   0.0003  0.0003  'Yy0'   0       0       0.0     0.0     0.0         0.0         31.5;
                5       6       0.0107  370.37  0.0007  0.0007  'Dy11'  0       1       0.0     0.0     0.0000046   0.0000063   2.5;
                7       8       0.0002  1000.0  0.00005 0.00005 'Dz0'   0       1       0.0     0.0     0.000025    0.000013    0.25; 
];