% Transformer Yy0 test file - harmonic power flow

function mpc = xfmr_yy0
mpc.version = '2';
mpc.baseMVA =  100.00;

%% bus data
%	bus_id	type    Pd      Qd	    Gs	    Bs	    area	Vm	    Va	    baseKV  zone	Vmax	Vmin
mpc.bus = [
    0       3       0.000   0.000   0.00    0.00    1       1.00    0.00    150.0   1       1.06    0.94;
    1       1       0.000   0.000   0.00    0.00    1       1.00    0.00    150.0   1       1.06    0.94;
    2       1       10.00   5.00    0.00    0.00    1       1.00	0.00    35.4    1       1.06    0.94;
];

%% bus harmonic data 
%column_names%  nh_1    nh_3    nh_5    nh_7    nh_9    nh_13
mpc.bus_harmonics = [
                0.000   0.0000  0.0000  0.0000  0.0000  0.0000; % 0
                0.000   0.0000  0.0000  0.0000  0.0000  0.0000; % 1
                1.000   0.2000  0.1000  0.0500  0.0600  0.0200; % 2
]

%% generator data
%   bus     Pg      Qg      Qmax    Qmin    Vg      mBase   status  Pmax    Pmin
mpc.gen = [
    0       0.0     0.0     100.00  -100.00 1.05    100.0   1       100.00  0;
];
mpc.gencost = [
	2	    0.0	    0.0	    3	    0.0     1.0     0.0;
];

%% branch data
%               f_bus	t_bus	r	        x	        b	        rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [  0       1       0.002377982 0.047559635 0.00000     1000    1000    1000    1       0       1       -60     60; % Ssc=2.1e9, XRr=20
];

%% transformer data
%column_names%  f_bus   t_bus   xsc     rsh     r1      r2      vg      gnd1    gnd2    re1     xe1     re2     xe2     rateA
mpc.xfmr = [    1       2       0.13    12.90   0.00211 0.00211 'Yy0'   0       0       0.0     0.0     0.0     0.0     125.0; % Snom=125e6 VA, uk=0.17, Pk=527e3 W, P0=77.5e3 W, N=500, A=0.5 m², l=11.4 m
];