% Case to test space based matlab matrix
% And other hard to parse cases
% also test data without a generator cost model

function mpc = case_xfmr_Yy0
mpc.version = '2';
mpc.baseMVA =  100.00;

%% bus data
%	bus_id	type    Pd      Qd	    Gs	    Bs	    area	Vm	    Va	    baseKV  zone	Vmax	Vmin
mpc.bus = [
    1       3       0.000   0.000   0.00    0.00    1       1.00    0.00    12.47   1       1.10    0.90;
    2       1       5.000   2.000   0.00    0.00    1       1.00	0.00    4.16    1       1.10    0.90;
];

%% bus harmonic data 
%column_names%  nh_1    nh_3    nh_5    nh_7    nh_9    nh_11   nh_13   nh_15 thdmax
mpc.bus_harmonics = [
                1.000   0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000 0.05
                1.000   0.2  0.3  0.4180  0.0148  0.0708  0.0312  0.0048 0.05
]

%% generator data
%   bus     Pg      Qg      Qmax    Qmin    Vg      mBase   status  Pmax    Pmin
mpc.gen = [
    1       0.0     0.0     200.00  -200.00 1.05    100.0   1       400.00  0.0;
];

mpc.gencost = [
	2	 0.0	 0.0	 3	   0.110000	   5.000000	   0.000000;
];

%% branch data
%   f_bus	t_bus	r	    x	    b	    rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [ 
    1       2     0.01    0.01       0       100    100    100           1   0       1       -60        60;
];

