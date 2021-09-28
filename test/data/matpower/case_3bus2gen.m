% Case to test space based matlab matrix
% And other hard to parse cases
% also test data without a generator cost model

function mpc = case_xfmr_Yy0
mpc.version = '2';
mpc.baseMVA =  1.00;

%% bus data
%	bus_id	type    Pd      Qd	    Gs	    Bs	    area	Vm	    Va	    baseKV  zone	Vmax	Vmin
mpc.bus = [
    1       3       0.000   0.000   0.00    0.00    1       1.00    0.00    10    1      1.10    0.90;
    2       1       5.000   2.000   0.00    0.00    1       1.00	0.00    10    1      1.10    0.90;
    3       1       0.000   0.000   0.00    0.00    1       1.00	0.00    10    1      1.10    0.90;
];

%% bus harmonic data 
%column_names%  nh_1    nh_3     thdmax
mpc.bus_harmonics = [
                1.000   0.0000      0.5
                1.000   0.1         0.5
                1.000   0.0000      0.5
]

%% generator data
%   bus     Pg      Qg      Qmax    Qmin    Vg      mBase   status  Pmax    Pmin
mpc.gen = [
    1       0.0     0.0     9.00  -9.00 1.05    100.0   1       10.00  -10.0;
    3       0.0     0.0     9.00  -9.00 1.05    100.0   1       10.00  -10.0;
];

mpc.gencost = [
	2	 0.0	 0.0	 3	   0.110000	   5.000000	   0.000000;
	2	 0.0	 0.0	 3	   0.100000	   4.000000	   0.000000;
];

%% branch data
%   f_bus	t_bus	r	    x	    b	    rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [ 
    1       2     0.004    0.008       0       10    10    10           1   0       1       -60        60;
    2       3     0.004    0.008       0       4    4    4           1   0       1       -60        60;
];
