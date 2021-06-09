% Case to test space based matlab matrix
% And other hard to parse cases
% also test data without a generator cost model

function mpc = case_xfmr_YNyn0
mpc.version = '2';
mpc.baseMVA =  100.00;

%% bus data
%	bus_id	type    Pd      Qd	    Gs	    Bs	    area	Vm	    Va	    baseKV  zone	Vmax	Vmin
mpc.bus = [
    1       3       0.00    0.00    0.00    0.00    1       1.00    0.00    150.00  1       1.10    0.90;
    2       1       100.00	50.00   0.00    0.00    1       1.00	0.00    35.00   1       1.10    0.90;
];

%% bus harmonic data 
%column_names%  nh_1    nh_5     nh_7    nh_11   nh_13
mpc.bus_harmonics = [
                1.000   0.000   0.000   0.000   0.000
                1.000   0.180   0.110   0.050   0.030
]

%% generator data
%   bus     Pg      Qg      Qmax    Qmin    Vg      mBase   status  Pmax    Pmin
mpc.gen = [
    1       100.0   50.0    200.00  -200.00 1.05    100.0   1       400.00  0.0;
];

%% branch data
%   fbus	tbus	r	    x	    b	    rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [ 
];

%% transformer data
%column_names%  fbus    tbus    xsc     r1      r2      vg      gnd1    gnd2    re1     xe1     re2     xe2     
mpc.xfmr = [
                1       2       0.001   0.12336 0.01644 'Yy0'   1       1       0.0     0.0     0.0     0.0
]
