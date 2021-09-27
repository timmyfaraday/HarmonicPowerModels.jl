% Case to test space based matlab matrix
% And other hard to parse cases
% also test data without a generator cost model

function mpc = case_xfmr_Yy0
mpc.version = '2';
mpc.baseMVA =  100.00;

%% bus data
%	bus_id	type    Pd      Qd	    Gs	    Bs	    area	Vm	    Va	    baseKV  zone	Vmax	Vmin
mpc.bus = [
    1       3       0.000   0.000   0.00    0.00    1       1.00    0.00    150.0   1       1.10    0.90;
    2       1       0.000   0.000   0.00    0.00    1       1.00	0.00    36.0    1       1.10    0.90;
    3       1       0.000   0.000   0.00    0.00    1       1.00    0.00    36.0    1       1.10    0.90;
    4       1       0.000   0.000   0.00    0.00    1       1.00	0.00    10.0    1       1.10    0.90;
    5       1       0.000   0.000   0.00    0.00    1       1.00    0.00    10.0    1       1.10    0.90;
    6       1       1.420   0.000   0.00    0.00    1       1.00	0.00    0.69    1       1.10    0.90;
    7       1       0.000   0.000   0.00    0.00    1       1.00    0.00    0.69    1       1.10    0.90;
    8       1       0.050   0.000   0.00    0.00    1       1.00	0.00    0.40    1       1.10    0.90;
];

% note that it is realistic to add additional fundamental power consumption to
% nodes 2 and 4, at 40% of the preceeding transformer power and cos(phi) = 0.95

%% bus harmonic data 
%column_names%  nh_1    nh_3    nh_5    nh_7    nh_9    nh_13   nh_17   nh_19
mpc.bus_harmonics = [
                1.000   0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000; % 1
                1.000   0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000; % 2
                1.000   0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000; % 3
                1.000   0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000; % 4
                1.000   0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000; % 5
                1.000   0.0000  0.4500  0.2500  0.0000  0.0800  0.0600  0.0400; % 6
                1.000   0.0000  0.0000  0.0000  0.0000  0.0000  0.0000  0.0000; % 7
                1.000   0.2000  0.1000  0.0500  0.0600  0.0000  0.0000  0.0000; % 8
]

%% Six pulse converter -- see document andreas 
%% TL Lights - https://www.dranetz.com/wp-content/uploads/2014/02/harmonics-understanding-thefacts-part3.pdf

%% generator data
%   bus     Pg      Qg      Qmax    Qmin    Vg      mBase   status  Pmax    Pmin
mpc.gen = [
    1       0.0     0.0     200.00  -200.00 1.05    100.0   1       400.00  0.0;
];

%% branch data
%   f_bus	t_bus	r	    x	    b	    rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [ 
    2       3       ; %% TODO FRE
    4       5       ;
    6       7       ;
];

%% transformer data
%column_names%  f_bus   t_bus   xsc     r1      r2      vg      gnd1    gnd2    re1     xe1     re2     xe2   
mpc.xfmr = [
                1       2       XXX     XXX     XXX     'Yd11'  0       0       0.0     0.0     0.0     0.0;
                3       4       XXX     XXX     XXX     'Yy0'   0       0       0.0     0.0     0.0     0.0;
                5       6       XXX     XXX     XXX     'Dy11'  0       1       XXX     XXX     XXX     XXX;
                7       8       XXX     XXX     XXX     'Dz0'   0       1       XXX     XXX     XXX     XXX;
]
