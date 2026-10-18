%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   3-node symmetric test case (recombination isolation experiment).     %
%                                                                       %
%        node 2 (load/harmonic unit)                                    %
%       /                                                               %
%   node 1 (gen, type 3)                                                %
%       \                                                               %
%        node 3 (load/harmonic unit)                                    %
%                                                                       %
%   Two IDENTICAL branches 1-2 and 1-3.  Nodes 2 and 3 host identical    %
%   harmonic units and are electrically interchangeable by symmetry.     %
%   The generator at node 1 reaches harmonic ground through its          %
%   short-circuit impedance (gen_sc); no line charging ⇒ no resonance.   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function mpc = case3_shhc
mpc.version = '2';
mpc.baseMVA = 100.0;
mpc.principle = 'absolute equality'

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	 3	 0.0	 0.0	 0.0	 0.0	 1	 1.00	 0.0	 100.0	 1	 1.10	 0.90;
	2	 1	 10.0	 5.0	 0.0	 0.0	 1	 1.00	 0.0	 100.0	 1	 1.10	 0.90;
	3	 1	 10.0	 5.0	 0.0	 0.0	 1	 1.00	 0.0	 100.0	 1	 1.10	 0.90;
];

%% bus harmonic data
%column_names%  ref_angle   standard
mpc.bus_harmonics = [
                0.0         'IEC61000-3-6:2008';       % 1  (gen / ref)
                0.0         'IEC61000-3-6:2008';       % 2  (harmonic unit)
                0.0         'IEC61000-3-6:2008';       % 3  (harmonic unit)
]

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin
mpc.gen = [
	1	 20.0	 10.0	 50.0	 -50.0	 1.0	 100.0	 1	 100.0	 0.0; % grid
];

%% generator cost data
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	 0.0	 0.0	 3	 0.0	 0.0	 0.0; % grid
];

%% generator short-circuit impedance data
%column_names%  rsc     xsc
mpc.gen_sc = [
	0.001	0.02;   % grid short-circuit impedance (Ssc = 5 GVA, XRr = 20)
]

%% branch data (two identical lines; b = 0)
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	 2	 0.01	 0.10	 0.0	 100	 100	 100	 0.0	 0.0	 1	 -360	 360;
	1	 3	 0.01	 0.10	 0.0	 100	 100	 100	 0.0	 0.0	 1	 -360	 360;
];
