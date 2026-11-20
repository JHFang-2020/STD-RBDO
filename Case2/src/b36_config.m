function c = b36_config()
%B36_CONFIG  Case 2 (beam under moving loads, v36) constants and protocol.
%   The physical problem is the Case 2 beam of c2_config:
%   simply supported Euler-Bernoulli beam, L = 12 m, b = 0.12 m, quadratic
%   Bernstein height with design d (3 control heights, [0.16,0.30] m),
%   17 Gaussian inputs xi (manufacturing 1:3, material 4:5, two moving-load
%   histories 6:17), full-DOF Newmark dynamics (HF = b36_hf, 80 elements,
%   dt = 1 ms).  Failure: max_{x,t} |sigma| >= 160 MPa.
%   RBDO: min expected mass (prop. to d1+d2+d3)  s.t.  Pf(d) <= 1 %.
c = struct();
c.version  = 'case2_v36';
c.hf       = c2_config();                     % HF simulator settings (unchanged)
c.L = c.hf.L; c.b = c.hf.b; c.E0 = c.hf.E0; c.rho = c.hf.rho;
c.sigmaH = c.hf.sigmaH; c.gammaD = c.hf.dampingGamma;
c.allow    = 160;                             % MPa
c.target   = 0.01;
c.lower    = [0.16 0.16 0.16];
c.upper    = [0.30 0.30 0.30];
c.massCoef = c.rho*c.b*c.L/3;                 % kg per metre of sum(d)
% populations (common random numbers for every method and seed)
c.Nwork    = 4096;  c.workSeed = 1222001;     % working population
c.Nval     = 16384; c.valSeed  = 1229001;     % independent validation bank
% working space-time grid shared by every surrogate (retained block)
c.nx   = 41;   c.xw = linspace(0,c.L,c.nx)';
c.dtw  = 0.005;
c.tLF  = 0:c.dtw:4.0;                         % LF integration grid (from rest)
c.tw   = 1.0:c.dtw:4.0;  c.ntw = numel(c.tw); % retained window (all peaks in [1.48,3.04] s)
c.twIdx  = round(c.tw/c.hf.dt) + 1;           % indices into the HF time grid
c.twInLF = round(c.tw/c.dtw) + 1;             % indices into c.tLF
% critical space-time windows (spatial max over the beam, 0.1 s in time)
c.wlen = 0.10; c.nW = round((c.tw(end)-c.tw(1))/c.wlen);
c.winOf = min(floor((c.tw - c.tw(1))/c.wlen + 1e-9) + 1, c.nW);
c.tWin  = c.tw(1) + c.wlen*((1:c.nW) - 0.5);
c.Kobs = 6;                                   % windows retained per HF history (observation block)
c.Kev  = 4;                                   % critical windows per trajectory event (Clark)
% simplified physics (LF): exact statics + J assumed sine modes
c.J   = 2;
c.ns  = 121; c.xs = linspace(0,c.L,c.ns)';    % load quadrature (trapezoid)
c.nq  = 121; c.xq = linspace(0,c.L,c.nq)';    % section quadrature for modal integrals
c.thetaActive = 0.20;                         % active histories: LF log margin < theta
c.thetaExt    = 0.35;                         % extended set (cached loads) of a full screening
c.dFullMove   = 0.004;                        % re-screen the full population after a larger design move
c.thetaTruth  = 0.25;                         % HF screening margin for the truth
end
