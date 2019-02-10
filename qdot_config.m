% This file is initially named: qdot_config.m
% The present config and output will be saved in directory: 
config_dir = 'config_M2=0.4_M2e=0.6';

% GammaL, GammaR: coupling strength to the normal leads [eV]
% kT: temperature [meV]
% Ea: energy of the first level [meV]
% biasV: bias voltage [mV]; gateVv: gate voltage
% Uc: charging energy [meV]

inelastic_cotunneling = 1; %1 to include inelastic cotunneling; 0 otherwise

Gamma = 0.001;
GammaL = 1 * Gamma;
GammaR = 1 * Gamma;
Gamma_rel = 0.001 * Gamma; % relax rate

Ea = 0.5;
Uc = 1; 
kT = 0.004;

E_spectrum = [  0 ...
                Ea ...
               (2*Ea+Uc) (2*Ea+Uc)+0.1*Uc+1e3*eps (2*Ea+Uc)+0.35*Uc ...
               3*Ea+3*Uc ];
E_electrons = [0  1  2 2 2  3];

%%% Function orb will return the E_spectrum index of an electron, from the number of electrons and the the orbital number (e.g. orb (2,2) refers to 2 electrons, 2nd excited state).
subindex = @(vector,i) vector(i);
% N = charge; orbital index is 0 for a ground state
orb = @(N,orbital_index) subindex(find(E_electrons==N),orbital_index+1);
% Groundstate for N electrons
GS = @(N) orb(N,0);

num_states = length(E_spectrum);
% Matrix of relaxation rates
% The upper part contains relaxation rates, the lower part contains excitation rates
Gamma_rel_mat = Gamma_rel * ones(num_states);

% Tunneling matrix elements (all 1 by default)
% Convention for M(n,m): from state m to state n (as opposed to n->m in Belzig)
ML = ones(num_states);
ML(GS(1),GS(2)) = 0.4; ML(GS(2),GS(1)) = 0.4;
ML(GS(3),GS(2)) = 0.4; ML(GS(2),GS(3)) = 0.4;
ML(orb(2,1),GS(1))=0.6; ML(GS(1),orb(2,1))=0.6;
ML(orb(2,1),GS(3))=0.6; ML(GS(3),orb(2,1))=0.6;
MR=ML;

biasVv = linspace(-1.2,1.2,120+1);
gateVv = linspace(1.3,2.7,140+1);