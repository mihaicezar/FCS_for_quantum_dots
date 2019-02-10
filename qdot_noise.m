%% Simulation of a single quantum dot based on a rate eq. model
% capable of noise calculation
% based on Kaasbjerg and Belzig, PRB 91 235413 (2015)
% Authors: Cezar Harabula (general case) and Gergo Fulop (1 excited level)

% This script uses psin.m from: https://www.mathworks.com/matlabcentral/fileexchange/
% 5687-polygamma-function-of-arbitrary-order-valid-complex-plane

clc;

% Spectrum and number of electrons in each orbital, to be defined in config file (qdot_config.m)
E_spectrum = [0];
E_electrons = [0];

qdot_config; % gives E_spectrum, E_electrons, num_states, Gamma*, kT...

hbar = 1;
chi_step=1e-3;
im=sqrt(-1);

E=E_spectrum;

same_charge   = @(i) find(E_electrons==E_electrons(i)  );
previous_charge = @(i) find(E_electrons==E_electrons(i)-1);
next_charge   = @(i) find(E_electrons==E_electrons(i)+1);
E_same_charge   = cellfun(same_charge,  num2cell(1:num_states),'UniformOutput',false);
E_previous    = cellfun(previous_charge,num2cell(1:num_states),'UniformOutput',false);
E_next   = cellfun(next_charge,  num2cell(1:num_states),'UniformOutput',false);

current = zeros(length(gateVv),length(biasVv));
current_old = zeros(length(gateVv),length(biasVv));
current_noise = zeros(length(gateVv),length(biasVv));
error_map = zeros(length(gateVv),length(biasVv));
trigamma_up = zeros(length(gateVv),length(biasVv));
pmatrix = zeros(num_states,length(gateVv)); % for storing the solution of the stationary p vector (zero bias cut)

Bose = @(en) 1/(exp(en/kT)-1);

I = @(E1,E2,Ea,Eb) Bose(E2-E1)./(Ea-Eb) .* real(...
  psin(  1/2 + im./(2*pi*kT) .* (E2-Ea))  ...
  -psin( 1/2 - im./(2*pi*kT) .* (E2-Eb)) ...
  -psin( 1/2 + im./(2*pi*kT) .* (E1-Ea))  ...
  +psin( 1/2 - im./(2*pi*kT) .* (E1-Eb)));

J = @(E1,E2,Ea)  1./(2*pi*kT) .*  Bose(E2-E1) .* imag(...
  psin(1, 1/2 + im./(2*pi*kT) .* (E2-Ea))  ...
  -psin(1, 1/2 + im./(2*pi*kT) .* (E1-Ea)));

mu = @(m,n) E(max(m,n))-E(min(m,n)); % because E spectrum is ordered upon charge state, this correctly returns negative result when the input consists of E(N,i)>E(N+1,j)

tic

for Vbi = 1:length(biasVv)
  biasV = biasVv(Vbi);
  parfor Vgi = 1:length(gateVv)
    % matrix of tunneling rates; convention: g(n,m) is the tunneling rate from state m to state n (as opposed to Belzig)
    g = zeros(num_states);
    g_chi = zeros(num_states);
    g_chi_minus = zeros(num_states);
    
    % matrices of cotunneling rates; same sense convention
    gCotLL = zeros(num_states);
    gCotRR = zeros(num_states);
    gCotLR = zeros(num_states);
    gCotRL = zeros(num_states);
    
    gL= zeros(num_states);
    gR= zeros(num_states);
    gLe= zeros(num_states); % Le means "exits the QD thru L"
    gRe= zeros(num_states);
    
    % Fermi functions, symmetric biasing:
    dEL=gateVv(Vgi)-biasV./2;
    dER=gateVv(Vgi)+biasV./2;
    
    FermiL = @(x)(exp((x-dEL)./kT)+1).^(-1);
    FermiR = @(x)(exp((x-dER)./kT)+1).^(-1);
    
    % for elastic cotunneling:
    ILR = @(Ea,Eb) I(dEL,dER,Ea,Eb);
    IRL = @(Ea,Eb) I(dER,dEL,Ea,Eb);
    ILL = @(Ea,Eb) I(dEL,dEL,Ea,Eb); % not in Gergo's code
    IRR = @(Ea,Eb) I(dER,dER,Ea,Eb); % not in Gergo's code
    
    JLR = @(Ea) J(dEL,dER,Ea);
    JRL = @(Ea) J(dER,dEL,Ea);
    JLL = @(Ea) J(dEL,dEL,Ea); % not in Gergo's code
    JRR = @(Ea) J(dER,dER,Ea); % not in Gergo's code
    
    % for inelastic cotunneling:
    
    % Delta = excitation energy = E(N,m)-E(N,n) with m,n orbital indices
    ILL_up = @(Ea,Eb,Delta) inelastic_cotunneling* I(dEL,dEL+Delta,Ea,Eb);
    IRR_up = @(Ea,Eb,Delta) inelastic_cotunneling* I(dER,dER+Delta,Ea,Eb);
    ILR_up = @(Ea,Eb,Delta) inelastic_cotunneling* I(dEL,dER+Delta,Ea,Eb);
    IRL_up = @(Ea,Eb,Delta) inelastic_cotunneling* I(dER,dEL+Delta,Ea,Eb);
    
    JLL_up = @(Ea,Delta) 0;%J(dEL,dEL+Delta,Ea);
    JRR_up = @(Ea,Delta) 0;%J(dER,dER+Delta,Ea);
    JLR_up = @(Ea,Delta) 0;%J(dEL,dER+Delta,Ea);
    JRL_up = @(Ea,Delta) 0;%J(dER,dEL+Delta,Ea);
    
    % Sequential tunneling
    % transition rates:
    
    for i = 1:num_states
      
      for ku = 1:length(E_next{i})
        u = E_next{i}(ku);
        % different charge states, indices i and u (i<u)
        % (filling the dot)
        gL(u,i)= GammaL/hbar .* abs(ML(u,i))^2 .* FermiL(mu(u,i)) ;
        gR(u,i)= GammaR/hbar .* abs(MR(u,i))^2 .* FermiR(mu(u,i));
        g(u,i) = gR(u,i);
        g_chi(u,i) = gL(u,i);
      end
      
      for kd = 1:length(E_previous{i})
        d = E_previous{i}(kd);
        % different charge states, indices i and d (i>d)
        % (emptying the dot)
        gLe(d,i) = GammaL/hbar .* abs(ML(d,i))^2 .* (1-FermiL(mu(d,i)));
        gRe(d,i) =  GammaR/hbar .* abs(MR(d,i))^2 .*(1-FermiR(mu(d,i)));
        g(d,i) = gRe(d,i);
        g_chi_minus(d,i) = gLe(d,i);
      end
      
    end
    
    % excitation/relaxation (i -> j or j -> i respectively )
    
    for i=1:num_states-1
      for j=i+1:num_states % j is always >i
        if(E_electrons(j)==E_electrons(i)) % tests if j is among the orbitals of i
          g(i,j)=Gamma_rel_mat(i,j)/hbar *abs(Bose(E(i)-E(j)));
          g(j,i)=Gamma_rel_mat(j,i)/hbar *abs(Bose(E(j)-E(i)));
          % no transport, thus g_chi(_minus) do not change
        end
      end
    end
    
    % rate equation matrix:
    
    % diagonal elements
    for i=1:num_states
      
      % diminishing, from emptying
      for j = 1:length(E_previous{i})
        d = E_previous{i}(j);
        g(i,i) = g(i,i) - GammaL/hbar * abs(ML(d,i))^2 .* (1 - FermiL(E(i)-E(d)));
        g(i,i) = g(i,i) - GammaR/hbar * abs(MR(d,i))^2 .* (1 - FermiR(E(i)-E(d)));
      end
      
      % from filling
      for j = 1:length(E_next{i})
        u = E_next{i}(j);
        g(i,i) = g(i,i) - GammaL/hbar * abs(ML(u,i))^2 .* FermiL(E(u)-E(i));
        g(i,i) = g(i,i) - GammaR/hbar * abs(MR(u,i))^2 .* FermiR(E(u)-E(i));
      end
      
      % from relaxation / thermal excitation
      for j = 1:length(E_same_charge{i})
        s = E_same_charge{i}(j);
        if(i~=s)
          % subtract relax/excitation from i to s
          g(i,i) = g(i,i) - g(s,i);
          % no transport, thus g_chi(_minus) do not change
        end
      end
      
    end

    % Elastic co-tunnneling:
    for i=1:num_states
      
      for j=1:length(E_previous{i})
        d = E_previous{i}(j);
        
        gCotLR(i,i) = gCotLR(i,i) + abs(ML(i,d)*MR(d,i))^2 * JLR(mu(i,d)); % eq.(45): |term|^2
        gCotRL(i,i) = gCotRL(i,i) + abs(MR(i,d)*ML(d,i))^2 * JRL(mu(i,d));
        gCotLL(i,i) = gCotLL(i,i) + abs(ML(i,d)*ML(d,i))^2 * JLL(mu(i,d)); % not in Gergo's code...
        gCotRR(i,i) = gCotRR(i,i) + abs(MR(i,d)*MR(d,i))^2 * JRR(mu(i,d));
        
        for j2=1:j-1
          d2 = E_previous{i}(j2);
          gCotLR(i,i) = gCotLR(i,i) + 2 * real( ML(i,d2)*MR(d2,i) * conj(ML(i,d)*MR(d,i)) * ILR(mu(i,d2),mu(i,d)));% eq.(45): products
          gCotRL(i,i) = gCotRL(i,i) + 2 * real( MR(i,d2)*ML(d2,i) * conj(MR(i,d)*ML(d,i)) * IRL(mu(i,d2),mu(i,d)));
          gCotLL(i,i) = gCotLL(i,i) + 2 * real( ML(i,d2)*ML(d2,i) * conj(ML(i,d)*ML(d,i)) * ILL(mu(i,d2),mu(i,d)));
          gCotRR(i,i) = gCotRR(i,i) + 2 * real( MR(i,d2)*MR(d2,i) * conj(MR(i,d)*MR(d,i)) * IRR(mu(i,d2),mu(i,d)));
        end
      end
      
      for j=1:length(E_next{i})
        u = E_next{i}(j);
        gCotLR(i,i) = gCotLR(i,i) + abs(MR(i,u)*ML(u,i))^2 * JLR(mu(i,u)); % eq.(42): |term|^2
        gCotRL(i,i) = gCotRL(i,i) + abs(ML(i,u)*MR(u,i))^2 * JRL(mu(i,u));
        gCotLL(i,i) = gCotLL(i,i) + abs(ML(i,u)*ML(u,i))^2 * JLL(mu(i,u));
        gCotRR(i,i) = gCotRR(i,i) + abs(MR(i,u)*MR(u,i))^2 * JRR(mu(i,u));
        
        for j2=1:j-1
          u2 = E_next{i}(j2);
          gCotLR(i,i) = gCotLR(i,i) + 2 * real(MR(i,u2)*ML(u2,i) * conj(MR(i,u)*ML(u,i)) * ILR(mu(i,u2),mu(i,u))); % eq.(42): products
          gCotRL(i,i) = gCotRL(i,i) + 2 * real(ML(i,u2)*MR(u2,i) * conj(ML(i,u)*MR(u,i)) * IRL(mu(i,u2),mu(i,u)));
          gCotLL(i,i) = gCotLL(i,i) + 2 * real(ML(i,u2)*ML(u2,i) * conj(ML(i,u)*ML(u,i)) * ILL(mu(i,u2),mu(i,u)));
          gCotRR(i,i) = gCotRR(i,i) + 2 * real(MR(i,u2)*MR(u2,i) * conj(MR(i,u)*MR(u,i)) * IRR(mu(i,u2),mu(i,u)));
        end
        
        % now combine up and down states to get the negative reals - eqs.(43,44)
        for k=1:length(E_previous{i})
          d = E_previous{i}(k);
          gCotLR(i,i) = gCotLR(i,i) - 2 * real(ML(i,d)*MR(d,i) * conj(MR(i,u)*ML(u,i)) * ILR(mu(i,d),mu(i,u)));
          gCotRL(i,i) = gCotRL(i,i) - 2 * real(MR(i,d)*ML(d,i) * conj(ML(i,u)*MR(u,i)) * IRL(mu(i,d),mu(i,u)));
          gCotLL(i,i) = gCotLL(i,i) - 2 * real(ML(i,d)*ML(d,i) * conj(ML(i,u)*ML(u,i)) * ILL(mu(i,d),mu(i,u)));
          gCotRR(i,i) = gCotRR(i,i) - 2 * real(MR(i,d)*MR(d,i) * conj(MR(i,u)*MR(u,i)) * IRR(mu(i,d),mu(i,u)));
        end
      end
      
      gCotLR(i,i) = 1/(2*pi*hbar) * GammaL * GammaR * gCotLR(i,i);
      gCotRL(i,i) = 1/(2*pi*hbar) * GammaL * GammaR * gCotRL(i,i);
      gCotLL(i,i) = 1/(2*pi*hbar) * GammaL * GammaL * gCotLL(i,i);
      gCotRR(i,i) = 1/(2*pi*hbar) * GammaR * GammaR * gCotRR(i,i);
      
      g(i,i) = g(i,i) - gCotLR(i,i) - gCotRL(i,i);% - gCotLL(i,i)- gCotRR(i,i);
      g_chi(i,i) = g_chi(i,i) + gCotLR(i,i);
      g_chi_minus(i,i) = g_chi_minus(i,i) + gCotRL(i,i);
      
    end
    
    % Inelastic cotunneling
    for i=1:num_states
      for j=1:length(E_same_charge{i})
        s = E_same_charge{i}(j);
        
        if(i==s)
          % this is the elastic case, already calculated
          continue;
        end
        
        % now, gCotXX(s,i)=0 from initialization
        
        if(i<s) % transition to a higher energy
          
          Delta = E(s)-E(i);
          
          for kd=1:length(E_previous{i})
            d = E_previous{i}(kd);
            
            gCotLL(s,i) = gCotLL(s,i) + abs(ML(s,d)*ML(d,i))^2 * JLL_up(mu(s,d),Delta);
            gCotRR(s,i) = gCotRR(s,i) + abs(MR(s,d)*MR(d,i))^2 * JRR_up(mu(s,d),Delta);
            gCotLR(s,i) = gCotLR(s,i) + abs(ML(s,d)*MR(d,i))^2 * JLR_up(mu(s,d),Delta); % eq.(46): |first_MM_term|^2
            gCotRL(s,i) = gCotRL(s,i) + abs(MR(s,d)*ML(d,i))^2 * JRL_up(mu(s,d),Delta);
            
          end
          
          for ku=1:length(E_next{i})
            u = E_next{i}(ku);
            
            gCotLL(s,i) = gCotLL(s,i) + abs(ML(s,u)*ML(u,i))^2 * JLL_up(mu(u,i),Delta);
            gCotRR(s,i) = gCotRR(s,i) + abs(MR(s,u)*MR(u,i))^2 * JRR_up(mu(u,i),Delta);
            gCotLR(s,i) = gCotLR(s,i) + abs(MR(s,u)*ML(u,i))^2 * JLR_up(mu(u,i),Delta); % eq.(46): |2nd_MM_term|^2
            gCotRL(s,i) = gCotRL(s,i) + abs(ML(s,u)*MR(u,i))^2 * JRL_up(mu(u,i),Delta);
            
          end
          
          for kd=1:length(E_previous{i})
            d = E_previous{i}(kd);
            for ku=1:length(E_next{i})
              u = E_next{i}(ku);
              
              gCotLL(s,i) = gCotLL(s,i) - 2 * real(ML(s,d)*ML(d,i) * conj(ML(s,u)*ML(u,i)) * ILL_up(mu(s,d),mu(u,i),Delta));
              gCotRR(s,i) = gCotRR(s,i) - 2 * real(MR(s,d)*MR(d,i) * conj(MR(s,u)*MR(u,i)) * IRR_up(mu(s,d),mu(u,i),Delta));
              gCotLR(s,i) = gCotLR(s,i) - 2 * real(ML(s,d)*MR(d,i) * conj(MR(s,u)*ML(u,i)) * ILR_up(mu(s,d),mu(u,i),Delta)); % eq.(46): product
              trigamma_up(Vgi,Vbi) = IRL_up(mu(s,d),mu(u,i),Delta); % for debugging
              gCotRL(s,i) = gCotRL(s,i) - 2 * real(MR(s,d)*ML(d,i) * conj(ML(s,u)*MR(u,i)) * IRL_up(mu(s,d),mu(u,i),Delta));
            end
          end
        end
        
        if(i>s) % transition to a lower energy
          Delta = -E(s)+E(i);
          
          for kd=1:length(E_previous{i})
            d = E_previous{i}(kd);
            
            gCotLL(s,i) = gCotLL(s,i) + abs(ML(s,d)*ML(d,i))^2 * JLL_up(mu(s,d),-Delta);
            gCotRR(s,i) = gCotRR(s,i) + abs(MR(s,d)*MR(d,i))^2 * JRR_up(mu(s,d),-Delta);
            gCotLR(s,i) = gCotLR(s,i) + abs(ML(s,d)*MR(d,i))^2 * JLR_up(mu(s,d),-Delta); % eq.(47): |first_MM_term|^2
            gCotRL(s,i) = gCotRL(s,i) + abs(MR(s,d)*ML(d,i))^2 * JRL_up(mu(s,d),-Delta);
            
          end
          
          for ku=1:length(E_next{i})
            u = E_next{i}(ku);
            
            gCotLL(s,i) = gCotLL(s,i) + abs(ML(s,u)*ML(u,i))^2 * JLL_up(mu(u,i),-Delta);
            gCotRR(s,i) = gCotRR(s,i) + abs(MR(s,u)*MR(u,i))^2 * JRR_up(mu(u,i),-Delta);
            gCotLR(s,i) = gCotLR(s,i) + abs(MR(s,u)*ML(u,i))^2 * JLR_up(mu(u,i),-Delta); % eq.(47): |2nd_MM_term|^2
            gCotRL(s,i) = gCotRL(s,i) + abs(ML(s,u)*MR(u,i))^2 * JRL_up(mu(u,i),-Delta);
            
          end
          
          for kd=1:length(E_previous{i})
            d = E_previous{i}(kd);
            for ku=1:length(E_next{i})
              u = E_next{i}(ku);
              
              gCotLL(s,i) = gCotLL(s,i) - 2 * real(ML(s,d)*ML(d,i) * conj(ML(s,u)*ML(u,i)) ...
                * ILL_up(mu(s,d),mu(u,i),-Delta));
              gCotRR(s,i) = gCotRR(s,i) - 2 * real(MR(s,d)*MR(d,i) * conj(MR(s,u)*MR(u,i)) ...
                * IRR_up(mu(s,d),mu(u,i),-Delta));
              gCotLR(s,i) = gCotLR(s,i) - 2 * real(ML(s,d)*MR(d,i) * conj(MR(s,u)*ML(u,i)) ...
                * ILR_up(mu(s,d),mu(u,i),-Delta)); % eq.(47): product
              gCotRL(s,i) = gCotRL(s,i) - 2 * real(MR(s,d)*ML(d,i) * conj(ML(s,u)*MR(u,i)) ...
                * IRL_up(mu(s,d),mu(u,i),-Delta));
            end
          end
        end
        
        gCotLL(s,i) = 1/(2*pi*hbar) * GammaL * GammaL * gCotLL(s,i);
        gCotRR(s,i) = 1/(2*pi*hbar) * GammaR * GammaR * gCotRR(s,i);
        gCotLR(s,i) = 1/(2*pi*hbar) * GammaL * GammaR * gCotLR(s,i);
        gCotRL(s,i) = 1/(2*pi*hbar) * GammaR * GammaL * gCotRL(s,i);
        
        g(s,i) = g(s,i) + gCotLL(s,i) ...
          + gCotRR(s,i) ;
        g_chi(s,i) = g_chi(s,i) + gCotLR(s,i);
        g_chi_minus(s,i) = g_chi_minus(s,i) + gCotRL(s,i);
        
        g(i,i) = g(i,i) - gCotLL(s,i) ...
          - gCotRR(s,i) ...
          - gCotLR(s,i) ...
          - gCotRL(s,i);
      end
    end
    
    for chi = [+chi_step -chi_step 0] % in this specific order
      
      P = g + g_chi*exp(im*chi) + g_chi_minus*exp(-im*chi);
      
      if (any(isnan(P(:))))
        sre = nan;
        error_map(Vgi,Vbi)=1;
      else
        if (any(isinf(P(:))))
          1
        end
        % check for inf and replace any infinity by 1e308 (with its sign)
        P_real = real(P);
        P_real(isinf(P_real)) = sign(P_real(isinf(P_real))) * 1e308; % larger exponent x 0 = NaN
        P_imag = imag(P);
        P_imag(isinf(P_imag)) = sign(P_imag(isinf(P_imag))) * 1e308; % larger exponent x 0 = NaN
        P = P_real + 1i * P_imag;
        
        % diagonalize
        [~,D] = eig(P);
        D=diag(D);
        [~,IND] = min(abs(real(D)));
        sre = D(IND) ;
      end
      
      if (chi > 0)
        % + chi_step;
        sre_plus = sre ;
        
      elseif (chi < 0)
        % -chi_step;
        sre_minus = sre ;
        
      elseif (chi == 0)
        sre_zero = sre ;
        
        % Old method: should be recovered for chi=0
        pdot_stationary = zeros(num_states,1); %the time derivative of the probability vector is nil mathematical trick: setting one pdot to 1 and a row in P to 1, we state that the sum of the probabilities is 1 (this trick makes the following equation solvable)
        pdot_stationary(1) = 1; %[1 0 0 0]';
        P(1,:) = ones(1,num_states);
        % solve the stationary probability vector, p
        p = P\pdot_stationary ;
        
        
        if (biasV<1e-4)
          pmatrix(:,Vgi) = p;
        end
        % calculating the current:
        crt=0;
        for i=1:num_states
          
          for ku = 1:length(E_next{i})
            u = E_next{i}(ku);
            crt = crt + p(i)*gL(u,i);
          end
          
          for kd = 1:length(E_previous{i})
            d = E_previous{i}(kd);
            crt = crt + p(i)*(-gLe(d,i));
          end
        end
        current_old(Vgi,Vbi) = -1 * crt;
      end
      
      
    end % chi end (calculated stuff for both +chi_step and -chi_step)
    %current(Vgi,Vbi) = (sre_plus)./(im*2*chi_step);
    current(Vgi,Vbi) = -(sre_plus-sre_minus)./(im*2*chi_step);
    current_noise(Vgi,Vbi) = (sre_plus-2*sre_zero+sre_minus)./((im*chi_step)^2);
  end
end

toc

%% display figures and save output

qdot_noise_display;

mkdir(config_dir);
save([config_dir '/qdot_noise.mat'],'current', 'current_noise', 'gateVv', 'biasVv','config_dir');
copyfile('qdot_config.m', [config_dir '/' config_dir '.m']);
savefig(fighandle, [config_dir '/' config_dir '.fig'], 'compact');
