%% Authored by Joshua Soneson 2018
function[u,U,Ppos,Pneg,I_td] = TDNL_mat(u,U,KK,JJ,c,cutoff,Ppos,Pneg,I_td)
% converts spectrum to one cycle of the time-domain waveform and
% integrates the invicid Burger's equation using upwind/downwind
% method with periodic boundary conditions. TDNL stands for Time
% Domain NonLinear.
% set peak and trough values to zero; in case they are not assigned later
% nonlinear case - enter loop



if(KK==1)	% linear case - do nothing
    for jj=1:JJ
        Ppos(jj) = abs(u(jj,1));
    end
    Pneg = -Ppos;
else		% nonlinear case - enter loop
    for jj=JJ:-1:1
        % execute nonlinear step only if amplitude is above cutoff;
        % row jj=1 is always computed so plots look nice
        I_td(jj) = 0;
        if(abs(u(jj,1)) < cutoff/20)	% if pressure is too low, skip nonlinear step
        else
            % convert from sin & cos representation to complex exponential
            U(2:KK+1) = conj(u(jj,:));
            U(2*KK:-1:KK+2) = u(jj,1:KK-1);
            U(1) = 0;
            % transform to time domain:
            U = KK*real(ifft(U));
            I_td(jj) = sum(U.^2);
            % determine how many steps necessary for CFL<1 (CFL<0.9 to be safe).
            PP = ceil(2*c*max(abs(U))/0.9);
            %if(j==1)
            %  plot(U);
            %  hold on
            %end
            % Nonlinear integration (upwind/downwind) algorithm
            for pp=1:PP

                % negels
                nU = U(U<0);
                Xn = zeros(size(nU));
                nei = find(U<0);

                if (min(nei)==1)
                    innern = nU.^2 - U(nei).^2;
                    Xn(1) = U(1) + c*(U(1)*U(1) - U(end)*U(end))/PP;
                    Xn = nU + c.*(innern)./PP;
                    Xn(1) = U(1) + c*(U(1)*U(1) - U(end)*U(end))/PP;

                else
                    innern = nU.^2 - U(nei-1).^2;
                    Xn = nU + c.*(innern)./PP;
                end


                % posels
                pU = U(U>0);
                Xp = zeros(size(pU));
                poi = find(U>0);
                Xp = zeros(length(pU)+1);
                if (max(poi)==length(U))
                    innerp = (U(poi).^2) - pU.^2;
                    Xp(end) = U(end) + c*(U(1)*U(1) - U(end)*U(end))/PP;

                    Xp = pU + c*(innerp)/PP;
                    Xp(end) = U(end) + c*(U(1)*U(1) - U(end)*U(end))/PP;
                else
                    innerp = (U(poi+1).^2) - pU.^2;
                    Xp = pU + c*(innerp)/PP;
                end

                U(nei) = Xn;
                U(poi) = Xp;
            end
            X = U;
            % account for nonlinear losses:
            I_td(jj) = I_td(jj) - sum(X.^2);
            % transform back to frequency domain:
            Ppos(jj) = max(X);
            Pneg(jj) = min(X);
            X = fft(X)/KK;
            % convert back to sin & cos representation:
            u(jj,:) = conj(X(2:KK+1));
        end
    end
end
end
