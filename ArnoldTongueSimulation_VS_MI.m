% Calculate vector strength in a deterministic model
% This uses the toy hair-bundle model from 2012 PNAS.
% jsalvi@rockefeller.edu

% Choose a stiffness range
ke = linspace(1.5,4,10);
% Choose noise levels
noiselevel = [0.05 0.1 0.2];
% Chose a force range
Fe = 0;
% Time vector
t = linspace(0,5e2,5e2);
% Choose amplitudes
ampl = linspace(0,1,500);
% Choose frequencies
freq = linspace(0,1,500);

% Generate data
for i = 1:length(ke)
    for j = 1:length(noiselevel)
        for k = 1:length(Fe)
            for l = 1:length(ampl)
                for m = 1:length(freq)
                    [Xdet{i,k,l,m},Xsto{i,j,k,l,m},Fext{l,m}] = hbtoymodel(Fe(k),ke(i),noiselevel(j),ampl(l),freq(m),t);
                end
            end
        end
    end
end

% Calculate vector strength in deterministic and stochastic cases
for i = 1:length(ke)
    for j = 1:length(noiselevel)
        for k = 1:length(Fe)
            for l = 1:length(ampl)
                for m = 1:length(freq)
                    [VSdet{i,k}(l,m), Rayleighdet{i,k}(l,m)] = vscalc2(Xdet{i,k,l,m}(1,floor(length(t)/2):end),Fext{l,m}(floor(length(t)/2):end),1);
                    [VSsto{i,j,k}(l,m), Rayleighdet{i,j,k}(l,m)] = vscalc2(Xsto{i,j,k,l,m}(1,floor(length(t)/2):end),Fext{l,m}(floor(length(t)/2):end),1);
                    %[~,pIdispldet{i,k}(l,m),~,Idispldet{i,k}(l,m)] = mutualinfostatkde(Xsto{i,j,k,l,m}(1,floor(length(t)/2):end),Fext{l,m}(floor(length(t)/2):end),[],1,0.01,1000,1,0,1);
                    %[~,pIdispldet{i,k}(l,m),~,Idispldet{i,k}(l,m)] = mutualinfostatkdephase(Xsto{i,j,k,l,m}(1,floor(length(t)/2):end),Fext{l,m}(floor(length(t)/2):end),[],1,0.01,1000,1,0,1);
                end
            end
        end
    end
end

                
