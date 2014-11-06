% Calculate vector strength in a deterministic model
% This uses the toy hair-bundle model from 2012 PNAS.
% jsalvi@rockefeller.edu

% Choose a stiffness range
ke = linspace(1.5,4,10);
%ke=3;
% Choose noise levels
noiselevel = [0.05 0.1 0.2];
% Chose a force range
Fe = 0;
% Time vector
tlengths=[1000];
for i = 1:length(tlengths)
    t{i} = linspace(0,1e2*tlengths(i),5e2*tlengths(i));
end
% Choose amplitudes
ampl = linspace(0,1,500);
%ampl = [0.01 0.1 0.5];
% Choose frequencies
%freq = [0.01 0.05 0.1];
freq = linspace(0,1,500);

% Generate data
for n = 1:length(t)
for i = 1:length(ke)
    for j = 1:length(noiselevel)
        for k = 1:length(Fe)
            for l = 1:length(ampl)
                for m = 1:length(freq)
                    [Xdet{i,k,l,m,n},Xsto{i,j,k,l,m,n},Fext{l,m,n}] = hbtoymodel(Fe(k),ke(i),noiselevel(j),ampl(l),freq(m),t{n});
                end
            end
        end
    end
end
end

%{
% Calculate vector strength in deterministic and stochastic cases
for n = 1:length(t)
for i = 1:length(ke)
    for j = 1:length(noiselevel)
        for k = 1:length(Fe)
            for l = 1:length(ampl)
                for m = 1:length(freq)
                    [VSdet{i,k}(l,m,n), Rayleighdet{i,k}(l,m,n)] = vscalc2(Xdet{i,k,l,m,n}(1,1:end),Fext{l,m,n}(1:end),1);
                    [VSsto{i,j,k}(l,m,n), Rayleighdet{i,j,k}(l,m,n)] = vscalc2(Xsto{i,j,k,l,m,n}(1,1:end),Fext{l,m,n}(1:end),1);
                    %[~,pIdispldet{i,k}(l,m),~,Idispldet{i,k}(l,m)] = mutualinfostatkde(Xsto{i,j,k,l,m}(1,1:end),Fext{l,m}(1:end),[],1,0.01,1000,1,0,1);
                    %[~,pIdispldet{i,k}(l,m),~,Idispldet{i,k}(l,m)] = mutualinfostatkdephase(Xsto{i,j,k,l,m}(1,1:end),Fext{l,m}(1:end),[],1,0.01,1000,1,0,1);
                end
            end
        end
    end
end
end
%}                
