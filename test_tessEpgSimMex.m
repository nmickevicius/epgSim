useMex = true;

if ~useMex
    
    fid = fopen('F0.bin','rb');
    raw = fread(fid, inf, 'float');
    F0 = raw(1:2:end) + 1i.*raw(2:2:end);
    F0 = reshape(F0,3,[])';
    figure; plot(abs(F0));
    
else

    nrf = 3*50;
    nz = 50;
    tr = single(0.018);
    te = single(0.006);
    tau = single(0.0028);
    tgap = single(0.0025);
    garea = single(70.455); % mT/m * ms
    T1 = single(0.8);
    T2 = single(0.08);
    D = 0.0;
    fa = single(25*ones(nrf,1)*pi/180);
    phi = single(zeros(nrf,1));
    sprof = single(ones(nz,1));
    
    while(exist('F0','var'))
        clear F0;
    end
    tic;
    F0_full = tessEpgSimMex(fa, phi, sprof, tr, te, tau, garea, tgap, T1, T2, D, 0);
    toc
    tic;
    F0 = tessEpgSimMex(fa, phi, sprof, tr, te, tau, garea, tgap, T1, T2, D, 20);
    toc
    F0_full = reshape(F0_full,3,[])';
    F0 = reshape(F0,3,[])';
    figure; 
    subplot(1,4,1); hold on; plot(abs(F0_full(:,1))); hold on; plot(abs(F0(:,1)));
    subplot(1,4,2); hold on; plot(abs(F0_full(:,2))); hold on; plot(abs(F0(:,2)));
    subplot(1,4,3); hold on; plot(abs(F0_full(:,3))); hold on; plot(abs(F0(:,3)));
    subplot(1,4,4); plot(abs(F0-F0_full)); ylim([-0.5,0.5]);
    
end