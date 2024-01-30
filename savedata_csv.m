clear all

nVp = [11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34]; % gib hier alle Vp-Nummern ein, deren Daten du als .csv speichern willst

for iVp = 1:length(nVp) 
    
    data = [];
    VpNo = num2str(nVp(iVp));
    filename = ['EC' VpNo];
    load(filename);
    savename = ['EC' VpNo '.csv'];
    s = struct2table(data);
    writetable(s,savename);
end

