% Erik Schluter
% 8/10/2015
% Shumway and Stoffer - 7.12

function fmriAnalysis()

    Shock = load('datasets\fmri_Awake_Shock.txt')';
    Heat  = load('datasets\fmri_Awake_Heat.txt')';
    t     = Shock(1,:);
    sdata = Shock(2:end,:); % separate time row and data
    hdata = Heat(2:end,:);
    
    % run PCA on shock data
    [s_PCs, s_var, s_pVar] = PCA(sdata);
    
    % run PCA on heat data
    [h_PCs, h_var, h_pVar] = PCA(hdata);
    
    % Use the matlab implementation of PCA and compare for kicks
    [m_sPCs, ~, m_sVar] = princomp(sdata');
    [m_hPCs, ~, m_hVar] = princomp(hdata');
    
    disp(['Shock data percentage variance explained = ' mat2str(s_pVar)])
    disp(['Heat data percentage variance explained  = ' mat2str(h_pVar)])
    
return