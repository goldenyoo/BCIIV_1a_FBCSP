% ----------------------------------------------------------------------- %
%    File_name: myLAP.m
%    Programmer: Seungjae Yoo                             
%                                           
%    Last Modified: 2020_01_26                           
%                                                            
 % ----------------------------------------------------------------------- %
%% 
function cnt_n = myLAP(cnt,nfo)

sap = {"F5",["F3","FC5"]; "F3",["F5","AF3","F1","FC3"]; "F1",["F3","Fz","FC1"]; "Fz",["F1","F2","FCz"]; "F2",["Fz","F4","FC2"]; ...
    "F4",["F2","AF4","F6","FC4"]; "F6",["F4","FC6"]; "FC5",["F5","FC3","CFC5"]; "FC3",["FC5","F3","FC1","CFC3"]; "FC1",["FC3","F1","FCz","CFC1"]; ...
    "FCz",["FC1","Fz","FC2","Cz"]; "FC2",["FCz","F2","FC4","CFC2"]; "FC4",["FC2","F4","FC6","CFC4"]; "FC6",["FC4","F6","CFC6"]; "CFC7",["CFC5","T7"]; ...
    "CFC5",["CFC7","FC5","CFC3","C5"]; "CFC3",["CFC5","FC3","CFC1","C3"]; "CFC1",["CFC3","FC1","C1"]; "CFC2",["FC2","CFC4","C2"]; "CFC4",["CFC2","FC4","CFC6","C4"]; ...
    "CFC6",["CFC4","FC6","CFC8","C6"]; "CFC8",["CFC6","T8"]; "T7",["CFC7","C5","CCP7"]; "C5",["T7","CFC5","C3","CCP5"]; "C3",["C5","CFC3","C1","CCP3"]; ...
    "C1",["C3","CFC1","Cz","CCP1"]; "Cz",["C1","FCz","C2","CPz"]; "C2",["Cz","CFC2","C4","CCP2"]; "C4",["C2","CFC4","C6","CCP4"]; "C6",["C4","CFC6","T8","CCP6"]; ...
    "T8",["C6","CFC8","CCP8"]; "CCP7",["T7","CCP5"]; "CCP5",["CCP7","C5","CCP3","CP5"]; "CCP3",["CCP5","C3","CCP1","CP3"]; "CCP1",["CCP3","C1","CP1"]; ...
    "CCP2",["C2","CCP4","CP2"]; "CCP4",["CCP2","C4","CCP6","CP4"]; "CCP6",["CCP4","C6","CCP8","CP6"]; "CCP8",["CCP6","T8"]; "CP5",["CCP5","CP3","P5"]; ...
    "CP3",["CP5","CCP3","CP1","P3"]; "CP1",["CP3","CCP1","CPz","P1"]; "CPz",["CP1","Cz","CP2","Pz"]; "CP2",["CPz","CCP2","CP4","P2"]; "CP4",["CP2","CCP4","CP6","P4"]; ...
    "CP6",["CP4","CCP6","P6"]; "P5",["CP5","P3"]; "P3",["P5","CP3","P1"]; "P1",["P3","CP1","Pz","PO1"]; "Pz",["P1","CPz","P2"]; "P2",["Pz","CP2","P4","PO2"]; ...
    "P4",["P2","CP4","P6"]; "P6",["P4","CP6"]; };

for i = 1:size(sap,1)
        index_er = find(nfo.clab == sap{i,1});        
        if isempty(index_er), error("Cant' find %s",sap{i,1}); end
        sap_n{i,1} = index_er;
        for j = 1:length(sap{i,2})
            index_sub = find(nfo.clab == sap{i,2}(j));            
            if isempty(index_sub), error("Cant' find %s",sap{i,2}(j)); end
            sap_n{i,2}(j) = index_sub;       
        end
end
%%

cnt_n = zeros(size(cnt,1),size(cnt,2));

for t = 1: size(cnt,2)
    
    for i = 1:size(sap_n,1)
        ind_er = sap_n{i,1};  
        V_er = cnt(ind_er,t);
        V_sub = 0;
        tmp = length(sap_n{i,2});
        for j = 1:length(sap_n{i,2})
            ind_sub = sap_n{i,2}(j);         
            V_sub = V_sub + cnt(ind_sub,t);
        end
        V_sub = V_sub/tmp;
        
        cnt_n(ind_er,t) = V_er - V_sub;
    end
end

end
% ----------------------------------------------------------------------- %
%                               EOF
% ----------------------------------------------------------------------- %

