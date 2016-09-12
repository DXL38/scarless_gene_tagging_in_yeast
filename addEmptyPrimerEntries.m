function primer = addEmptyPrimerEntries(gene_name)
%
% ADDEMPTYPRIMERENTRIES adds an empty primer entry.
% (e.g. for genes which termini are too close to the telomers)
%

%primer 1
primer(1).name=[gene_name,'_F1']; 
primer(1).seq=''; 
primer(1).length_total=''; 
primer(1).Tm='';

%primer 2
primer(2).name=[gene_name,'_R1']; 
primer(2).seq=''; 
primer(2).length_total=''; 
primer(2).Tm='';

%primer 3
primer(3).name=[gene_name,'_F2']; 
primer(3).seq=''; 
primer(3).length_total=''; 
primer(3).Tm='';

%primer 4
primer(4).name=[gene_name,'_R2']; 
primer(4).seq=''; 
primer(4).length_total=''; 
primer(4).Tm='';

%primer 5
primer(5).name=[gene_name,'_Fc']; 
primer(5).seq=''; 
primer(5).length_total=''; 
primer(5).Tm='';
