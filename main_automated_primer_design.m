%
% Main file for the automated primer design
%
clear all, close all, clc

%Choose organism/genome
organism = 'Saccharomyces_cerevisiae';

%Choose method (e.g. C-terminal, N-terminal, deletion):
method = 'C-terminal';
%method = 'C-terminal_with_marker';
%method = 'N-terminal';

%Write primers to a single .CSV file (CSC = comma-separated values/characters)
if strcmpi(method,'C-terminal_with_marker')
    filename_out1 = [datestr(now,'yyyymmdd_THHMMSS'),'_C-terminal_with_marker_gene_tagging_primers.csv']; display(filename_out1)
elseif strcmpi(method,'C-terminal')
    filename_out1 = [datestr(now,'yyyymmdd_THHMMSS'),'_C-terminal_gene_tagging_primers.csv']; display(filename_out1)
elseif strcmpi(method,'N-terminal')
    filename_out1 = [datestr(now,'yyyymmdd_THHMMSS'),'_N-terminal_gene_tagging_primers.csv']; display(filename_out1)
else
end

fid1 = fopen(filename_out1, 'w','n','US-ASCII'); %write mode
fprintf(fid1, '%s,%s,%s\n','Name','Sequence','Notes');
fclose(fid1);

%List with yeast gene names
%gene_list = {'YER184C'}
%gene_list = {'YRF1-2'}; %gene is too close to the telomere region
%gene_list = {'RSP5'}
%gene_list = {'PRE1'}
gene_list = {'CRM1'}
%gene_list = {'NOP56'}
%gene_list = {'TDH1';'YRF1-2';'TOG1';'CRZ1'}
%gene_list = {'CSC1'}


fid1 = fopen(filename_out1, 'a','n','US-ASCII'); %append mode

%Design primers
for i = 1:size(gene_list,1) %loop over all genes in the gene list
        
    [g,h] = findGeneOfInterest(organism,gene_list{i}); 
    
    if ~isempty(g) && ~isempty(h) %gene entry was found        
        primer = designPrimers_Saccharomyces_cerevisiae(method,g,h); %Saccharomyces cerevisiae
    else
        primer = addEmptyPrimerEntries(gene_list{i}); %gene entry was not found
    end
            
    for j = 1:5 %write the five primers to the CSV file        
        fprintf(fid1, '%s,%s,%s\n',[primer(j).name],primer(j).seq,['Length = ',num2str(primer(j).length_total),' nt ; Tm = ', num2str(sprintf('%.2f',primer(j).Tm)),' C']);            
    end    
    
end

fclose(fid1); 
