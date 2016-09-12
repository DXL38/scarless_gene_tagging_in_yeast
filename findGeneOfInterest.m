function [g, h] = findGeneOfInterest(organism,gene_name)
%
% FINDGENEOFINTERST finds the gene with gene_name and loads the DNA sequence and other
% metadata into MATLAB from the genome.gb file, which is stored locally.
%

%Load the yeast genome into matlab
switch organism
    case 'Saccharomyces_cerevisiae'
        R = rdir2('~/Lab/DNA/genomes/Saccharomyces_cerevisiae_S288c/*.gb'); %budding yeast
    case 'Schizosaccharomyces_pombe'
        R = rdir2('~/Lab/DNA/genomes/Schizosaccharomyces_pombe/*.gb'); %fission yeast
    otherwise
        display('Organismus is not implemented yet.')
end
            
%Search the standard name entries (e.g. ASH1)
for i = 1:size(R,1) %loop over all chromosomes
    g = genbankread(R(i).name);    
    for j = 1:size(g.CDS,2) %loop over all genes (of given chromosome)                
        if strcmpi(g.CDS(j).gene,(gene_name)) == 1
            g.CDS(j)
            h = j; display(h) %index of the hit            
            filename = R(i).name; display(filename)            
            return        
        end        
    end
end


if exist('filename','var') == 0 %filename does not exist
%Search the systematic name entries (e.g. YKL185W)
    display('Search the systematic name entries:')
    search_str = ['\/locus_tag','="',gene_name,'"'];
    display(search_str)
    for i = 1:size(R,1) %loop over all chromosomes
        g = genbankread(R(i).name);    
        for j = 1:size(g.CDS,2) %loop over all genes (of given chromosome)                                                        
            text = g.CDS(j).text;            
            for ii = 1:size(text,1) %loop over all the text lines                          
                if regexp(text(ii,:),search_str) == 1; %idx of match should be 1                                        
                    g.CDS(j)
                    h = j; display(h) %index of the hit            
                    filename = R(i).name; display(filename)                       
                    g.CDS(h).gene = gene_name; %use the systematic name as the gene name                    
                    return        
                end
            end                
        end
    end

else
end

if exist('filename')==1
    g = genbankread(filename);
else
    g='';h='';
end

