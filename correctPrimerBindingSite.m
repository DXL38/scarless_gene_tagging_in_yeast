function primer = correctPrimerBindingSite(g,h,primer)
%
% CORRECTPRIMERBINDINGSITE makes sure that index_5p and index_3p have only one entry.
% If a primer has multiple binding sites (e.g. because of gene duplication), index_5p 
% and index_3p would have multile entries. The solution implement here is to check 
% with index_5p is closest to the 5p positon of the ORF and then use this index.
%


if size(primer.index_5p,2) == 1
else           
    %g.CDS(h).indices(1) %5p index of the ORF; it should be sufficient to just focus on the 5p index
	dmin = abs(primer.index_5p - g.CDS(h).indices(1)); %calculate absolute distance of 5p index of primer to 5p index of ORF
	k = find(dmin == min(dmin)); %find index of the entry that corresponds to the primer that is closest to the ORF
	primer.index_5p = primer.index_5p(k);
	primer.index_3p = primer.index_3p(k);            
end

