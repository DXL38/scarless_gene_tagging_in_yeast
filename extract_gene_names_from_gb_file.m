filename_out1 = 'Output_file_of_gene_names_from_gb_files.txt';
fid1 = fopen(filename_out1, 'a','n','US-ASCII'); %append mode

R = rdir2('/Users/dirk/Lab/DNA/genomes/Saccharomyces_cerevisiae_S288c/*.gb'); %budding yeast

search_str = ['\/locus_tag="','[^ ]*','"']; %[^ ]* this is all characters except spaces

for i = 1:size(R,1)
%for i = 1
    g = genbankread(R(i).name);

    for j = 1:size(g.CDS,2) %loop over all genes j on chromosome i
        gene_name = g.CDS(j).gene;
        %display(gene_name)
        
        if isempty(gene_name)
            text = g.CDS(j).text;
            text_str = reshape(text',1,size(text,1)*size(text,2));
            [k1, k2]=regexp(text_str,search_str);
            gene_name=text_str(k1+12:k2-1);      
        end
        display(gene_name)
         fprintf(fid1, '%s\n',gene_name);            
    end

end

!cat Output_file_of_gene_names_from_gb_files.txt | wc -l
