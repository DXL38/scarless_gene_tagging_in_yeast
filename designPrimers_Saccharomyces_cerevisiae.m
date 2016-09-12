function primer = designPrimers_Saccharomyces_cerevisiae(method,g,h)
%
% Function for automated primer design. Three different methods are implemented:
% Method #1: C-terminal tagging with marker using e.g. pDML152 --> CM_F1, CM_R1, CM_F2, CM_R2, and CM_Fc
% Method #2: C-terminal tagging (scarless) using e.g. pDML219 or pDML223 --> C_F1, C_R1, C_F2, C_R2, and C_Fc
% Method #3: N-terminal tagging (scarless) using e.g. pDML190 or pDML222 --> N_F1, N_R1, N_F2, N_R2, and N_Fc
%
%
% Method #1: C-terminal tagging with marker using e.g. pDML152
% - design 4 primers (CM_F1, CM_R1, CM_F2, and CM_R2) to generate >300 bp upstream homology (H1) and downstream homology (H2)
% - 500-bp window for searching for F1 and R2 primers
% - try to design primers with GC clamps; F1 and R2 primers have always a GC clamp (because of the design of the search algorithm)
% - code calculates the chromosomal indices of the 5' and 3' positions of the primer binding sites 
% - code also caculates the size of the H1 and H2 homologies (= PCR products)
%
% Method #2: C-terminal tagging w/o marker (scarless) using e.g. pDML219 or pDML223
% - no disruption of downstream regulatory sequences
% - similar to method #1 but markerless tagging w/o a scar
%
% Method #3: N-terminal tagging (scarless) using e.g. pDML190 or pDML222
% - similar to method #2 but N-terminal
%
%
% Work flow of the algorithm:
% - start with 20-mer and extend the primer until Tm > 58 C
% - try to add GC clamp by extending the primer by one or two more nt
% - exclude bad primers (based on hairpin formation propensity)
%
%
% Variables used:
% g : structure with information about chromsome 
%     (e.g. g.Definition = Saccharomyces cerevisiae S288c chromosome VII, complete sequence.)
% h : gene index number on chromosome
%     --> get gene information: g.CDS(h)
%




switch method
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------ METHOD #1 ----- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'C-terminal_with_marker'
        
    %First determine whether gene is on the upper (sense) or lower (anti-sense) strand
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- DNA sense strand (gene has normal orientation) ---
    if g.CDS(h).indices(end) > g.CDS(h).indices(1)         
        display([g.CDS(h).gene, ' gene has normal orientation.'])    
        ORF = g.Sequence(g.CDS(h).indices(1) : g.CDS(h).indices(end)-3); %exclude stop codon
        stop_codon = g.Sequence(g.CDS(h).indices(end)-2 : g.CDS(h).indices(end)); display(['Stop codon: ',upper(stop_codon)])
        
        try
        %500-bp search window for CM F1 primer; search window starts 280 bp upstream of stop codon
        upstream = g.Sequence(g.CDS(h).indices(end)-782 : g.CDS(h).indices(end)-283); 
        %500-bp search window for CM R2 primer; search window starts 280 bp upstream of stop codon
        downstream = g.Sequence(g.CDS(h).indices(end)+281 : g.CDS(h).indices(end)+780);
        catch ME
            display(['Upstream or downstream sequence(s) are not long enough. Skipped primer design for gene ',g.CDS(h).gene,'.'])
            primer = addEmptyPrimerEntries(g.CDS(h).gene);
            return
        end

        % ---------- First primer (CM F1 primer) ----------
        display('----- Designing CM F1 primer -----')    
        k = regexp(upstream, '[c|g]{2}'); %find all GC clamps
        k = sort(k,'descend'); %start with the primer closest to stop codon (i.e. the insertion site) and move further upstream from there

        for i = 1:size(k,2) %loop over all potential primers
            pl = 20; %minimal primer length
            p = upstream(k(i)-pl+2 : k(i)+1); 
            Tm = calculate_Tm_like_PrimerSelect(p);
            
            while Tm < 58.0 %Tm of primer needs to be above 58.0 C
                pl = pl +1; %increase primer length by 1       
                p = upstream(k(i)-pl+2 : k(i)+1);
                Tm = calculate_Tm_like_PrimerSelect(p);
            end

            s = oligoprop(p); %calculate properties of primer
            if isempty(s.Hairpins)                        
                primer(1).name = [g.CDS(h).gene,'_CM_F1'];
                primer(1).seq = p;
                primer(1).Tm = Tm;                                
                primer(1).overhang = '';
                primer(1).length_binding = length(p);                        
                primer(1).length_overhang = 0;
                primer(1).length_total = primer(1).length_overhang + primer(1).length_binding;                
                primer(1).index_5p = regexp(g.Sequence,primer(1).seq);   
                primer(1).index_3p = primer(1).index_5p + primer(1).length_binding - 1;      
                orderfields(primer(1),{'name','seq','Tm','overhang','length_binding','length_overhang','length_total','index_5p','index_3p'});
                break
            else
            end      
        end
        primer(1) = correctPrimerBindingSite(g,h,primer(1));
        clear k i pl p Tm s
        %Note, not necessary to check whether F1 primer is larger than 60 nt since a primer with 35x As or 35x Ts has a Tm of 58.1436 C              


        % ---------- Second primer (CM R1 primer) ----------
        display('----- Designing CM R1 primer -----')    
        pl = 20; %minimal primer length
        p = g.Sequence(g.CDS(h).indices(end)-pl-3+1 : g.CDS(h).indices(end)-3);   
        Tm = calculate_Tm_like_PrimerSelect(p);
        
        while Tm < 58.0
            pl = pl +1; %increase primer lenght by 1
            p = g.Sequence(g.CDS(h).indices(end)-pl-3+1 : g.CDS(h).indices(end)-3);
            Tm = calculate_Tm_like_PrimerSelect(p);
        end

        %Try to include GC clamp if possible
        display(['Original CM R1 primer: ',seqrcomplement(p)])

        %Case #1: CM R1 primer has GC clamp
        if regexp(p(1:2), '^[gc]{2}')
            display('CM R1 primer has GC clamp.') % --> do nothing since primer has GC clamp

        %Case #2: CM R1 primer has G/C at first position but A/T at second position    
        elseif regexp(p(1:2), '^[gc][at]')
            display('CM R1 primer has G/C at first position but A/T at second position.')                
            p_new = g.Sequence(g.CDS(h).indices(end)-(pl+1)-3+1 : g.CDS(h).indices(end)-3); %increase primer lenght by 1        
            display(['CM R1 primer (+1): ',seqrcomplement(p_new)])
            
            if regexp(p_new(1), '^[gc]') % new CM R1 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %use old seqeunce; CM R1 primer has G/C at first position
            end

        %Case #3: CM R1 primer has A/T at first position but G/C at second last position    
        elseif regexp(p(1:2), '^[at][gc]')
            display('CM R1 primer has A/T at first position but G/C at second position.')               
            p_new = g.Sequence(g.CDS(h).indices(end)-(pl+2)-3+1 : g.CDS(h).indices(end)-3); %increase primer lenght by 2        
            display(['CM R1 primer (+2): ',seqrcomplement(p_new)])
            
            if regexp(p_new(1:2), '^[gc]{2}') % new CM R1 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(1:2), '^[gc][at]') % new CM R1 primer has G/C at first position and A/T at second first position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(1:2), '^[at][gc]') % new CM R1 primer has A/T at first position and G/C at second first position
                %increase lenght of OLD CM R1 primer by 1 (this is same as decreasing length of NEW CM R1 primer by 1)
                p_new = g.Sequence(g.CDS(h).indices(end)-(pl+1)-3+1 : g.CDS(h).indices(end)-3);                     
                display(['CM R1 primer (+1): ',seqrcomplement(p_new)])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else % new CM R1 primer has A/T at first and second position --> use old CM R1 primer
            end 

        %Case #4: CM R1 primer has A/T at first and second position
        elseif regexp(p(1:2), '^[at]{2}')        
            display('CM R1 primer has A/T at first and second position.')                
            p_new = g.Sequence(g.CDS(h).indices(end)-(pl+2)-3+1 : g.CDS(h).indices(end)-3); %increase primer lenght by 2        
            display(['CM R1 primer (+2): ',seqrcomplement(p_new)])
            
            if regexp(p_new(1:2), '^[gc]{2}') % new CM R1 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(1:2), '^[gc][at]') % new CM R1 primer has G/C at first position and A/T at second position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(1:2), '^[at][gc]') % new CM R1 primer has A/T at first position and G/C at second position
                %increase lenght of OLD CM R1 primer by 1 (this is same as decreasing length of NEW CM R1 primer by 1)
                p_new = g.Sequence(g.CDS(h).indices(end)-(pl+1)-3+1 : g.CDS(h).indices(end)-3); 
                display(['CM R1 primer (+1): ',seqrcomplement(p_new)])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %new CM R1 primer has A/T at first and second positions
                p = p_new; %this extends the primer ending with e.g. 5'-AA to 5'-AAAA 
                Tm = calculate_Tm_like_PrimerSelect(p); 
            end                    
        end

        
        %CM R1 is a reverse primer --> take reverse complement of primer sequence        
        primer(2).name = [g.CDS(h).gene,'_CM_R1'];             
        primer(2).overhang = seqrcomplement('GGAGCAGGTGCTGGTGCTGG');        
        %check length of CM R1 primer
        %primer cannot be longer than 40 nucleotide since the overhang is already 20 nucleotide and the maximal primer length is 60 nucleotide
        if size(p,2) > 60-length(primer(2).overhang)             
            display(['Need to truncate CM R1 primer to ',num2str(60-length(primer(2).overhang)),' nucleotides (--> total length is 60 nucleotides).'])            
            p = p(end-(60-length(primer(2).overhang))+1:end); % --> truncate to 40 nucleotides           
            Tm = calculate_Tm_like_PrimerSelect(p); %calculate the new Tm
        else
        end                       
        primer(2).seq = [primer(2).overhang,seqrcomplement(p)]; %add 5' appendix to CM R1 primer
        primer(2).Tm = Tm;        
        primer(2).length_binding = length(p);
        primer(2).length_overhang = length(primer(2).overhang);        
        primer(2).length_total = primer(2).length_binding + primer(2).length_overhang; %add appendix to primer length                     
        primer(2).index_3p = regexp(g.Sequence,p); %index of 5'
        primer(2).index_5p = primer(2).index_3p + primer(2).length_binding - 1; %index of 3'                        
        primer(2) = correctPrimerBindingSite(g,h,primer(2));
        clear k i pl p p_new Tm s


        % ---------- Third primer (CM F2 primer) ----------
        display('----- Designing CM F2 primer -----')
        pl = 20;
        p = g.Sequence(g.CDS(h).indices(end)+1 : g.CDS(h).indices(end)+pl);
        Tm = calculate_Tm_like_PrimerSelect(p);

        while Tm < 58.0
            pl = pl +1; %increase primer lenght by 1
            p = g.Sequence(g.CDS(h).indices(end)+1 : g.CDS(h).indices(end)+pl);
            Tm = calculate_Tm_like_PrimerSelect(p);
        end

        %Try to include GC clamp in CM F2 primer if possible; Dirk 20140510    
        display(['Original CM F2 primer: ',p])

        %Case #1: CM F2 primer has already GC clamp
        if regexp(p(end-1:end), '^[gc]{2}')
            display('CM F2 primer has GC clamp.') % --> do nothing

        %Case #2: CM F2 primer has G/C at first position but A/T at second position        
        elseif regexp(p(end-1:end), '^[at][gc]')
            display('CM F2 primer has G/C at first position but A/T at second position.')               
            p_new = g.Sequence(g.CDS(h).indices(end)+1 : g.CDS(h).indices(end)+(pl+1)); %increase primer lenght by 1                
            display(['CM F2 primer (+1): ',p_new])
            
            if regexp(p_new(end), '^[gc]') % new CM F2 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %use old primer seqeunce; CM F2 primer has G/C at first position
            end

        %Case #3: CM F2 primer has A/T at first position but G/C at second last position        
        elseif regexp(p(end-1:end), '^[gc][at]')
            display('CM F2 Primer has A/T at first position but G/C at second position.')                               
            p_new = g.Sequence(g.CDS(h).indices(end)+1 : g.CDS(h).indices(end)+(pl+2)); %increase primer lenght by 2                
            display(['CM F2 primer (+2): ',p_new])
            
            if regexp(p_new(end-1:end), '^[gc]{2}') % new CM F2 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[at][gc]') % new CM F2 primer has G/C at first position and A/T at second first position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[gc][at]') % new CM F2 primer has A/T at first position and G/C at second first position            
                %increase lenght of OLD CM F2 primer by 1 (this is same as decreasing length of NEW CM F2 primer by 1)
                p_new = g.Sequence(g.CDS(h).indices(end)+1 : g.CDS(h).indices(end)+(pl+1));            
                display(['CM F2 primer (+1): ',p_new])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else % new CM F2 primer has A/T at first and second position --> use old primer
            end 

        %Case #4: CM F2 primer has A/T at first and second position
        elseif regexp(p(end-1:end), '^[at]{2}')        
            display('CM F2 primer has A/T at first and second position.')                        
            p_new = g.Sequence(g.CDS(h).indices(end)+1 : g.CDS(h).indices(end)+(pl+2)); %increase primer lenght by 2  
            display(['CM F2 primer (+2): ',p_new])
            
            if regexp(p_new(end-1:end), '^[gc]{2}') % new CM F2 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[at][gc]') % new CM F2 primer has G/C at first position and A/T at second position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[gc][at]') % new CM F2 primer has A/T at first position and G/C at second  position            
                %increase lenght of OLD CM F2 primer by 1 (this is same as decreasing length of NEW CM F2 primer by 1)
                p_new = g.Sequence(g.CDS(h).indices(end)+1 : g.CDS(h).indices(end)+(pl+1)); 
                display(['CM F2 primer (+1): ',p_new])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %new CM F2 primer has A/T at first and second positions
                p = p_new; %this extends the primer ending with e.g. 5'-AA to 5'-AAAA 
                Tm = calculate_Tm_like_PrimerSelect(p);
            end                    
        end

        primer(3).name = [g.CDS(h).gene,'_CM_F2'];                        
        primer(3).overhang = seqrcomplement('TCGATGAATTCGAGCTCGTTTAAAC');
        %check length of CM F2 primer
        %primer cannot be longer than 35 nucleotide since the overhang is already 25 nucleotide and the maximal primer length is 60 nucleotide
        if size(p,2) > 60-length(primer(3).overhang)
            display(['Need to truncate CM F2 primer to ',num2str(60-length(primer(3).overhang)),' nucleotides (--> total length is 60 nucleotides).'])                
            p = p(1:60-length(primer(3).overhang));
            Tm = calculate_Tm_like_PrimerSelect(p);
        else
        end                               
        primer(3).seq = [primer(3).overhang,p]; %add 5' appendix
        primer(3).Tm = Tm;
        primer(3).length_binding = length(p); %add appendix to primer length                
        primer(3).length_overhang = length(primer(3).overhang); %add appendix to primer length                
        primer(3).length_total = primer(3).length_overhang + primer(3).length_binding;
        primer(3).index_5p = regexp(g.Sequence,p);
        primer(3).index_3p = primer(3).index_5p + primer(3).length_binding - 1;                        
        primer(3) = correctPrimerBindingSite(g,h,primer(3));
        clear k i pl p p_new Tm s


        % ---------- Fourth primer (CM R2 primer) ----------
        display('----- Designing CM R2 primer -----')    
        k = regexp(downstream, '[c|g]{2}'); %find all GC clamps
        k = sort(k,'ascend'); %start with the primer closest to stop codon and move further downstream from there    

        for i = 1:size(k,2)
            pl = 20;
            p = downstream(k(i) : k(i)+pl-1);                
            Tm = calculate_Tm_like_PrimerSelect(p);
            
            while Tm < 58.0        
                pl = pl + 1;                        
                p = downstream(k(i) : k(i)+pl-1);
                Tm = calculate_Tm_like_PrimerSelect(p);                
            end    

            s = oligoprop(p);
            if isempty(s.Hairpins)        
                primer(4).name = [g.CDS(h).gene,'_CM_R2'];                
                primer(4).seq = seqrcomplement(p);
                primer(4).Tm = Tm;
                primer(4).overhang = '';
                primer(4).length_binding = length(p);            
                primer(4).length_overhang = length(primer(4).overhang);
                primer(4).length_total = primer(4).length_overhang + primer(4).length_binding;            
                primer(4).index_3p = regexp(g.Sequence,seqrcomplement(primer(4).seq)); %index of 5'
                primer(4).index_5p = primer(4).index_3p + primer(4).length_binding - 1;                                
                break
            else
            end      
        end
        primer(4) = correctPrimerBindingSite(g,h,primer(4));
        clear k i pl p Tm s


        % ---------- Calculate length of PCR product H1_CM ----------        
        H1_length = primer(2).index_5p + 1 + primer(2).length_overhang - primer(1).index_5p;
        display(['Length of PCR product H1_CM is ',num2str(H1_length),' bp.'])
                
         % ---------- Calculate length of PCR product H2_CM ----------        
        H2_length = primer(4).index_5p - primer(3).index_5p + 1 + primer(3).length_overhang;
        display(['Length of PCR product H2_CM is ',num2str(H2_length),' bp.'])
                
        %write PCR product lenght of H1_CM and H2_CM to a text file        
        filename_out3 = [datestr(now,'yyyymmdd_THHMMSS'),'_length_of_H1_CM_and_H2_CM_PCR_products_for_',g.CDS(h).gene,'.txt']; display(filename_out3);
        fid = fopen(filename_out3, 'w','n','US-ASCII'); %write mode
        fprintf(fid, '%s\t%s\t%s\t\n','Gene','Name','Length (bp)');        
        fprintf(fid, '%s\t%s\t%d\t\n',g.CDS(h).gene,['H1_CM_',g.CDS(h).gene], H1_length);
        fprintf(fid, '%s\t%s\t%d\t\n',g.CDS(h).gene,['H2_CM_',g.CDS(h).gene], H2_length);
        fclose (fid);
        
        
        % ---------- Fifth primer: forward primer for colony PCR (CM Fc primer) ----------
        display('----- Designing CM Fc primer -----')           
        k = strfind(upstream,primer(1).seq);
        upstream2 = upstream(1:k-1); %region upstream of the CM F1 primer;
        k = regexp(upstream2, '[c|g]{2}'); %find all GC clamps in the region upstream of the CM F1 primer
        k = sort(k,'descend'); %start with the primer closest to the CM F1 primer and move further upstream from there

        for i = 1:size(k,2) %loop over all potential primers    
            pl = 20; %minimal primer length
            p = upstream2(k(i)-pl+2 : k(i)+1); %size(p,2) = 20   
            Tm = calculate_Tm_like_PrimerSelect(p);
            
            while Tm < 58.0 %Tm of primer needs to be above 58.0 C
                pl = pl +1; %increase primer length by 1       
                p = upstream2(k(i)-pl+2 : k(i)+1);
                Tm = calculate_Tm_like_PrimerSelect(p);
            end

            s = oligoprop(p); %calculate properties of primer
            if isempty(s.Hairpins)                        
                primer(5).name = [g.CDS(h).gene,'_CM_Fc'];
                primer(5).seq = p;
                primer(5).Tm = Tm;
                primer(5).overhang = '';
                primer(5).length_binding = length(p);    
                primer(5).length_overhang = length(primer(5).overhang);    
                primer(5).length_total = primer(5).length_overhang + primer(5).length_binding;                                
                primer(5).index_5p = regexp(g.Sequence,p);                                   
                primer(5).index_3p = primer(5).index_5p + primer(5).length_binding - 1;                                
                break
            else
            end      
        end
        primer(5) = correctPrimerBindingSite(g,h,primer(5));
        clear k i pl p Tm s

                
                

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % --- DNA antisense strand (gene has reverse orientation) ---
    elseif g.CDS(h).indices(end) < g.CDS(h).indices(1)         
        display([g.CDS(h).gene, ' gene has reverse orientation.'])    

        ORF = seqrcomplement(g.Sequence(g.CDS(h).indices(end)+3 : g.CDS(h).indices(1)));    
        stop_codon = seqrcomplement(g.Sequence(g.CDS(h).indices(end) : g.CDS(h).indices(end)+3-1));
        display(['Stop codon: ',upper(stop_codon)])

        try
        %500-bp search window for CM F1 primer; search window starts 280 bp upstream of stop codon
        upstream = seqrcomplement(g.Sequence(g.CDS(h).indices(end)+283 : g.CDS(h).indices(end)+782)); 
        %500-bp search window for CM R2 primer; search window starts 280 bp upstream of stop codon
        downstream = seqrcomplement(g.Sequence(g.CDS(h).indices(end)-780 : g.CDS(h).indices(end)-281));
        catch ME
            display(['Upstream or downstream sequence(s) are not long enough. Skipped primer design for gene ',g.CDS(h).gene,'.'])
            primer = addEmptyPrimerEntries(g.CDS(h).gene);
            return
        end
        
        % ---------- First primer (CM F1 primer) ----------
        display('----- Designing CM F1 primer -----')    
        k = regexp(upstream, '[c|g]{2}'); %find all GC clamps
        k = sort(k,'descend'); %start with the primer closest to stop codon and move further upstream from there
        
        for i = 1:size(k,2) %loop over all potential primers    
            pl = 20; %minimal primer length
            p = upstream(k(i)-pl+2 : k(i)+1); %size(p,2) = 20       
            Tm = calculate_Tm_like_PrimerSelect(p);
        
            while Tm < 58.0 %Tm of primer needs to be above 58.0 C
                pl = pl + 1; %increase primer length by 1       
                p = upstream(k(i)-pl+2 : k(i)+1);
                Tm = calculate_Tm_like_PrimerSelect(p);
            end

            s = oligoprop(p); %calculate properties of primer
            if isempty(s.Hairpins)                        
                primer(1).name = [g.CDS(h).gene,'_CM_F1'];
                primer(1).seq = p;
                primer(1).Tm = Tm;
                primer(1).overhang = '';
                primer(1).length_binding = length(p);                        
                primer(1).length_overhang = 0;
                primer(1).length_total = primer(1).length_overhang + primer(1).length_binding;
                primer(1).index_3p = regexp(g.Sequence,seqrcomplement(primer(1).seq)); 
                primer(1).index_5p = primer(1).index_3p + primer(1).length_binding - 1;
                orderfields(primer(1),{'name','seq','Tm','overhang','length_binding','length_overhang','length_total','index_5p','index_3p'});                
                break
            else
            end      
        end
        primer(1) = correctPrimerBindingSite(g,h,primer(1));
        clear k i pl p Tm s
        %Note, not necessary to check whether CM F1 primer is larger than 60 nt since primer with 35x As or 35x Ts has a Tm of 58.1436 C  


        % ---------- Second primer (CM R1 primer) ----------
        display('----- Designing CM R1 primer -----')                    
        pl = 20; %minimal primer length    
        p = g.Sequence(g.CDS(h).indices(end)+3 : g.CDS(h).indices(end)+3-1+pl);
        Tm = calculate_Tm_like_PrimerSelect(p);
        
        while Tm < 58.0
            pl = pl + 1; %increase primer lenght by 1        
            p = g.Sequence(g.CDS(h).indices(end)+3 : g.CDS(h).indices(end)+3-1+pl);
            Tm = calculate_Tm_like_PrimerSelect(p);
        end

        %Try to include GC clamp if possible
        display(['Original CM R1 primer: ',p])
        %no need to take the reverse complement here since the gene is already in reverse orientation

        %Case #1: CM R1 primer has GC clamp
        if regexp(p(end-1:end), '^[gc]{2}')
            display('CM R1 primer has GC clamp.') % --> do nothing

        %Case #2: CM R1 primer has G/C at first position but A/T at second position    
        elseif regexp(p(end-1:end), '^[at][gc]')
            display('CM R1 primer has G/C at first position but A/T at second position.')                        
            p_new = g.Sequence(g.CDS(h).indices(end)+3 : g.CDS(h).indices(end)+3-1+(pl+1));%increase primer lenght by 1
            display(['CM R1 primer (+1): ',p_new])
            
            if regexp(p_new(end), '^[gc]') % new CM R1 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %use old seqeunce; CM R1 primer has G/C at first position
            end

        %Case #3: CM R1 primer has A/T at first position but G/C at second last position    
        elseif regexp(p(end-1:end), '^[gc][at]')
            display('CM R1 primer has A/T at first position but G/C at second position.')                           
            p_new = g.Sequence(g.CDS(h).indices(end)+3 : g.CDS(h).indices(end)+3-1+(pl+2)); %increase primer lenght by 2
            display(['CM R1 primer (+2): ',p_new])
            
            if regexp(p_new(end-1:end), '^[gc]{2}') % new CM R1 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[at][gc]') % new CM R1 primer has G/C at first position and A/T at second first position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[gc][at]') % new CM R1 primer has A/T at first position and G/C at second first position
                %increase lenght of OLD CM R1 primer by 1 (this is same as decreasing length of NEW CM R1 primer by 1)
                p_new = g.Sequence(g.CDS(h).indices(end)+3 : g.CDS(h).indices(end)+3-1+(pl+1)); 
                display(['CM R1 primer (+1): ',p_new])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else % new CM R1 primer has A/T at first and second position --> use old CM R1 primer
            end 

        %Case #4: CM R1 primer has A/T at first and second position
        elseif regexp(p(end-1:end), '^[at]{2}')        
            display('CM R1 primer has A/T at first and second position.')                        
            p_new = g.Sequence(g.CDS(h).indices(end)+3 : g.CDS(h).indices(end)+3-1+(pl+2)); %increase primer lenght by 2
            display(['CM R1 primer (+2): ',p_new])
            
            if regexp(p_new(end-1:end), '^[gc]{2}') % new CM R1 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[at][gc]') % new CM R1 primer has G/C at first position and A/T at second position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[gc][at]') % new CM R1 primer has A/T at first position and G/C at second position
                %increase lenght of OLD CM R1 primer by 1 (this is same as decreasing length of NEW CM R1 primer by 1)
                p_new = g.Sequence(g.CDS(h).indices(end)+3 : g.CDS(h).indices(end)+3-1+(pl+1));
                display(['CM R1 primer (+1): ',p_new])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %new CM R1 primer has A/T at first and second positions
                p = p_new; %this extends the primer ending with e.g. 5'-AA to 5'-AAAA 
                Tm = calculate_Tm_like_PrimerSelect(p); 
            end                    
        end

        primer(2).name = [g.CDS(h).gene,'_CM_R1'];
        %no need to take the reverse complement here since the gene is already in reverse orientation        
        primer(2).overhang = seqrcomplement('GGAGCAGGTGCTGGTGCTGG');
        %check length of CM R1 primer
        %primer cannot be longer than 40 nucleotide since the overhang is already 20 nucleotide and the maximal primer length is 60 nucleotide
        if size(p,2) > 60-length(primer(2).overhang) 
            display(['Need to truncate CM R1 primer to ',num2str(60-length(primer(2).overhang)),' nucleotides (--> total length is 60 nucleotides).'])
            p = p(1:60-length(primer(2).overhang)); % --> truncate to 40 nucleotides           
            Tm = calculate_Tm_like_PrimerSelect(p);
        else
        end                    
        primer(2).seq = [primer(2).overhang,p]; %add 5' appendix to CM R1 primer
        primer(2).Tm = Tm;        
        primer(2).length_binding = length(p);
        primer(2).length_overhang = length(primer(2).overhang);        
        primer(2).length_total = primer(2).length_binding + primer(2).length_overhang; %add appendix to primer length             
        primer(2).index_5p = regexp(g.Sequence,p); %index of 5'
        primer(2).index_3p = primer(2).index_5p + primer(2).length_binding - 1; %index of 3'
        primer(2) = correctPrimerBindingSite(g,h,primer(2));
        clear k i pl p p_new Tm s
                               
        
        % ---------- Third primer (CM F2 primer) ----------
        display('----- Designing CM F2 primer -----')
        pl = 20;
        p = seqrcomplement(g.Sequence(g.CDS(h).indices(end)-pl : g.CDS(h).indices(end)-1)); %this is now a reverse primer
        Tm = calculate_Tm_like_PrimerSelect(p);
        
        while Tm < 58.0
            pl = pl + 1; %increase primer lenght by 1
            p = seqrcomplement(g.Sequence(g.CDS(h).indices(end)-pl : g.CDS(h).indices(end)-1));
            Tm = calculate_Tm_like_PrimerSelect(p);
        end

        %Try to include GC clamp in CM F2 primer if possible; Dirk 20140513    
        display(['Original CM F2 primer: ',p])

        %Case #1: CM F2 primer has already GC clamp
        if regexp(p(end-1:end), '^[gc]{2}')
            display('CM F2 primer has GC clamp.') % --> do nothing

        %Case #2: CM F2 primer has G/C at first position but A/T at second position            
        elseif regexp(p(end-1:end), '^[at][gc]')
            display('CM F2 primer has G/C at first position but A/T at second position.')                   
            p_new = seqrcomplement(g.Sequence(g.CDS(h).indices(end)-(pl+1) : g.CDS(h).indices(end)-1)); %increase primer lenght by 1
            display(['CM F2 primer (+1): ',p_new])
            
            if regexp(p_new(end), '^[gc]') % new primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %use old primer seqeunce; CM F2 primer has G/C at first position
            end

        %Case #3: CM F2 primer has A/T at first position but G/C at second last position            
        elseif regexp(p(end-1:end), '^[gc][at]')
            display('CM F2 primer has A/T at first position but G/C at second position.')                               
            p_new = seqrcomplement(g.Sequence(g.CDS(h).indices(end)-(pl+2) : g.CDS(h).indices(end)-1)); %increase primer lenght by 2
            display(['CM F2 primer (+2): ',p_new])
            
            if regexp(p_new(end-1:end), '^[gc]{2}') % new primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[at][gc]') % new primer has G/C at first position and A/T at second first position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[gc][at]') % new primer has A/T at first position and G/C at second first position            
                %increase lenght of OLD CM F2 primer by 1 (this is same as decreasing length of NEW CM F2 primer by 1)
                p_new = seqrcomplement(g.Sequence(g.CDS(h).indices(end)-(pl+1) : g.CDS(h).indices(end)-1)); 
                display(['CM F2 primer (+1): ',p_new])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else % new CM F2 primer has A/T at first and second position --> use old CM F2 primer
            end         

        %Case #4: CM F2 primer has A/T at first and second position
        elseif regexp(p(end-1:end), '^[at]{2}')        
            display('CM F2 primer has A/T at first and second position.')                                
            p_new = seqrcomplement(g.Sequence(g.CDS(h).indices(end)-(pl+2) : g.CDS(h).indices(end)-1)); %increase primer lenght by 2
            display(['CM F2 primer (+2): ',p_new])
            
            if regexp(p_new(end-1:end), '^[gc]{2}') % new primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[at][gc]') % new primer has G/C at first position and A/T at second position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[gc][at]') % new primer has A/T at first position and G/C at second position            
            %increase lenght of OLD CM F2 primer by 1 (this is same as decreasing length of NEW CM F2 primer by 1) 
            p_new = seqrcomplement(g.Sequence(g.CDS(h).indices(end)-(pl+1) : g.CDS(h).indices(end)-1));            
            display(['CM F2 primer (+1): ',p_new])
            p = p_new;
            Tm = calculate_Tm_like_PrimerSelect(p);
            else %new CM F2 primer has A/T at first and second positions
                p = p_new; %this extends the primer ending with e.g. 5'-AA to 5'-AAAA 
                Tm = calculate_Tm_like_PrimerSelect(p);
            end                    
        end

        primer(3).name = [g.CDS(h).gene,'_CM_F2'];        
        primer(3).overhang = seqrcomplement('TCGATGAATTCGAGCTCGTTTAAAC');
        %check length of CM F2 primer
        %primer cannot be longer than 35 nucleotide since the overhang is already 25 nucleotide and the maximal primer length is 60 nucleotide
        if size(p,2) > 60-length(primer(3).overhang)
            display(['Need to truncate CM F2 primer to ',num2str(60-length(primer(3).overhang)),' nucleotides (--> total length is 60 nucleotides).'])                 
            p = p(1:60-length(primer(3).overhang)); % --> truncate to 35 nucleotides           
            Tm = calculate_Tm_like_PrimerSelect(p);
        else
        end                     
        primer(3).seq = [primer(3).overhang,p]; %add 5' appendix
        primer(3).Tm = Tm;
        primer(3).length_binding = length(p); %add appendix to primer length                
        primer(3).length_overhang = length(primer(3).overhang); %add appendix to primer length                
        primer(3).length_total = primer(3).length_overhang + primer(3).length_binding;
        primer(3).index_3p = regexp(g.Sequence,seqrcomplement(p));
        primer(3).index_5p = primer(3).index_3p + primer(3).length_binding - 1;
        primer(3) = correctPrimerBindingSite(g,h,primer(3));
        clear k i pl p p_new Tm s


        % ---------- Fourth primer (CM R2 primer) ----------
        display('----- Designing CM R2 primer -----')    
        k = regexp(downstream, '[c|g]{2}'); %find all GC clamps
        k = sort(k,'ascend'); %start with the primer closest to stop codon and move further downstream from there    
        
        for i = 1:size(k,2)   
            pl = 20;
            p = downstream(k(i) : k(i)+pl-1);    
            Tm = calculate_Tm_like_PrimerSelect(p); 
            
            while Tm < 58.0
                pl = pl +1;                   
                p = downstream(k(i) : k(i)+pl-1);
                Tm = calculate_Tm_like_PrimerSelect(p);                        
            end    

            s = oligoprop(p);
            if isempty(s.Hairpins)                        
                primer(4).name = [g.CDS(h).gene,'_CM_R2'];
                primer(4).seq = seqrcomplement(p);
                primer(4).Tm = Tm;
                primer(4).overhang = '';
                primer(4).length_binding = length(p);            
                primer(4).length_overhang = length(primer(4).overhang);
                primer(4).length_total = primer(4).length_overhang + primer(4).length_binding;            
                primer(4).index_5p = regexp(g.Sequence,primer(4).seq); %index of 5'
                primer(4).index_3p = primer(4).index_5p + primer(4).length_binding - 1;
                break
            else
            end      
        end
        primer(4) = correctPrimerBindingSite(g,h,primer(4));
        clear k i pl p Tm s

        
         % ---------- Calculate length of PCR product H1_CM ----------        
        H1_length = primer(1).index_5p - primer(2).index_5p + 1 + primer(2).length_overhang;
        display(['Length of PCR product H1_CM is ',num2str(H1_length),' bp.'])
                
         % ---------- Calculate length of PCR product H2_CM ----------        
        H2_length = primer(3).index_5p - primer(4).index_5p + 1 + primer(3).length_overhang;
        display(['Length of PCR product H2_CM is ',num2str(H2_length),' bp.'])
                        
        %write PCR product lenght of H1_CM and H2_CM to a text file        
        filename_out3 = [datestr(now,'yyyymmdd_THHMMSS'),'_length_of_H1_CM_and_H2_CM_PCR_products_for_',g.CDS(h).gene,'.txt']; display(filename_out3);
        fid = fopen(filename_out3, 'w','n','US-ASCII'); %write mode
        fprintf(fid, '%s\t%s\t%s\t\n','Gene','Name','Length (bp)');        
        fprintf(fid, '%s\t%s\t%d\t\n',g.CDS(h).gene,['H1_CM_',g.CDS(h).gene], H1_length);
        fprintf(fid, '%s\t%s\t%d\t\n',g.CDS(h).gene,['H2_CM_',g.CDS(h).gene], H2_length);
        fclose (fid);
        

        % ---------- Fifth primer: forward primer for colony PCR (CM Fc primer) ----------
        display('----- Designing CM Fc primer -----')           
        k = strfind(upstream,primer(1).seq);
        upstream2 = upstream(1:k-1); %region upstream of the CM F1 primer; k(i) comes from the calculation of the CM F1 primer        
        k = regexp(upstream2, '[c|g]{2}'); %find all GC clamps in the region upstream of the CM F1 primer
        k = sort(k,'descend'); %start with the primer closest to the CM F1 primer and move further upstream from there
        
        for i = 1:size(k,2) %loop over all potential primers    
            pl = 20; %minimal primer length
            p = upstream2(k(i)-pl+2 : k(i)+1); %size(p,2) = 20   
            Tm = calculate_Tm_like_PrimerSelect(p);
            
            while Tm < 58.0 %Tm of primer needs to be above 58.0 C
                pl = pl + 1; %increase primer length by 1       
                p = upstream2(k(i)-pl+2 : k(i)+1);
                Tm = calculate_Tm_like_PrimerSelect(p);
            end

            s = oligoprop(p); %calculate properties of primer
            if isempty(s.Hairpins)                        
                primer(5).name = [g.CDS(h).gene,'_CM_Fc'];
                primer(5).seq = p;
                primer(5).Tm = Tm;
                primer(5).overhang = '';
                primer(5).length_binding = length(p);    
                primer(5).length_overhang = length(primer(5).overhang);    
                primer(5).length_total = primer(5).length_overhang + primer(5).length_binding;                
                primer(5).index_3p = regexp(g.Sequence,seqrcomplement(p));
                primer(5).index_5p = primer(5).index_3p + primer(5).length_binding -1;                                
                break
            else
            end      
        end
        primer(5) = correctPrimerBindingSite(g,h,primer(5));
        clear k i pl p Tm s
                      
    else
    end
    

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------ METHOD #2 ----- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'C-terminal'
            
    %First determine whether gene is on the upper (sense) or lower (anti-sense) strand
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- DNA sense strand (gene has normal orientation) ---
    if g.CDS(h).indices(end) > g.CDS(h).indices(1)         
        display([g.CDS(h).gene, ' gene has normal orientation.'])    
        ORF = g.Sequence(g.CDS(h).indices(1) : g.CDS(h).indices(end)-3); %exclude stop codon
        stop_codon = g.Sequence(g.CDS(h).indices(end)-2 : g.CDS(h).indices(end)); display(['Stop codon: ',upper(stop_codon)])

        try
        %500-bp search window for C F1 primer; search window starts 280 bp upstream of stop codon
        upstream = g.Sequence(g.CDS(h).indices(end)-782 : g.CDS(h).indices(end)-283);
        %500-bp search window for C R2 primer; search window starts 280 bp upstream of stop codon        
        downstream = g.Sequence(g.CDS(h).indices(end)+281-3 : g.CDS(h).indices(end)+780-3);                
        catch ME
            display(['Upstream or downsream sequences are not long enough. Skipped primer design for gene ',g.CDS(h).gene,'.'])
            primer = addEmptyPrimerEntries(g.CDS(h).gene);
            return
        end
        
        % ---------- First primer (C F1 primer) ----------
        display('----- Designing C F1 primer -----')    
        k = regexp(upstream, '[c|g]{2}'); %find all GC clamps
        k = sort(k,'descend'); %start with the primer closest to stop codon (i.e. the insertion site) and move further upstream from there

        for i = 1:size(k,2) %loop over all potential primers
            pl = 20; %minimal primer length
            p = upstream(k(i)-pl+2 : k(i)+1); %size(p,2) = 20   
            Tm = calculate_Tm_like_PrimerSelect(p);
            
            while Tm < 58.0 %Tm of primer needs to be above 58.0 C
                pl = pl +1; %increase primer length by 1       
                p = upstream(k(i)-pl+2 : k(i)+1);
                Tm = calculate_Tm_like_PrimerSelect(p);
            end

            s = oligoprop(p); %calculate properties of primer
            if isempty(s.Hairpins)                        
                primer(1).name = [g.CDS(h).gene,'_C_F1'];
                primer(1).seq = p;
                primer(1).Tm = Tm;                                
                primer(1).overhang = '';
                primer(1).length_binding = length(p);                        
                primer(1).length_overhang = 0;
                primer(1).length_total = primer(1).length_overhang + primer(1).length_binding;                
                primer(1).index_5p = regexp(g.Sequence,primer(1).seq);   
                primer(1).index_3p = primer(1).index_5p + primer(1).length_binding - 1;      
                orderfields(primer(1),{'name','seq','Tm','overhang','length_binding','length_overhang','length_total','index_5p','index_3p'});
                break
            else
            end      
        end               
        primer(1) = correctPrimerBindingSite(g,h,primer(1));        
        clear k i pl p Tm s
        %Note, not necessary to check whether C F1 primer is larger than 60 nt since a primer with 35x As or 35x Ts has a Tm of 58.1436 C              


        % ---------- Second primer (C R1 primer) ----------
        display('----- Designing C R1 primer -----')    
        pl = 20; %minimal primer length
        p = g.Sequence(g.CDS(h).indices(end)-pl-3+1 : g.CDS(h).indices(end)-3);   
        Tm = calculate_Tm_like_PrimerSelect(p);
        
        while Tm < 58.0
            pl = pl +1; %increase primer lenght by 1
            p = g.Sequence(g.CDS(h).indices(end)-pl-3+1 : g.CDS(h).indices(end)-3);
            Tm = calculate_Tm_like_PrimerSelect(p);
        end

        %Try to include GC clamp if possible
        display(['Original C R1 primer: ',seqrcomplement(p)])

        %Case #1: C R1 primer has GC clamp
        if regexp(p(1:2), '^[gc]{2}')
            display('C R1 primer has GC clamp.') % --> do nothing since primer has GC clamp

        %Case #2: C R1 primer has G/C at first position but A/T at second position    
        elseif regexp(p(1:2), '^[gc][at]')
            display('C R1 primer has G/C at first position but A/T at second position.')                
            p_new = g.Sequence(g.CDS(h).indices(end)-(pl+1)-3+1 : g.CDS(h).indices(end)-3); %increase primer lenght by 1        
            display(['C R1 primer (+1): ',seqrcomplement(p_new)])
            
            if regexp(p_new(1), '^[gc]') % new C R1 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %use old seqeunce; C R1 primer has G/C at first position
            end

        %Case #3: C R1 primer has A/T at first position but G/C at second last position    
        elseif regexp(p(1:2), '^[at][gc]')
            display('C R1 primer has A/T at first position but G/C at second position.')               
            p_new = g.Sequence(g.CDS(h).indices(end)-(pl+2)-3+1 : g.CDS(h).indices(end)-3); %increase primer lenght by 2        
            display(['C R1 primer (+2): ',seqrcomplement(p_new)])
            
            if regexp(p_new(1:2), '^[gc]{2}') % new C R1 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(1:2), '^[gc][at]') % new C R1 primer has G/C at first position and A/T at second first position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(1:2), '^[at][gc]') % new C R1 primer has A/T at first position and G/C at second first position
                %increase lenght of OLD C R1 primer by 1 (this is same as decreasing length of NEW C R1 primer by 1)
                p_new = g.Sequence(g.CDS(h).indices(end)-(pl+1)-3+1 : g.CDS(h).indices(end)-3);                     
                display(['C R1 primer (+1): ',seqrcomplement(p_new)])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else % new C R1 primer has A/T at first and second position --> use old C R1 primer
            end 

        %Case #4: C R1 primer has A/T at first and second position
        elseif regexp(p(1:2), '^[at]{2}')        
            display('C R1 primer has A/T at first and second position.')                
            p_new = g.Sequence(g.CDS(h).indices(end)-(pl+2)-3+1 : g.CDS(h).indices(end)-3); %increase primer lenght by 2        
            display(['C R1 primer (+2): ',seqrcomplement(p_new)])
            
            if regexp(p_new(1:2), '^[gc]{2}') % new C R1 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(1:2), '^[gc][at]') % new C R1 primer has G/C at first position and A/T at second position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(1:2), '^[at][gc]') % new C R1 primer has A/T at first position and G/C at second position
                %increase lenght of OLD C R1 primer by 1 (this is same as decreasing length of NEW C R1 primer by 1)
                p_new = g.Sequence(g.CDS(h).indices(end)-(pl+1)-3+1 : g.CDS(h).indices(end)-3); 
                display(['C R1 primer (+1): ',seqrcomplement(p_new)])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %new C R1 primer has A/T at first and second positions
                p = p_new; %this extends the primer ending with e.g. 5'-AA to 5'-AAAA 
                Tm = calculate_Tm_like_PrimerSelect(p); 
            end                    
        end

        %C R1 is a reverse primer --> take reverse complement of primer sequence        
        primer(2).name = [g.CDS(h).gene,'_C_R1'];                
        primer(2).overhang = seqrcomplement('GGAGCAGGTGCTGGTGCTGG');
        %check length of C R1 primer
        %primer cannot be longer than 40 nucleotide since the overhang is already 20 nucleotide and the maximal primer length is 60 nucleotide
        if size(p,2) > 60-length(primer(2).overhang) 
            display(['Need to truncate C R1 primer to ',num2str(60-length(primer(2).overhang)),' nucleotides (--> total length is 60 nucleotides).'])                 
            p = p(end-(60-length(primer(2).overhang))+1:end); % --> truncate to 40 nucleotides           
            Tm = calculate_Tm_like_PrimerSelect(p); %calculate the new Tm
        else
        end            
        primer(2).seq = [primer(2).overhang,seqrcomplement(p)]; %add 5' appendix to C R1 primer
        primer(2).Tm = Tm;        
        primer(2).length_binding = length(p);
        primer(2).length_overhang = length(primer(2).overhang);        
        primer(2).length_total = primer(2).length_binding + primer(2).length_overhang; %add appendix to primer length                     
        primer(2).index_3p = regexp(g.Sequence,p); %index of 5'
        primer(2).index_5p = primer(2).index_3p + primer(2).length_binding - 1; %index of 3'                        
        primer(2) = correctPrimerBindingSite(g,h,primer(2));
        clear k i pl p p_new Tm s


        % ---------- Third primer (C F2 primer) ----------
        display('----- Designing C F2 primer -----')
        pl = 20;        
        p = g.Sequence(g.CDS(h).indices(end)+1-3 : g.CDS(h).indices(end)+pl-3); 
        Tm = calculate_Tm_like_PrimerSelect(p);

        while Tm < 58.0
            pl = pl +1; %increase primer lenght by 1            
            p = g.Sequence(g.CDS(h).indices(end)+1-3 : g.CDS(h).indices(end)+pl-3);
            Tm = calculate_Tm_like_PrimerSelect(p);
        end

        %Try to include GC clamp in C F2 primer if possible
        display(['Original C F2 primer: ',p])

        %Case #1: C F2 primer has already GC clamp
        if regexp(p(end-1:end), '^[gc]{2}')
            display('C F2 primer has GC clamp.') % --> do nothing

        %Case #2: C F2 primer has G/C at first position but A/T at second position        
        elseif regexp(p(end-1:end), '^[at][gc]')
            display('C F2 primer has G/C at first position but A/T at second position.')                           
            p_new = g.Sequence(g.CDS(h).indices(end)+1-3 : g.CDS(h).indices(end)+(pl+1)-3); %increase primer lenght by 1                
            display(['C F2 primer (+1): ',p_new])
            
            if regexp(p_new(end), '^[gc]') % new C F2 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %use old primer seqeunce; C F2 primer has G/C at first position
            end

        %Case #3: C F2 primer has A/T at first position but G/C at second last position        
        elseif regexp(p(end-1:end), '^[gc][at]')
            display('C F2 Primer has A/T at first position but G/C at second position.')                                           
            p_new = g.Sequence(g.CDS(h).indices(end)+1-3 : g.CDS(h).indices(end)+(pl+2)-3); %increase primer lenght by 2                        
            display(['C F2 primer (+2): ',p_new])
            
            if regexp(p_new(end-1:end), '^[gc]{2}') % new C F2 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[at][gc]') % new C F2 primer has G/C at first position and A/T at second first position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[gc][at]') % new C F2 primer has A/T at first position and G/C at second first position            
                %increase lenght of OLD C F2 primer by 1 (this is same as decreasing length of NEW C F2 primer by 1)                
                p_new = g.Sequence(g.CDS(h).indices(end)+1-3 : g.CDS(h).indices(end)+(pl+1)-3);      
                display(['C F2 primer (+1): ',p_new])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else % new C F2 primer has A/T at first and second position --> use old primer
            end 

        %Case #4: C F2 primer has A/T at first and second position
        elseif regexp(p(end-1:end), '^[at]{2}')        
            display('C F2 primer has A/T at first and second position.')                                    
            p_new = g.Sequence(g.CDS(h).indices(end)+1-3 : g.CDS(h).indices(end)+(pl+2)-3); %increase primer lenght by 2  
            display(['C F2 primer (+2): ',p_new])
            
            if regexp(p_new(end-1:end), '^[gc]{2}') % new C F2 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[at][gc]') % new C F2 primer has G/C at first position and A/T at second position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[gc][at]') % new C F2 primer has A/T at first position and G/C at second  position            
                %increase lenght of OLD C F2 primer by 1 (this is same as decreasing length of NEW C F2 primer by 1)                
                p_new = g.Sequence(g.CDS(h).indices(end)+1-3 : g.CDS(h).indices(end)+(pl+1)-3); 
                display(['C F2 primer (+1): ',p_new])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %new C F2 primer has A/T at first and second positions
                p = p_new; %this extends the primer ending with e.g. 5'-AA to 5'-AAAA 
                Tm = calculate_Tm_like_PrimerSelect(p);
            end                    
        end

        primer(3).name = [g.CDS(h).gene,'_C_F2'];                                
        primer(3).overhang = seqrcomplement('CTTGTACAATTCGTCCATACCCATAACG'); %e.g. using pDML219
        %check length of C F2 primer
        %primer cannot be longer than 32 nucleotide since the overhang is already 28 nucleotide and the maximal primer length is 60 nucleotide
        if size(p,2) > 60-length(primer(3).overhang)
            display(['Need to truncate C F2 primer to ',num2str(60-length(primer(3).overhang)),' nucleotides (--> total length is 60 nucleotides).'])             
            p = p(1:60-length(primer(3).overhang)); % --> truncate to 32 nucleotides           
            Tm = calculate_Tm_like_PrimerSelect(p);
        else
        end                               
        primer(3).seq = [primer(3).overhang,p]; %add 5' appendix
        primer(3).Tm = Tm;
        primer(3).length_binding = length(p); %add appendix to primer length                
        primer(3).length_overhang = length(primer(3).overhang); %add appendix to primer length                
        primer(3).length_total = primer(3).length_overhang + primer(3).length_binding;
        primer(3).index_5p = regexp(g.Sequence,p);
        primer(3).index_3p = primer(3).index_5p + primer(3).length_binding - 1;                        
        primer(3) = correctPrimerBindingSite(g,h,primer(3));
        clear k i pl p p_new Tm s


        % ---------- Fourth primer (C R2 primer) ----------
        display('----- Designing C R2 primer -----')    
        k = regexp(downstream, '[c|g]{2}'); %find all GC clamps
        k = sort(k,'ascend'); %start with the primer closest to stop codon and move further downstream from there    

        for i = 1:size(k,2)
            pl = 20;
            p = downstream(k(i) : k(i)+pl-1);                
            Tm = calculate_Tm_like_PrimerSelect(p); 
            
            while Tm < 58.0        
                pl = pl + 1;                        
                p = downstream(k(i) : k(i)+pl-1);
                Tm = calculate_Tm_like_PrimerSelect(p);                
            end    

            s = oligoprop(p);
            if isempty(s.Hairpins)        
                primer(4).name = [g.CDS(h).gene,'_C_R2'];                
                primer(4).seq = seqrcomplement(p);
                primer(4).Tm = Tm;
                primer(4).overhang = '';
                primer(4).length_binding = length(p);            
                primer(4).length_overhang = length(primer(4).overhang);
                primer(4).length_total = primer(4).length_overhang + primer(4).length_binding;            
                primer(4).index_3p = regexp(g.Sequence,seqrcomplement(primer(4).seq)); %index of 5'
                primer(4).index_5p = primer(4).index_3p + primer(4).length_binding - 1;                                
                break
            else
            end      
        end
        primer(4) = correctPrimerBindingSite(g,h,primer(4));
        clear k i pl p Tm s


        % ---------- Calculate length of PCR product H1_C ----------        
        H1_length = primer(2).index_5p + 1 + primer(2).length_overhang - primer(1).index_5p;
        display(['Length of PCR product H1_C is ',num2str(H1_length),' bp.'])
                
         % ---------- Calculate length of PCR product H2_C ----------        
        H2_length = primer(4).index_5p - primer(3).index_5p + 1 + primer(3).length_overhang;
        display(['Length of PCR product H2_C is ',num2str(H2_length),' bp.'])
                
        %write PCR product lenght of H1_C and H2_C to a text file        
        filename_out3 = [datestr(now,'yyyymmdd_THHMMSS'),'_length_of_H1_C_and_H2_C_PCR_products_for_',g.CDS(h).gene,'.txt']; display(filename_out3);
        fid = fopen(filename_out3, 'w','n','US-ASCII'); %write mode
        fprintf(fid, '%s\t%s\t%s\t\n','Gene','Name','Length (bp)');        
        fprintf(fid, '%s\t%s\t%d\t\n',g.CDS(h).gene,['H1_C_',g.CDS(h).gene], H1_length);
        fprintf(fid, '%s\t%s\t%d\t\n',g.CDS(h).gene,['H2_C_',g.CDS(h).gene], H2_length);
        fclose (fid);
        
        
        % ---------- Fifth primer: forward primer for colony PCR (C Fc primer) ----------
        display('----- Designing C Fc primer -----')           
        k = strfind(upstream,primer(1).seq);
        upstream2 = upstream(1:k-1); %region upstream of the C F1 primer;
        k = regexp(upstream2, '[c|g]{2}'); %find all GC clamps in the region upstream of the C F1 primer
        k = sort(k,'descend'); %start with the primer closest to the C F1 primer and move further upstream from there
        
        for i = 1:size(k,2) %loop over all potential primers    
            pl = 20; %minimal primer length
            p = upstream2(k(i)-pl+2 : k(i)+1); %size(p,2) = 20   
            Tm = calculate_Tm_like_PrimerSelect(p);
        
            while Tm < 58.0 %Tm of primer needs to be above 58.0 C
                pl = pl+1; %increase primer length by 1       
                p = upstream2(k(i)-pl+2 : k(i)+1);
                Tm = calculate_Tm_like_PrimerSelect(p);
            end

            s = oligoprop(p); %calculate properties of primer
            if isempty(s.Hairpins)                        
                primer(5).name = [g.CDS(h).gene,'_C_Fc'];
                primer(5).seq = p;
                primer(5).Tm = Tm;
                primer(5).overhang = '';
                primer(5).length_binding = length(p);    
                primer(5).length_overhang = length(primer(5).overhang);    
                primer(5).length_total = primer(5).length_overhang + primer(5).length_binding;                                
                primer(5).index_5p = regexp(g.Sequence,p);                                   
                primer(5).index_3p = primer(5).index_5p + primer(5).length_binding - 1;                                
                break
            else
            end      
        end
        primer(5) = correctPrimerBindingSite(g,h,primer(5));
        clear k i pl p Tm s
                        
        


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % --- DNA antisense strand (gene has reverse orientation) ---
    elseif g.CDS(h).indices(end) < g.CDS(h).indices(1)         
        display([g.CDS(h).gene, ' gene has reverse orientation.'])    

        ORF = seqrcomplement(g.Sequence(g.CDS(h).indices(end)+3 : g.CDS(h).indices(1)));    
        stop_codon = seqrcomplement(g.Sequence(g.CDS(h).indices(end) : g.CDS(h).indices(end)+3-1));
        display(['Stop codon: ',upper(stop_codon)])

        try
        %500-bp search window for CM F1 primer; search window starts 280 bp upstream of stop codon
        upstream = seqrcomplement(g.Sequence(g.CDS(h).indices(end)+283 : g.CDS(h).indices(end)+782));
        %500-bp search window for CM R2 primer; search window starts 280 bp upstream of stop codon        
        downstream = seqrcomplement(g.Sequence(g.CDS(h).indices(end)-780+3 : g.CDS(h).indices(end)-281+3)); 
        catch ME
            display(['Upstream or downstream sequence(s) are not long enough. Skipped primer design for gene ',g.CDS(h).gene,'.'])
            primer = addEmptyPrimerEntries(g.CDS(h).gene);
            return
        end

        % ---------- First primer (C F1 primer) ----------
        display('----- Designing C F1 primer -----')    
        k = regexp(upstream, '[c|g]{2}'); %find all GC clamps
        k = sort(k,'descend'); %start with the primer closest to stop codon and move further upstream from there
        
        for i = 1:size(k,2) %loop over all potential primers    
            pl = 20; %minimal primer length
            p = upstream(k(i)-pl+2 : k(i)+1); %size(p,2) = 20       
            Tm = calculate_Tm_like_PrimerSelect(p);
        
            while Tm < 58.0 %Tm of primer needs to be above 58.0 C
                pl = pl + 1; %increase primer length by 1       
                p = upstream(k(i)-pl+2 : k(i)+1);
                Tm = calculate_Tm_like_PrimerSelect(p);
            end

            s = oligoprop(p); %calculate properties of primer
            if isempty(s.Hairpins)                        
                primer(1).name = [g.CDS(h).gene,'_C_F1'];
                primer(1).seq = p;
                primer(1).Tm = Tm;
                primer(1).overhang = '';
                primer(1).length_binding = length(p);                        
                primer(1).length_overhang = 0;
                primer(1).length_total = primer(1).length_overhang + primer(1).length_binding;
                primer(1).index_3p = regexp(g.Sequence,seqrcomplement(primer(1).seq)); 
                primer(1).index_5p = primer(1).index_3p + primer(1).length_binding - 1;
                orderfields(primer(1),{'name','seq','Tm','overhang','length_binding','length_overhang','length_total','index_5p','index_3p'});                
                break
            else
            end      
        end
        primer(1) = correctPrimerBindingSite(g,h,primer(1));
        clear k i pl p Tm s
        %Note, not necessary to check whether C F1 primer is larger than 60 nt since primer with 35x As or 35x Ts has a Tm of 58.1436 C  


        % ---------- Second primer (C R1 primer) ----------
        display('----- Designing C R1 primer -----')                    
        pl = 20; %minimal primer length    
        p = g.Sequence(g.CDS(h).indices(end)+3 : g.CDS(h).indices(end)+3-1+pl);
        Tm = calculate_Tm_like_PrimerSelect(p);
        
        while Tm < 58.0
            pl = pl + 1; %increase primer lenght by 1        
            p = g.Sequence(g.CDS(h).indices(end)+3 : g.CDS(h).indices(end)+3-1+pl);
            Tm = calculate_Tm_like_PrimerSelect(p);
        end

        %Try to include GC clamp if possible;
        display(['Original C R1 primer: ',p])
        %no need to take the reverse complement here since the gene is already in reverse orientation

        %Case #1: C R1 primer has GC clamp
        if regexp(p(end-1:end), '^[gc]{2}')
            display('C R1 primer has GC clamp.') % --> do nothing

        %Case #2: C R1 primer has G/C at first position but A/T at second position    
        elseif regexp(p(end-1:end), '^[at][gc]')
            display('C R1 primer has G/C at first position but A/T at second position.')                        
            p_new = g.Sequence(g.CDS(h).indices(end)+3 : g.CDS(h).indices(end)+3-1+(pl+1));%increase primer lenght by 1
            display(['C R1 primer (+1): ',p_new])
            
            if regexp(p_new(end), '^[gc]') % new C R1 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %use old seqeunce; C R1 primer has G/C at first position
            end

        %Case #3: C R1 primer has A/T at first position but G/C at second last position    
        elseif regexp(p(end-1:end), '^[gc][at]')
            display('C R1 primer has A/T at first position but G/C at second position.')                           
            p_new = g.Sequence(g.CDS(h).indices(end)+3 : g.CDS(h).indices(end)+3-1+(pl+2)); %increase primer lenght by 2
            display(['C R1 primer (+2): ',p_new])
            
            if regexp(p_new(end-1:end), '^[gc]{2}') % new C R1 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[at][gc]') % new C R1 primer has G/C at first position and A/T at second first position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[gc][at]') % new C R1 primer has A/T at first position and G/C at second first position
                %increase lenght of OLD C R1 primer by 1 (this is same as decreasing length of NEW C R1 primer by 1)
                p_new = g.Sequence(g.CDS(h).indices(end)+3 : g.CDS(h).indices(end)+3-1+(pl+1)); 
                display(['C R1 primer (+1): ',p_new])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else % new C R1 primer has A/T at first and second position --> use old C R1 primer
            end 

        %Case #4: C R1 primer has A/T at first and second position
        elseif regexp(p(end-1:end), '^[at]{2}')        
            display('C R1 primer has A/T at first and second position.')                        
            p_new = g.Sequence(g.CDS(h).indices(end)+3 : g.CDS(h).indices(end)+3-1+(pl+2)); %increase primer lenght by 2
            display(['C R1 primer (+2): ',p_new])
            
            if regexp(p_new(end-1:end), '^[gc]{2}') % new C R1 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[at][gc]') % new C R1 primer has G/C at first position and A/T at second position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[gc][at]') % new C R1 primer has A/T at first position and G/C at second position
                %increase lenght of OLD C R1 primer by 1 (this is same as decreasing length of NEW C R1 primer by 1)
                p_new = g.Sequence(g.CDS(h).indices(end)+3 : g.CDS(h).indices(end)+3-1+(pl+1));
                display(['C R1 primer (+1): ',p_new])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %new C R1 primer has A/T at first and second positions
                p = p_new; %this extends the primer ending with e.g. 5'-AA to 5'-AAAA 
                Tm = calculate_Tm_like_PrimerSelect(p);
            end                    
        end

        primer(2).name = [g.CDS(h).gene,'_C_R1'];
        %no need to take the reverse complement here since the gene is already in reverse orientation        
        primer(2).overhang = seqrcomplement('GGAGCAGGTGCTGGTGCTGG');
        %check length of C R1 primer
        %primer cannot be longer than 40 nucleotide since the overhang is already 20 nucleotide and the maximal primer length is 60 nucleotide
        if size(p,2) > 60-length(primer(2).overhang)
            display(['Need to truncate CM R1 primer to ',num2str(60-length(primer(2).overhang)),' nucleotides (--> total length is 60 nucleotides).'])               
            p = p(1:60-length(primer(2).overhang)); % --> truncate to 40 nucleotides           
            Tm = calculate_Tm_like_PrimerSelect(p);
        else
        end                    
        primer(2).seq = [primer(2).overhang,p]; %add 5' appendix to CM R1 primer
        primer(2).Tm = Tm;        
        primer(2).length_binding = length(p);
        primer(2).length_overhang = length(primer(2).overhang);        
        primer(2).length_total = primer(2).length_binding + primer(2).length_overhang; %add appendix to primer length             
        primer(2).index_5p = regexp(g.Sequence,p); %index of 5'
        primer(2).index_3p = primer(2).index_5p + primer(2).length_binding - 1; %index of 3'
        primer(2) = correctPrimerBindingSite(g,h,primer(2));
        clear k i pl p p_new Tm s
                               
        
        % ---------- Third primer (C F2 primer) ----------
        display('----- Designing C F2 primer -----')
        pl = 20;        
        p = seqrcomplement(g.Sequence(g.CDS(h).indices(end)-pl+3 : g.CDS(h).indices(end)-1+3)); %this is now a reverse primer 
        Tm = calculate_Tm_like_PrimerSelect(p);
        
        while Tm < 58.0
            pl = pl + 1; %increase primer lenght by 1            
            p = seqrcomplement(g.Sequence(g.CDS(h).indices(end)-pl+3 : g.CDS(h).indices(end)-1+3));
            Tm = calculate_Tm_like_PrimerSelect(p);
        end

        %Try to include GC clamp in C F2 primer if possible
        display(['Original C F2 primer: ',p])

        %Case #1: C F2 primer has already GC clamp
        if regexp(p(end-1:end), '^[gc]{2}')
            display('C F2 primer has GC clamp.') % --> do nothing

        %Case #2: C F2 primer has G/C at first position but A/T at second position            
        elseif regexp(p(end-1:end), '^[at][gc]')
            display('C F2 primer has G/C at first position but A/T at second position.')                               
            p_new = seqrcomplement(g.Sequence(g.CDS(h).indices(end)-(pl+1)+3 : g.CDS(h).indices(end)-1+3)); %increase primer lenght by 1
            display(['C F2 primer (+1): ',p_new])
            
            if regexp(p_new(end), '^[gc]') % new primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %use old primer seqeunce; C F2 primer has G/C at first position
            end

        %Case #3: C F2 primer has A/T at first position but G/C at second last position            
        elseif regexp(p(end-1:end), '^[gc][at]')
            display('C F2 primer has A/T at first position but G/C at second position.')                                           
            p_new = seqrcomplement(g.Sequence(g.CDS(h).indices(end)-(pl+2)+3 : g.CDS(h).indices(end)-1+3)); %increase primer lenght by 2
            display(['C F2 primer (+2): ',p_new])
            
            if regexp(p_new(end-1:end), '^[gc]{2}') % new primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[at][gc]') % new primer has G/C at first position and A/T at second first position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[gc][at]') % new primer has A/T at first position and G/C at second first position            
                %increase lenght of OLD C F2 primer by 1 (this is same as decreasing length of NEW C F2 primer by 1)                
                p_new = seqrcomplement(g.Sequence(g.CDS(h).indices(end)-(pl+1)+3 : g.CDS(h).indices(end)-1+3)); 
                display(['CM F2 primer (+1): ',p_new])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else % new C F2 primer has A/T at first and second position --> use old C F2 primer
            end         

        %Case #4: C F2 primer has A/T at first and second position
        elseif regexp(p(end-1:end), '^[at]{2}')        
            display('C F2 primer has A/T at first and second position.')                                            
            p_new = seqrcomplement(g.Sequence(g.CDS(h).indices(end)-(pl+2)+3 : g.CDS(h).indices(end)-1+3)); %increase primer lenght by 2
            display(['C F2 primer (+2): ',p_new])
            
            if regexp(p_new(end-1:end), '^[gc]{2}') % new primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[at][gc]') % new primer has G/C at first position and A/T at second position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[gc][at]') % new primer has A/T at first position and G/C at second position            
            %increase lenght of OLD C F2 primer by 1 (this is same as decreasing length of NEW C F2 primer by 1)             
            p_new = seqrcomplement(g.Sequence(g.CDS(h).indices(end)-(pl+1)+3 : g.CDS(h).indices(end)-1+3));
            display(['C F2 primer (+1): ',p_new])
            p = p_new;
            Tm = calculate_Tm_like_PrimerSelect(p);
            else %new C F2 primer has A/T at first and second positions
                p = p_new; %this extends the primer ending with e.g. 5'-AA to 5'-AAAA 
                Tm = calculate_Tm_like_PrimerSelect(p);
            end                    
        end

        primer(3).name = [g.CDS(h).gene,'_C_F2'];                
        primer(3).overhang = seqrcomplement('CTTGTACAATTCGTCCATACCCATAACG'); %using e.g. pDML219
        %check length of C F2 primer
        %primer cannot be longer than 35 nucleotide since the overhang is already 28 nucleotide and the maximal primer length is 60 nucleotide
        if size(p,2) > 60-length(primer(3).overhang)
            display(['Need to truncate C F2 primer to ',num2str(60-length(primer(3).overhang)),' nucleotides (--> total length is 60 nucleotides).'])                    
            p = p(1:60-length(primer(3).overhang)); % --> truncate to 32 nucleotides           
            Tm = calculate_Tm_like_PrimerSelect(p);
        else
        end                     
        primer(3).seq = [primer(3).overhang,p]; %add 5' appendix
        primer(3).Tm = Tm;
        primer(3).length_binding = length(p); %add appendix to primer length                
        primer(3).length_overhang = length(primer(3).overhang); %add appendix to primer length                
        primer(3).length_total = primer(3).length_overhang + primer(3).length_binding;
        primer(3).index_3p = regexp(g.Sequence,seqrcomplement(p));
        primer(3).index_5p = primer(3).index_3p + primer(3).length_binding - 1;
        primer(3) = correctPrimerBindingSite(g,h,primer(3));
        clear k i pl p p_new Tm s


        % ---------- Fourth primer (C R2 primer) ----------
        display('----- Designing C R2 primer -----')    
        k = regexp(downstream, '[c|g]{2}'); %find all GC clamps
        k = sort(k,'ascend'); %start with the primer closest to stop codon and move further downstream from there    
        
        for i = 1:size(k,2)   
            pl = 20;
            p = downstream(k(i) : k(i)+pl-1);    
            Tm = calculate_Tm_like_PrimerSelect(p); 
            
            while Tm < 58.0
                pl = pl +1;                   
                p = downstream(k(i) : k(i)+pl-1);
                Tm = calculate_Tm_like_PrimerSelect(p);                        
            end    

            s = oligoprop(p);
            if isempty(s.Hairpins)                        
                primer(4).name = [g.CDS(h).gene,'_C_R2'];
                primer(4).seq = seqrcomplement(p);
                primer(4).Tm = Tm;
                primer(4).overhang = '';
                primer(4).length_binding = length(p);            
                primer(4).length_overhang = length(primer(4).overhang);
                primer(4).length_total = primer(4).length_overhang + primer(4).length_binding;            
                primer(4).index_5p = regexp(g.Sequence,primer(4).seq); %index of 5'
                primer(4).index_3p = primer(4).index_5p + primer(4).length_binding - 1;
                break
            else
            end      
        end
        primer(4) = correctPrimerBindingSite(g,h,primer(4));
        clear k i pl p Tm s

        
         % ---------- Calculate length of PCR product H1_C ----------        
        H1_length = primer(1).index_5p - primer(2).index_5p + 1 + primer(2).length_overhang;
        display(['Length of PCR product H1_C is ',num2str(H1_length),' bp.'])
                
         % ---------- Calculate length of PCR product H2_C ----------        
        H2_length = primer(3).index_5p - primer(4).index_5p + 1 + primer(3).length_overhang;
        display(['Length of PCR product H2_C is ',num2str(H2_length),' bp.'])
                        
        %write PCR product lenght of H1_C and H2_C to a text file        
        filename_out3 = [datestr(now,'yyyymmdd_THHMMSS'),'_length_of_H1_C_and_H2_C_PCR_products_for_',g.CDS(h).gene,'.txt']; display(filename_out3);
        fid = fopen(filename_out3, 'w','n','US-ASCII'); %write mode
        fprintf(fid, '%s\t%s\t%s\t\n','Gene','Name','Length (bp)');        
        fprintf(fid, '%s\t%s\t%d\t\n',g.CDS(h).gene,['H1_C_',g.CDS(h).gene], H1_length);
        fprintf(fid, '%s\t%s\t%d\t\n',g.CDS(h).gene,['H2_C_',g.CDS(h).gene], H2_length);
        fclose (fid);
        

        % ---------- Fifth primer: forward primer for colony PCR (C Fc primer) ----------
        display('----- Designing C Fc primer -----')           
        k = strfind(upstream,primer(1).seq);
        upstream2 = upstream(1:k-1); %region upstream of the CM F1 primer; k(i) comes from the calculation of the CM F1 primer        
        k = regexp(upstream2, '[c|g]{2}'); %find all GC clamps in the region upstream of the CM F1 primer
        k = sort(k,'descend'); %start with the primer closest to the CM F1 primer and move further upstream from there
        
        for i = 1:size(k,2) %loop over all potential primers    
            pl = 20; %minimal primer length
            p = upstream2(k(i)-pl+2 : k(i)+1); %size(p,2) = 20   
            Tm = calculate_Tm_like_PrimerSelect(p);
        
            while Tm < 58.0 %Tm of primer needs to be above 58.0 C
                pl = pl + 1; %increase primer length by 1       
                p = upstream2(k(i)-pl+2 : k(i)+1);
                Tm = calculate_Tm_like_PrimerSelect(p);
            end

            s = oligoprop(p); %calculate properties of primer
            if isempty(s.Hairpins)                        
                primer(5).name = [g.CDS(h).gene,'_C_Fc'];
                primer(5).seq = p;
                primer(5).Tm = Tm;
                primer(5).overhang = '';
                primer(5).length_binding = length(p);    
                primer(5).length_overhang = length(primer(5).overhang);    
                primer(5).length_total = primer(5).length_overhang + primer(5).length_binding;                
                primer(5).index_3p = regexp(g.Sequence,seqrcomplement(p));
                primer(5).index_5p = primer(5).index_3p + primer(5).length_binding -1;                                
                break
            else
            end      
        end
        primer(5) = correctPrimerBindingSite(g,h,primer(5));
        clear k i pl p Tm s
                      
    else
    end 
        
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------ METHOD #3 ----- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'N-terminal'    
    
    
    %First determine whether the gene is on the upper (sense) or lower (anti-sense) strand
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- DNA sense strand (gene has normal orientation) ---
    if g.CDS(h).indices(end) > g.CDS(h).indices(1)         
        display([g.CDS(h).gene, ' gene has normal orientation.'])    
        ORF = g.Sequence(g.CDS(h).indices(1) : g.CDS(h).indices(end)-3); %does not contain the STOP codon
        start_codon = g.Sequence(g.CDS(h).indices(1) : g.CDS(h).indices(1)+2); display(['Start codon: ',upper(start_codon)])
        stop_codon = g.Sequence(g.CDS(h).indices(end)-2 : g.CDS(h).indices(end)); display(['Stop codon: ',upper(stop_codon)])       
        
        try
        %500-bp search window for F1 primer; search window starts 280 bp upstream of start codon
        upstream = g.Sequence(g.CDS(h).indices(1)-780 : g.CDS(h).indices(1)-281); 
        %500-bp search window for R2 primer; search window starts 280 bp downstream of start codon        
        downstream = g.Sequence(g.CDS(h).indices(1)+280 : g.CDS(h).indices(1)+779);
        catch ME
            display(['Upstream or downstream sequence(s) are not long enough. Skipped primer design for gene ',g.CDS(h).gene,'.'])
            primer = addEmptyPrimerEntries(g.CDS(h).gene);
            return
        end

        % ---------- First primer (N F1 primer) ----------
        display('----- Designing N F1 primer -----')    
        k = regexp(upstream, '[c|g]{2}'); %find all GC clamps
        k = sort(k,'descend'); %start with the primer closest to start codon (i.e. insertion site) and move further upstream from there

        for i = 1:size(k,2) %loop over all potential primers
            pl = 20; %minimal primer length
            p = upstream(k(i)-pl+2 : k(i)+1); %size(p,2) = 20   
            Tm = calculate_Tm_like_PrimerSelect(p);
            
            while Tm < 58.0 %Tm of primer needs to be above 58.0 C
                pl = pl +1; %increase primer length by 1       
                p = upstream(k(i)-pl+2 : k(i)+1);
                Tm = calculate_Tm_like_PrimerSelect(p);
            end

            s = oligoprop(p); %calculate properties of primer
            if isempty(s.Hairpins)                        
                primer(1).name = [g.CDS(h).gene,'_N_F1'];
                primer(1).seq = p;
                primer(1).Tm = Tm;                                
                primer(1).overhang = '';
                primer(1).length_binding = length(p);                        
                primer(1).length_overhang = 0;
                primer(1).length_total = primer(1).length_overhang + primer(1).length_binding;                
                primer(1).index_5p = regexp(g.Sequence,primer(1).seq);   
                primer(1).index_3p = primer(1).index_5p + primer(1).length_binding - 1;      
                orderfields(primer(1),{'name','seq','Tm','overhang','length_binding','length_overhang','length_total','index_5p','index_3p'});
                break
            else
            end      
        end
        primer(1) = correctPrimerBindingSite(g,h,primer(1));
        clear k i pl p Tm s
        %Note, not necessary to check whether N F1 primer is larger than 60 nt since primer with 35x As or 35x Ts has a Tm of 58.1436 C              


        % ---------- Second primer (N R1 primer) ----------
        display('----- Designing N R1 primer -----')    
        pl = 20; %minimal primer length           
        p = g.Sequence(g.CDS(h).indices(1)-pl+3 : g.CDS(h).indices(1)+3-1); 
        Tm = calculate_Tm_like_PrimerSelect(p);
        
        while Tm < 58.0
            pl = pl +1; %increase primer lenght by 1            
            p = g.Sequence(g.CDS(h).indices(1)-pl+3 : g.CDS(h).indices(1)+3-1); 
            Tm = calculate_Tm_like_PrimerSelect(p);
        end

        %Try to include GC clamp if possible
        display(['Original N R1 primer: ',seqrcomplement(p)])

        %Case #1: N R1 primer has GC clamp
        if regexp(p(1:2), '^[gc]{2}')
            display('N R1 primer has GC clamp.') % --> do nothing

        %Case #2: N R1 primer has G/C at first position but A/T at second position    
        elseif regexp(p(1:2), '^[gc][at]')
            display('N R1 primer has G/C at first position but A/T at second position.')                            
             p_new = g.Sequence(g.CDS(h).indices(1)-(pl+1)+3 : g.CDS(h).indices(1)+3-1); %increase primer lenght by 1
            display(['N R1 primer (+1): ',seqrcomplement(p_new)])
            
            if regexp(p_new(1), '^[gc]') % new N R1 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %use old seqeunce; N R1 primer has G/C at first position
            end

        %Case #3: N R1 primer has A/T at first position but G/C at second last position    
        elseif regexp(p(1:2), '^[at][gc]')
            display('N R1 primer has A/T at first position but G/C at second position.')                           
            p_new = g.Sequence(g.CDS(h).indices(1)-(pl+2)+3 : g.CDS(h).indices(1)+3-1); %increase primer lenght by 2
            display(['N R1 primer (+2): ',seqrcomplement(p_new)])
            
            if regexp(p_new(1:2), '^[gc]{2}') % new N R1 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(1:2), '^[gc][at]') % new N R1 primer has G/C at first position and A/T at second first position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(1:2), '^[at][gc]') % new N R1 primer has A/T at first position and G/C at second first position
                %increase lenght of OLD N R1 primer by 1 (this is same as decreasing length of NEW N R1 primer by 1)                
                p_new = g.Sequence(g.CDS(h).indices(1)-(pl+1)+3 : g.CDS(h).indices(1)+3-1);
                display(['N R1 primer (+1): ',seqrcomplement(p_new)])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else % new N R1 primer has A/T at first and second position --> use old N R1 primer
            end 

        %Case #4: N R1 primer has A/T at first and second position
        elseif regexp(p(1:2), '^[at]{2}')        
            display('N R1 primer has A/T at first and second position.')                            
            p_new = g.Sequence(g.CDS(h).indices(1)-(pl+2)+3 : g.CDS(h).indices(1)+3-1); %increase primer lenght by 2
            display(['N R1 primer (+2): ',seqrcomplement(p_new)])
            
            if regexp(p_new(1:2), '^[gc]{2}') % new N R1 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(1:2), '^[gc][at]') % new N R1 primer has G/C at first position and A/T at second position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(1:2), '^[at][gc]') % new N R1 primer has A/T at first position and G/C at second  position
                %increase lenght of OLD N R1 primer by 1 (this is same as decreasing length of NEW N R1 primer by 1)                
                p_new = g.Sequence(g.CDS(h).indices(1)-(pl+1)+3 : g.CDS(h).indices(1)+3-1); 
                display(['N R1 primer (+1): ',seqrcomplement(p_new)])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %new N R1 primer has A/T at first and second positions
                p = p_new; %this extends the primer ending with e.g. 5'-AA to 5'-AAAA 
                Tm = calculate_Tm_like_PrimerSelect(p);               
            end                    
        end

        %N R1 is a reverse primer --> take reverse complement of primer sequence        
        primer(2).name = [g.CDS(h).gene,'_N_R1'];                
        primer(2).overhang = 'GTTTCTAAGGGTGAAGAAGACAACATG'; %works with pDML190
        %check length of N R1 primer
        %primer cannot be longer than 33 nucleotide since the overhang is already 27 nucleotide and the maximal primer length is 60 nucleotide
        if size(p,2) > 60-length(primer(2).overhang)
            display(['Need to truncate N R1 primer to ',num2str(60-length(primer(2).overhang)),' nucleotides (--> total length is 60 nucleotides).'])                      
            p = p(end-(60-length(primer(2).overhang))+1:end); % --> truncate to 33 nucleotides           
            Tm = calculate_Tm_like_PrimerSelect(p); %calculate the new Tm
        else
        end            
        primer(2).seq = [seqrcomplement(primer(2).overhang),seqrcomplement(p)]; %add 5' appendix to N R1 primer
        primer(2).Tm = Tm;        
        primer(2).length_binding = length(p);
        primer(2).length_overhang = length(primer(2).overhang);        
        primer(2).length_total = primer(2).length_binding + primer(2).length_overhang; %add appendix to primer length                     
        primer(2).index_3p = regexp(g.Sequence,p); %index of 5'
        primer(2).index_5p = primer(2).index_3p + primer(2).length_binding - 1; %index of 3'                        
        primer(2) = correctPrimerBindingSite(g,h,primer(2));
        clear k i pl p p_new Tm s


        % ---------- Third primer (N F2 primer) ----------
        display('----- Designing N F2 primer -----')
        pl = 20;        
        p = g.Sequence(g.CDS(h).indices(1)+3 : g.CDS(h).indices(1)+pl+3-1); 
        Tm = calculate_Tm_like_PrimerSelect(p);

        while Tm < 58.0
            pl = pl +1; %increase primer lenght by 1            
            p = g.Sequence(g.CDS(h).indices(1)+3 : g.CDS(h).indices(1)+pl+3-1);
            Tm = calculate_Tm_like_PrimerSelect(p);
        end

        %Try to include GC clamp in N F2 primer if possible
        display(['Original N F2 primer: ',p])

        %Case #1: N F2 primer has already GC clamp
        if regexp(p(end-1:end), '^[gc]{2}')
            display('N F2 primer has GC clamp.') % --> do nothing

        %Case #2: N F2 primer has G/C at first position but A/T at second position        
        elseif regexp(p(end-1:end), '^[at][gc]')
            display('N F2 primer has G/C at first position but A/T at second position.')                           
            p_new = g.Sequence(g.CDS(h).indices(1)+3 : g.CDS(h).indices(1)+(pl+1)+3-1); %increase primer lenght by 1
            display(['N F2 primer (+1): ',p_new])
            
            if regexp(p_new(end), '^[gc]') % new N F2 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %use old primer seqeunce; N F2 primer has G/C at first position
            end

        %Case #3: N F2 primer has A/T at first position but G/C at second last position        
        elseif regexp(p(end-1:end), '^[gc][at]')
            display('N F2 Primer has A/T at first position but G/C at second position.')                                           
            p_new = g.Sequence(g.CDS(h).indices(1)+3 : g.CDS(h).indices(1)+(pl+2)+3-1); %increase primer lenght by 2
            display(['N F2 primer (+2): ',p_new])
            
            if regexp(p_new(end-1:end), '^[gc]{2}') % new N F2 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[at][gc]') % new N F2 primer has G/C at first position and A/T at second first position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[gc][at]') % new N F2 primer has A/T at first position and G/C at second first position            
                %increase lenght of OLD N F2 primer by 1 (this is same as decreasing length of NEW N F2 primer by 1)                       
                p_new = g.Sequence(g.CDS(h).indices(1)+3 : g.CDS(h).indices(1)+(pl+1)+3-1); 
                display(['N F2 primer (+1): ',p_new])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else % new N F2 primer has A/T at first and second position --> use old primer
            end 

        %Case #4: N F2 primer has A/T at first and second position
        elseif regexp(p(end-1:end), '^[at]{2}')        
            display('N F2 primer has A/T at first and second position.')                                    
            p_new = g.Sequence(g.CDS(h).indices(1)+3 : g.CDS(h).indices(1)+(pl+2)+3-1); %increase primer lenght by 2
            display(['N F2 primer (+2): ',p_new])
            
            if regexp(p_new(end-1:end), '^[gc]{2}') % new N F2 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[at][gc]') % new N F2 primer has G/C at first position and A/T at second position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[gc][at]') % new N F2 primer has A/T at first position and G/C at second  position            
                %increase lenght of OLD N F2 primer by 1 (this is same as decreasing length of NEW N F2 primer by 1)                
                p_new = g.Sequence(g.CDS(h).indices(1)+3 : g.CDS(h).indices(1)+(pl+1)+3-1); 
                display(['N F2 primer (+1): ',p_new])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %new N F2 primer has A/T at first and second positions
                p = p_new; %this extends the primer ending with e.g. 5'-AA to 5'-AAAA 
                Tm = calculate_Tm_like_PrimerSelect(p); 
            end                    
        end

        primer(3).name = [g.CDS(h).gene,'_N_F2'];                
        primer(3).overhang = 'TGCACCAGCACCAGCACCAGC';        
        %check length of N F2 primer
        %primer cannot be longer than 39 nucleotide since the overhang is already 21 nucleotide and the maximal primer length is 60 nucleotide
        if size(p,2) > 60-length(primer(3).overhang)
            display(['Need to truncate N F2 primer to ',num2str(60-length(primer(3).overhang)),' nucleotides (--> total length is 60 nucleotides).'])
            p = p(1:60-length(primer(3).overhang)); % --> truncate to 39 nucleotides           
            Tm = calculate_Tm_like_PrimerSelect(p);
        else
        end                       
        primer(3).seq = [seqrcomplement(primer(3).overhang),p]; %add 5' appendix
        primer(3).Tm = Tm;
        primer(3).length_binding = length(p); %add appendix to primer length                
        primer(3).length_overhang = length(primer(3).overhang); %add appendix to primer length                
        primer(3).length_total = primer(3).length_overhang + primer(3).length_binding;
        primer(3).index_5p = regexp(g.Sequence,p);
        primer(3).index_3p = primer(3).index_5p + primer(3).length_binding - 1;                        
        primer(3) = correctPrimerBindingSite(g,h,primer(3));
        clear k i pl p p_new Tm s


        % ---------- Fourth primer (N R2 primer) ----------
        display('----- Designing N R2 primer -----')    
        k = regexp(downstream, '[c|g]{2}'); %find all GC clamps
        k = sort(k,'ascend'); %start with the primer closest to stop codon and move further downstream from there    

        for i = 1:size(k,2)
            pl = 20;
            p = downstream(k(i) : k(i)+pl-1);                
            Tm = calculate_Tm_like_PrimerSelect(p);    
            
            while Tm < 58.0        
                pl = pl +1;                        
                p = downstream(k(i) : k(i)+pl-1);
                Tm = calculate_Tm_like_PrimerSelect(p);                
            end    

            s = oligoprop(p);
            if isempty(s.Hairpins)        
                primer(4).name = [g.CDS(h).gene,'_N_R2'];                
                primer(4).seq = seqrcomplement(p);
                primer(4).Tm = Tm;
                primer(4).overhang = '';
                primer(4).length_binding = length(p);            
                primer(4).length_overhang = length(primer(4).overhang);
                primer(4).length_total = primer(4).length_overhang + primer(4).length_binding;            
                primer(4).index_3p = regexp(g.Sequence,seqrcomplement(primer(4).seq)); %index of 5'
                primer(4).index_5p = primer(4).index_3p + primer(4).length_binding - 1;                                
                break
            else
            end      
        end
        primer(4) = correctPrimerBindingSite(g,h,primer(4));
        clear k i pl p Tm s


        % ---------- Calculate length of PCR product H1_N ----------        
        H1_length = primer(2).index_5p + 1 + primer(2).length_overhang - primer(1).index_5p;
        display(['Length of PCR product H1_N is ',num2str(H1_length),' bp.'])
                
         % ---------- Calculate length of PCR product H2_N ----------        
        H2_length = primer(4).index_5p - primer(3).index_5p + 1 + primer(3).length_overhang;
        display(['Length of PCR product H2_N is ',num2str(H2_length),' bp.'])
                
        %write PCR product lenght of H1_N and H2_N to a text file        
        filename_out3 = [datestr(now,'yyyymmdd_THHMMSS'),'_length_of_H1_N_and_H2_N_PCR_products_for_',g.CDS(h).gene,'.txt']; display(filename_out3);
        fid = fopen(filename_out3, 'w','n','US-ASCII'); %write mode
        fprintf(fid, '%s\t%s\t%s\t\n','Gene','Name','Length (bp)');        
        fprintf(fid, '%s\t%s\t%d\t\n',g.CDS(h).gene,['H1_N_',g.CDS(h).gene], H1_length);
        fprintf(fid, '%s\t%s\t%d\t\n',g.CDS(h).gene,['H2_N_',g.CDS(h).gene], H2_length);
        fclose (fid);
        
        
        % ---------- Fifth primer: forward primer for colony PCR (N Fc primer) ----------
        display('----- Designing N Fc primer -----')           
        k = strfind(upstream,primer(1).seq);
        upstream2 = upstream(1:k-1); %region upstream of the N F1 primer;
        k = regexp(upstream2, '[c|g]{2}'); %find all GC clamps in the region upstream of the N F1 primer
        k = sort(k,'descend'); %start with the primer closest to the N F1 primer and move further upstream from there

        for i = 1:size(k,2) %loop over all potential primers    
            pl = 20; %minimal primer length
            p = upstream2(k(i)-pl+2 : k(i)+1); %size(p,2) = 20   
            Tm = calculate_Tm_like_PrimerSelect(p);
            
            while Tm < 58.0 %Tm of primer needs to be above 58.0 C
                pl = pl +1; %increase primer length by 1       
                p = upstream2(k(i)-pl+2 : k(i)+1);
                Tm = calculate_Tm_like_PrimerSelect(p);
            end

            s = oligoprop(p); %calculate properties of primer
            if isempty(s.Hairpins)                        
                primer(5).name = [g.CDS(h).gene,'_N_Fc'];
                primer(5).seq = p;
                primer(5).Tm = Tm;
                primer(5).overhang = '';
                primer(5).length_binding = length(p);    
                primer(5).length_overhang = length(primer(5).overhang);    
                primer(5).length_total = primer(5).length_overhang + primer(5).length_binding;                                
                primer(5).index_5p = regexp(g.Sequence,p);                                   
                primer(5).index_3p = primer(5).index_5p + primer(5).length_binding - 1;                                
                break
            else
            end      
        end
        primer(5) = correctPrimerBindingSite(g,h,primer(5));
        clear k i pl p Tm s

                               
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % --- DNA antisense strand (gene has reverse orientation) ---
    elseif g.CDS(h).indices(end) < g.CDS(h).indices(1)         
        display([g.CDS(h).gene, ' gene has reverse orientation.'])    

        ORF = seqrcomplement(g.Sequence(g.CDS(h).indices(end)+3 : g.CDS(h).indices(1)));    
        start_codon = seqrcomplement(g.Sequence(g.CDS(h).indices(1)-2 : g.CDS(h).indices(1))); display(['Start codon: ',upper(start_codon)])
        stop_codon = seqrcomplement(g.Sequence(g.CDS(h).indices(end) : g.CDS(h).indices(end)+3-1)); display(['Stop codon: ',upper(stop_codon)])                

        try
        %500-bp search window for N F1 primer; search window starts 280 bp upstream of stop codon              
        upstream = seqrcomplement(g.Sequence(g.CDS(h).indices(1)+278 : g.CDS(h).indices(1)+777));
        %500-bp search window for N R2 primer; search window starts 280 bp upstream of stop codon        
        downstream = seqrcomplement(g.Sequence(g.CDS(h).indices(1)-782 : g.CDS(h).indices(1)-283));
        catch ME
            display(['Upstream or downstream sequence(s) are not long enough. Skipped primer design for gene ',g.CDS(h).gene,'.'])
            primer = addEmptyPrimerEntries(g.CDS(h).gene);
            return
        end

        % ---------- First primer (N F1 primer) ----------
        display('----- Designing N F1 primer -----')    
        k = regexp(upstream, '[c|g]{2}'); %find all GC clamps
        k = sort(k,'descend'); %start with the primer closest to stop codon and move further upstream from there
        
        for i = 1:size(k,2) %loop over all potential primers    
            pl = 20; %minimal primer length
            p = upstream(k(i)-pl+2 : k(i)+1); %size(p,2) = 20       
            Tm = calculate_Tm_like_PrimerSelect(p);
            
            while Tm < 58.0 %Tm of primer needs to be above 58.0 C
                pl = pl +1; %increase primer length by 1       
                p = upstream(k(i)-pl+2 : k(i)+1);
                Tm = calculate_Tm_like_PrimerSelect(p);
            end

            s = oligoprop(p); %calculate properties of primer
            if isempty(s.Hairpins)                        
                primer(1).name = [g.CDS(h).gene,'_N_F1'];
                primer(1).seq = p;
                primer(1).Tm = Tm;
                primer(1).overhang = '';
                primer(1).length_binding = length(p);                        
                primer(1).length_overhang = 0;
                primer(1).length_total = primer(1).length_overhang + primer(1).length_binding;
                primer(1).index_3p = regexp(g.Sequence,seqrcomplement(primer(1).seq)); 
                primer(1).index_5p = primer(1).index_3p + primer(1).length_binding - 1;
                orderfields(primer(1),{'name','seq','Tm','overhang','length_binding','length_overhang','length_total','index_5p','index_3p'});                
                break
            else
            end      
        end
        primer(1) = correctPrimerBindingSite(g,h,primer(1));
        clear k i pl p Tm s
        %Note, not necessary to check whether N F1 primer is larger than 60 nt since primer with 35x As or 35x Ts has a Tm of 58.1436 C  


        % ---------- Second primer (N R1 primer) ----------
        display('----- Designing N R1 primer -----')                    
        pl = 20; %minimal primer length            
        p = g.Sequence(g.CDS(h).indices(1)-3+1 : g.CDS(h).indices(1)-3+pl); 
        Tm = calculate_Tm_like_PrimerSelect(p);
        
        while Tm < 58.0
            pl = pl + 1; %increase primer lenght by 1                  
            p = g.Sequence(g.CDS(h).indices(1)-3+1 : g.CDS(h).indices(1)-3+pl);
            Tm = calculate_Tm_like_PrimerSelect(p);
        end

        %Try to include GC clamp if possible
        display(['Original N R1 primer: ',p])
        %no need to take the reverse complement here since the gene is already in reverse orientation

        %Case #1: N R1 primer has GC clamp
        if regexp(p(end-1:end), '^[gc]{2}')
            display('N R1 primer has GC clamp.') % --> do nothing

        %Case #2: N R1 primer has G/C at first position but A/T at second position    
        elseif regexp(p(end-1:end), '^[at][gc]')
            display('N R1 primer has G/C at first position but A/T at second position.')                                    
            p_new = g.Sequence(g.CDS(h).indices(1)-3+1 : g.CDS(h).indices(1)-3+(pl+1)); %increase primer lenght by 1
            display(['N R1 primer (+1): ',p_new])
            
            if regexp(p_new(end), '^[gc]') % new N R1 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %use old seqeunce; N R1 primer has G/C at first position
            end

        %Case #3: N R1 primer has A/T at first position but G/C at second last position    
        elseif regexp(p(end-1:end), '^[gc][at]')
            display('N R1 primer has A/T at first position but G/C at second position.')                                       
            p_new = g.Sequence(g.CDS(h).indices(1)-3+1 : g.CDS(h).indices(1)-3+(pl+2)); %increase primer lenght by 2
            display(['N R1 primer (+2): ',p_new])
            
            if regexp(p_new(end-1:end), '^[gc]{2}') % new N R1 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[at][gc]') % new N R1 primer has G/C at first position and A/T at second first position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[gc][at]') % new N R1 primer has A/T at first position and G/C at second first position
                %increase lenght of OLD N R1 primer by 1 (this is same as decreasing length of NEW C R1 primer by 1)                
                p_new = g.Sequence(g.CDS(h).indices(1)-3+1 : g.CDS(h).indices(1)-3+(pl+1)); 
                display(['N R1 primer (+1): ',p_new])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else % new N R1 primer has A/T at first and second position --> use old N R1 primer
            end 

        %Case #4: N R1 primer has A/T at first and second position
        elseif regexp(p(end-1:end), '^[at]{2}')        
            display('N R1 primer has A/T at first and second position.')                                    
            p_new = g.Sequence(g.CDS(h).indices(1)-3+1 : g.CDS(h).indices(1)-3+(pl+2)); %increase primer lenght by 2
            display(['N R1 primer (+2): ',p_new])
            
            if regexp(p_new(end-1:end), '^[gc]{2}') % new N R1 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[at][gc]') % new N R1 primer has G/C at first position and A/T at second position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[gc][at]') % new N R1 primer has A/T at first position and G/C at second  position
                %increase lenght of OLD N R1 primer by 1 (this is same as decreasing length of NEW N R1 primer by 1)                
                p_new = g.Sequence(g.CDS(h).indices(1)-3+1 : g.CDS(h).indices(1)-3+(pl+1)); 
                display(['N R1 primer (+1): ',p_new])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %new N R1 primer has A/T at first and second positions
                p = p_new; %this extends the primer ending with e.g. 5'-AA to 5'-AAAA 
                Tm = calculate_Tm_like_PrimerSelect(p);
            end                    
        end

        primer(2).name = [g.CDS(h).gene,'_N_R1'];
        %no need to take the reverse complement here since the gene is already in reverse orientation
        primer(2).overhang = 'GTTTCTAAGGGTGAAGAAGACAACATG'; %works with pDML190
        %check length of N R1 primer
        %primer cannot be longer than 33 nucleotide since the overhang is already 27 nucleotide and the maximal primer length is 60 nucleotide
        if size(p,2) > 60-length(primer(2).overhang)
            display(['Need to truncate N R1 primer to ',num2str(60-length(primer(2).overhang)),' nucleotides (--> total length is 60 nucleotides).'])             
            p = p(1:60-length(primer(2).overhang)); % --> truncate to 33 nucleotides           
            Tm = calculate_Tm_like_PrimerSelect(p);
        else
        end                    
        primer(2).seq = [seqrcomplement(primer(2).overhang),p]; %add 5' appendix to N R1 primer
        primer(2).Tm = Tm;        
        primer(2).length_binding = length(p);
        primer(2).length_overhang = length(primer(2).overhang);        
        primer(2).length_total = primer(2).length_binding + primer(2).length_overhang; %add appendix to primer length             
        primer(2).index_5p = regexp(g.Sequence,p); %index of 5'
        primer(2).index_3p = primer(2).index_5p + primer(2).length_binding - 1; %index of 3'
        primer(2) = correctPrimerBindingSite(g,h,primer(2));
        clear k i pl p p_new Tm s
                               
        
        % ---------- Third primer (N F2 primer) ----------
        display('----- Designing N F2 primer -----')
        pl = 20;        
        p = seqrcomplement(g.Sequence(g.CDS(h).indices(1)-3+1-pl : g.CDS(h).indices(1)-3)); %this is now a reverse primer
        Tm = calculate_Tm_like_PrimerSelect(p);
        
        while Tm < 58.0
            pl = pl + 1; %increase primer lenght by 1            
            p = seqrcomplement(g.Sequence(g.CDS(h).indices(1)-3+1-pl : g.CDS(h).indices(1)-3));
            Tm = calculate_Tm_like_PrimerSelect(p);
        end

        %Try to include GC clamp in N F2 primer if possible
        display(['Original N F2 primer: ',p])

        %Case #1: N F2 primer has already GC clamp
        if regexp(p(end-1:end), '^[gc]{2}')
            display('N F2 primer has GC clamp.') % --> do nothing

        %Case #2: N F2 primer has G/C at first position but A/T at second position            
        elseif regexp(p(end-1:end), '^[at][gc]')
            display('N F2 primer has G/C at first position but A/T at second position.')                               
            p_new = seqrcomplement(g.Sequence(g.CDS(h).indices(1)-3+1-(pl+1) : g.CDS(h).indices(1)-3)); %increase primer lenght by 1
            display(['N F2 primer (+1): ',p_new])
            
            if regexp(p_new(end), '^[gc]') % new primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %use old primer seqeunce; N F2 primer has G/C at first position
            end

        %Case #3: N F2 primer has A/T at first position but G/C at second last position            
        elseif regexp(p(end-1:end), '^[gc][at]')
            display('N F2 primer has A/T at first position but G/C at second position.')                                           
            p_new = seqrcomplement(g.Sequence(g.CDS(h).indices(1)-3+1-(pl+2) : g.CDS(h).indices(1)-3)); %increase primer lenght by 2
            display(['N F2 primer (+2): ',p_new])
            
            if regexp(p_new(end-1:end), '^[gc]{2}') % new primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[at][gc]') % new primer has G/C at first position and A/T at second first position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[gc][at]') % new primer has A/T at first position and G/C at second first position            
                %increase lenght of OLD N F2 primer by 1 (this is same as decreasing length of NEW N F2 primer by 1)                 
                p_new = seqrcomplement(g.Sequence(g.CDS(h).indices(1)-3+1-(pl+1) : g.CDS(h).indices(1)-3));
                display(['N F2 primer (+1): ',p_new])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else % new N F2 primer has A/T at first and second position --> use old N F2 primer
            end         

        %Case #4: N F2 primer has A/T at first and second position
        elseif regexp(p(end-1:end), '^[at]{2}')        
            display('N F2 primer has A/T at first and second position.')                                            
            p_new = seqrcomplement(g.Sequence(g.CDS(h).indices(1)-3+1-(pl+2) : g.CDS(h).indices(1)-3)); %increase primer lenght by 2
            display(['N F2 primer (+2): ',p_new])
            
            if regexp(p_new(end-1:end), '^[gc]{2}') % new primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[at][gc]') % new primer has G/C at first position and A/T at second position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[gc][at]') % new primer has A/T at first position and G/C at second position            
            %increase lenght of OLD N F2 primer by 1 (this is same as decreasing length of NEW N F2 primer by 1)                     
            p_new = seqrcomplement(g.Sequence(g.CDS(h).indices(1)-3+1-(pl+1) : g.CDS(h).indices(1)-3));
            display(['N F2 primer (+1): ',p_new])
            p = p_new;
            Tm = calculate_Tm_like_PrimerSelect(p);
            else %new N F2 primer has A/T at first and second positions
                p = p_new; %this extends the primer ending with e.g. 5'-AA to 5'-AAAA 
                Tm = calculate_Tm_like_PrimerSelect(p);
            end                    
        end

        primer(3).name = [g.CDS(h).gene,'_N_F2'];
        primer(3).overhang = 'TGCACCAGCACCAGCACCAGC'; %e.g. pDML190
        %check length of N F2 primer
        %primer cannot be longer than 39 nucleotide since the overhang is already 21 nucleotide and the maximal primer length is 60 nucleotide
        if size(p,2) > 60-length(primer(3).overhang)
            display(['Need to truncate N F2 primer to ',num2str(60-length(primer(3).overhang)),' nucleotides (--> total length is 60 nucleotides).'])
            p = p(1:60-length(primer(3).overhang)); % --> truncate to 39 nucleotides           
            Tm = calculate_Tm_like_PrimerSelect(p);
        else
        end                     
        primer(3).seq = [seqrcomplement(primer(3).overhang),p]; %add 5' appendix
        primer(3).Tm = Tm;
        primer(3).length_binding = length(p); %add appendix to primer length                
        primer(3).length_overhang = length(primer(3).overhang); %add appendix to primer length                
        primer(3).length_total = primer(3).length_overhang + primer(3).length_binding;
        primer(3).index_3p = regexp(g.Sequence,seqrcomplement(p));
        primer(3).index_5p = primer(3).index_3p + primer(3).length_binding - 1;
        primer(3) = correctPrimerBindingSite(g,h,primer(3));
        clear k i pl p p_new Tm s


        % ---------- Fourth primer (N R2 primer) ----------
        display('----- Designing N R2 primer -----')    
        k = regexp(downstream, '[c|g]{2}'); %find all GC clamps
        k = sort(k,'ascend'); %start with the primer closest to start codon (i.e. insertion site) and move further downstream from there 
        
        for i = 1:size(k,2)   
            pl = 20;
            p = downstream(k(i) : k(i)+pl-1);    
            Tm = calculate_Tm_like_PrimerSelect(p); 
            
            while Tm < 58.0
                pl = pl +1;                   
                p = downstream(k(i) : k(i)+pl-1);
                Tm = calculate_Tm_like_PrimerSelect(p);                        
            end    

            s = oligoprop(p);
            if isempty(s.Hairpins)                        
                primer(4).name = [g.CDS(h).gene,'_N_R2'];
                primer(4).seq = seqrcomplement(p);
                primer(4).Tm = Tm;
                primer(4).overhang = '';
                primer(4).length_binding = length(p);            
                primer(4).length_overhang = length(primer(4).overhang);
                primer(4).length_total = primer(4).length_overhang + primer(4).length_binding;            
                primer(4).index_5p = regexp(g.Sequence,primer(4).seq); %index of 5'
                primer(4).index_3p = primer(4).index_5p + primer(4).length_binding - 1;
                break
            else
            end      
        end
        primer(4) = correctPrimerBindingSite(g,h,primer(4));
        clear k i pl p Tm s

        
         % ---------- Calculate length of PCR product H1_N ----------        
        H1_length = primer(1).index_5p - primer(2).index_5p + 1 + primer(2).length_overhang;
        display(['Length of PCR product H1_N is ',num2str(H1_length),' bp.'])
                
         % ---------- Calculate length of PCR product H2_N ----------        
        H2_length = primer(3).index_5p - primer(4).index_5p + 1 + primer(3).length_overhang;
        display(['Length of PCR product H2_N is ',num2str(H2_length),' bp.'])
                        
        %write PCR product lenght of H1_N and H2_N to a text file        
        filename_out3 = [datestr(now,'yyyymmdd_THHMMSS'),'_length_of_H1_N_and_H2_N_PCR_products_for_',g.CDS(h).gene,'.txt']; display(filename_out3);
        fid = fopen(filename_out3, 'w','n','US-ASCII'); %write mode
        fprintf(fid, '%s\t%s\t%s\t\n','Gene','Name','Length (bp)');        
        fprintf(fid, '%s\t%s\t%d\t\n',g.CDS(h).gene,['H1_N_',g.CDS(h).gene], H1_length);
        fprintf(fid, '%s\t%s\t%d\t\n',g.CDS(h).gene,['H2_N_',g.CDS(h).gene], H2_length);
        fclose (fid);
        

        % ---------- Fifth primer: forward primer for colony PCR (N Fc primer) ----------
        display('----- Designing N Fc primer -----')           
        k = strfind(upstream,primer(1).seq);
        upstream2 = upstream(1:k-1); %region upstream of the N F1 primer; k(i) comes from the calculation of the N F1 primer        
        k = regexp(upstream2, '[c|g]{2}'); %find all GC clamps in the region upstream of the N F1 primer
        k = sort(k,'descend'); %start with the primer closest to the N F1 primer and move further upstream from there
        
        for i = 1:size(k,2) %loop over all potential primers    
            pl = 20; %minimal primer length
            p = upstream2(k(i)-pl+2 : k(i)+1); %size(p,2) = 20   
            Tm = calculate_Tm_like_PrimerSelect(p);
            
            while Tm < 58.0 %Tm of primer needs to be above 58.0 C
                pl = pl +1; %increase primer length by 1       
                p = upstream2(k(i)-pl+2 : k(i)+1);
                Tm = calculate_Tm_like_PrimerSelect(p);
            end

            s = oligoprop(p); %calculate properties of primer
            if isempty(s.Hairpins)                        
                primer(5).name = [g.CDS(h).gene,'_N_Fc'];
                primer(5).seq = p;
                primer(5).Tm = Tm;
                primer(5).overhang = '';
                primer(5).length_binding = length(p);    
                primer(5).length_overhang = length(primer(5).overhang);    
                primer(5).length_total = primer(5).length_overhang + primer(5).length_binding;                
                primer(5).index_3p = regexp(g.Sequence,seqrcomplement(p));
                primer(5).index_5p = primer(5).index_3p + primer(5).length_binding -1;                                
                break
            else
            end      
        end
        primer(5) = correctPrimerBindingSite(g,h,primer(5));
        clear k i pl p Tm s
                      
    else
    end                                  
    
    
    
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------ METHOD #3 ----- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'Deletion' %gene knockout/deletion from start(ATG) to stop(TAA, TAG, TGA) 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- DNA sense strand (gene has normal orientation) ---    
    if g.CDS(h).indices(end) > g.CDS(h).indices(1)         
        display([g.CDS(h).gene, ' gene has normal orientation.'])    
        ORF = g.Sequence(g.CDS(h).indices(1) : g.CDS(h).indices(end)-3); %does not contain the STOP codon
        start_codon = g.Sequence(g.CDS(h).indices(1) : g.CDS(h).indices(1)+2); display(['Start codon: ',upper(start_codon)])
        stop_codon = g.Sequence(g.CDS(h).indices(end)-2 : g.CDS(h).indices(end)); display(['Stop codon: ',upper(stop_codon)])

        try
        %500-bp search window for KO F1 primer; search window starts 280 bp upstream of start codon
        upstream = g.Sequence(g.CDS(h).indices(1)-780 : g.CDS(h).indices(1)-281); display(upstream)
        %500-bp search window for KO R2 primer; search window starts 280 bp upstream of stop codon
        downstream = g.Sequence(g.CDS(h).indices(end)+281 : g.CDS(h).indices(end)+780); display(downstream)
        catch ME
            display(['Upstream or downstream sequence(s) are not long enough. Skipped primer design for gene ',g.CDS(h).gene,'.'])
            primer = addEmptyPrimerEntries(g.CDS(h).gene);
            return
        end
        
        % ---------- First primer (KO F1 primer) ----------
        display('----- Designing KO F1 primer -----')    
        k = regexp(upstream, '[c|g]{2}'); %find all GC clamps
        k = sort(k,'descend'); %start with the primer closest to stop codon and move further upstream from there

        for i = 1:size(k,2) %loop over all potential primers
            pl = 20; %minimal primer length
            p = upstream(k(i)-pl+2 : k(i)+1); %size(p,2) = 20   
            Tm = calculate_Tm_like_PrimerSelect(p);
            
            while Tm < 58.0 %Tm of primer needs to be above 58.0 C
                pl = pl + 1; %increase primer length by 1       
                p = upstream(k(i)-pl+2 : k(i)+1);
                Tm = calculate_Tm_like_PrimerSelect(p);
            end

            s = oligoprop(p); %calculate properties of primer
            if isempty(s.Hairpins)                        
                primer(1).name = [g.CDS(h).gene,'_KO_F1'];
                primer(1).seq = p;
                primer(1).Tm = Tm;                                
                primer(1).overhang = '';
                primer(1).length_binding = length(p);                        
                primer(1).length_overhang = 0;
                primer(1).length_total = primer(1).length_overhang + primer(1).length_binding;                
                primer(1).index_5p = regexp(g.Sequence,primer(1).seq);   
                primer(1).index_3p = primer(1).index_5p + primer(1).length_binding - 1;      
                orderfields(primer(1),{'name','seq','Tm','overhang','length_binding','length_overhang','length_total','index_5p','index_3p'});
                break
            else
            end      
        end
        primer(1) = correctPrimerBindingSite(g,h,primer(1));
        clear k i pl p Tm s
        %Note, not necessary to check whether KO F1 primer is larger than 60 nt since primer with 35x As or 35x Ts has a Tm of 58.1436 C              


        % ---------- Second primer (KO R1 primer) ----------
        display('----- Designing KO R1 primer -----')    
        pl = 20; %minimal primer length        
        p = g.Sequence(g.CDS(h).indices(1)-pl : g.CDS(h).indices(1)-1);   
        Tm = calculate_Tm_like_PrimerSelect(p);
        
        while Tm < 58.0
            pl = pl +1; %increase primer lenght by 1
            p = g.Sequence(g.CDS(h).indices(1)-pl : g.CDS(h).indices(1)-1);
            Tm = calculate_Tm_like_PrimerSelect(p);
        end

        %Try to include GC clamp if possible; Dirk 20140510    
        display(['Original KO R1 primer: ',seqrcomplement(p)])

        %Case #1: KO R1 primer has GC clamp
        if regexp(p(1:2), '^[gc]{2}')
            display('KO R1 primer has GC clamp.') % --> do nothing

        %Case #2: KO R1 primer has G/C at first position but A/T at second position    
        elseif regexp(p(1:2), '^[gc][at]')
            display('KO R1 primer has G/C at first position but A/T at second position.')                
            p_new = g.Sequence(g.CDS(h).indices(1)-(pl+1) : g.CDS(h).indices(1)-1); %increase primer lenght by 1        
            display(['KO R1 primer (+1): ',seqrcomplement(p_new)])
            
            if regexp(p_new(1), '^[gc]') % new KO R1 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %use old seqeunce; KO R1 primer has G/C at first position
            end

        %Case #3: KO R1 primer has A/T at first position but G/C at second last position    
        elseif regexp(p(1:2), '^[at][gc]')
            display('KO R1 primer has A/T at first position but G/C at second position.')               
            p_new = g.Sequence(g.CDS(h).indices(1)-(pl+2) : g.CDS(h).indices(1)-1); %increase primer lenght by 2        
            display(['KO R1 primer (+2): ',seqrcomplement(p_new)])
            
            if regexp(p_new(1:2), '^[gc]{2}') % new KO R1 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(1:2), '^[gc][at]') % new KO R1 primer has G/C at first position and A/T at second first position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(1:2), '^[at][gc]') % new KO R1 primer has A/T at first position and G/C at second first position
                %increase lenght of OLD KO R1 primer by 1 (this is same as decreasing length of NEW KO R1 primer by 1)
                p_new = g.Sequence(g.CDS(h).indices(1)-(pl+1) : g.CDS(h).indices(1)-1);                     
                display(['KO R1 primer (+1): ',seqrcomplement(p_new)])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else % new KO R1 primer has A/T at first and second position --> use old R1 primer
            end 

        %Case #4: KO R1 primer has A/T at first and second position
        elseif regexp(p(1:2), '^[at]{2}')        
            display('KO R1 primer has A/T at first and second position.')                            
            p_new = g.Sequence(g.CDS(h).indices(1)-(pl+2) : g.CDS(h).indices(1)-1); %increase primer lenght by 2        
            display(['KO R1 primer (+2): ',seqrcomplement(p_new)])
            
            if regexp(p_new(1:2), '^[gc]{2}') % new KO R1 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(1:2), '^[gc][at]') % new KO R1 primer has G/C at first position and A/T at second position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(1:2), '^[at][gc]') % new KO R1 primer has A/T at first position and G/C at second  position
                %increase lenght of OLD KO R1 primer by 1 (this is same as decreasing length of NEW KO R1 primer by 1)
                p_new = g.Sequence(g.CDS(h).indices(1)-(pl+1) : g.CDS(h).indices(1)-1); 
                display(['KO R1 primer (+1): ',seqrcomplement(p_new)])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %new KO R1 primer has A/T at first and second positions
                p = p_new; %NEW, Dirk 20151107; this extends the primer ending with e.g. 5'-AA to 5'-AAAA 
                Tm = calculate_Tm_like_PrimerSelect(p); %NEW, Dirk 20151107
            end                    
        end

        %R1 is a reverse primer --> take reverse complement of primer sequence        
        primer(2).name = [g.CDS(h).gene,'_KO_R1'];                
        primer(2).overhang = 'GAGGGTATTCTGGGCCTCCATGT';  %this is the reverse complement of the longer F primer
        %check length of KO R1 primer
        %primer cannot be longer than 37 nucleotide since the overhang is already 23 nucleotide and the maximal primer length is 60 nucleotide
        if size(p,2) > 60-length(primer(2).overhang)
            display(['Need to truncate KO R1 primer to ',num2str(60-length(primer(2).overhang)),' nucleotides (--> total length is 60 nucleotides).'])            
            %p = p(end-37+1:end); % --> truncate to 40 nucleotides           
            p = p(end-(60-length(primer(2).overhang))+1:end); % --> truncate to 40 nucleotides           
            Tm = calculate_Tm_like_PrimerSelect(p); %calculate the new Tm
        else
        end            
        primer(2).seq = [primer(2).overhang,seqrcomplement(p)]; %add 5' appendix to R1 primer        
        primer(2).Tm = Tm;        
        primer(2).length_binding = length(p);
        primer(2).length_overhang = length(primer(2).overhang);        
        primer(2).length_total = primer(2).length_binding + primer(2).length_overhang; %add appendix to primer length                     
        primer(2).index_3p = regexp(g.Sequence,p); %index of 5'
        primer(2).index_5p = primer(2).index_3p + primer(2).length_binding - 1; %index of 3'                        
        primer(2) = correctPrimerBindingSite(g,h,primer(2));
        clear k i pl p p_new Tm s

        
        % ---------- Third primer (KO F2 primer) ----------
        display('----- Designing KO F2 primer -----')
        pl = 20;
        p = g.Sequence(g.CDS(h).indices(end)+1 : g.CDS(h).indices(end)+pl);
        Tm = calculate_Tm_like_PrimerSelect(p);

        while Tm < 58.0
            pl = pl + 1; %increase primer lenght by 1
            p = g.Sequence(g.CDS(h).indices(end)+1 : g.CDS(h).indices(end)+pl);
            Tm = calculate_Tm_like_PrimerSelect(p);
        end

        %Try to include GC clamp in KO F2 primer if possible; Dirk 20140510    
        display(['Original KO F2 primer: ',p])

        %Case #1: KO F2 primer has already GC clamp
        if regexp(p(end-1:end), '^[gc]{2}')
            display('KO F2 primer has GC clamp.') % --> do nothing

        %Case #2: KO F2 primer has G/C at first position but A/T at second position        
        elseif regexp(p(end-1:end), '^[at][gc]')
            display('KO F2 primer has G/C at first position but A/T at second position.')               
            p_new = g.Sequence(g.CDS(h).indices(end)+1 : g.CDS(h).indices(end)+(pl+1)); %increase primer lenght by 1                
            display(['KO F2 primer (+1): ',p_new])
            
            if regexp(p_new(end), '^[gc]') % new KO F2 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %use old primer seqeunce; KO F2 primer has G/C at first position
            end

        %Case #3: KO F2 primer has A/T at first position but G/C at second last position        
        elseif regexp(p(end-1:end), '^[gc][at]')
            display('KO F2 Primer has A/T at first position but G/C at second position.')                               
            p_new = g.Sequence(g.CDS(h).indices(end)+1 : g.CDS(h).indices(end)+(pl+2)); %increase primer lenght by 2                
            display(['KO F2 primer (+2): ',p_new])
            
            if regexp(p_new(end-1:end), '^[gc]{2}') % new KO F2 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[at][gc]') % new KO F2 primer has G/C at first position and A/T at second first position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[gc][at]') % new KO F2 primer has A/T at first position and G/C at second first position            
                %increase lenght of OLD KO F2 primer by 1 (this is same as decreasing length of NEW KO F2 primer by 1)
                p_new = g.Sequence(g.CDS(h).indices(end)+1 : g.CDS(h).indices(end)+(pl+1));            
                display(['KO F2 primer (+1): ',p_new])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else % new KO F2 primer has A/T at first and second position --> use old primer
            end 

        %Case #4: KO F2 primer has A/T at first and second position
        elseif regexp(p(end-1:end), '^[at]{2}')        
            display('KO F2 primer has A/T at first and second position.')                        
            p_new = g.Sequence(g.CDS(h).indices(end)+1 : g.CDS(h).indices(end)+(pl+2)); %increase primer lenght by 2  
            display(['KO F2 primer (+2): ',p_new])
            
            if regexp(p_new(end-1:end), '^[gc]{2}') % new KO F2 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[at][gc]') % new KO F2 primer has G/C at first position and A/T at second position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[gc][at]') % new KO F2 primer has A/T at first position and G/C at second  position            
                %increase lenght of OLD KO F2 primer by 1 (this is same as decreasing length of NEW KO F2 primer by 1)
                p_new = g.Sequence(g.CDS(h).indices(end)+1 : g.CDS(h).indices(end)+(pl+1)); 
                display(['KO F2 primer (+1): ',p_new])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %new KO F2 primer has A/T at first and second positions
                p = p_new; %NEW, Dirk 20151107; this extends the primer ending with e.g. 5'-AA to 5'-AAAA 
                Tm = calculate_Tm_like_PrimerSelect(p); %NEW, Dirk 20151107
            end                    
        end

        primer(3).name = [g.CDS(h).gene,'_KO_F2'];                        
        primer(3).overhang = 'CGTATGTGAATGCTGGTCGCTATACTG'; %this is the reverse complememnt of the longer R primer
        %check length of F2 primer           
        %primer cannot be longer than 33 nucleotide since the overhang is already 27 nucleotide and the maximal primer length is 60 nucleotide
        if size(p,2) > 60-length(primer(3).overhang)
            display(['Need to truncate KO F2 primer to ',num2str(60-length(primer(3).overhang)),' nucleotides (--> total length is 60 nucleotides).'])
            %p = p(1:33); % --> truncate to 33 nucleotides           
            p = p(1:60-length(primer(3).overhang));
            Tm = calculate_Tm_like_PrimerSelect(p);
        else
        end                               
        primer(3).seq = [primer(3).overhang,p]; %add 5' appendix
        primer(3).Tm = Tm;
        primer(3).length_binding = length(p); %add appendix to primer length                
        primer(3).length_overhang = length(primer(3).overhang); %add appendix to primer length                
        primer(3).length_total = primer(3).length_overhang + primer(3).length_binding;
        primer(3).index_5p = regexp(g.Sequence,p);
        primer(3).index_3p = primer(3).index_5p + primer(3).length_binding - 1;                        
        primer(3) = correctPrimerBindingSite(g,h,primer(3));
        clear k i pl p p_new Tm s
    
        
        % ---------- Fourth primer (KO R2 primer) ----------
        display('----- Designing KO R2 primer -----')    
        k = regexp(downstream, '[c|g]{2}'); %find all GC clamps
        k = sort(k,'ascend'); %start with the primer closest to stop codon and move further downstream from there    

        for i = 1:size(k,2)
            pl = 20;
            p = downstream(k(i) : k(i)+pl-1);                
            Tm = calculate_Tm_like_PrimerSelect(p); 
            
            while Tm < 58.0        
                pl = pl +1;                        
                p = downstream(k(i) : k(i)+pl-1);
                Tm = calculate_Tm_like_PrimerSelect(p);                
            end    

            s = oligoprop(p);
            if isempty(s.Hairpins)        
                primer(4).name = [g.CDS(h).gene,'_KO_R2'];                
                primer(4).seq = seqrcomplement(p);
                primer(4).Tm = Tm;
                primer(4).overhang = '';
                primer(4).length_binding = length(p);            
                primer(4).length_overhang = length(primer(4).overhang);
                primer(4).length_total = primer(4).length_overhang + primer(4).length_binding;            
                primer(4).index_3p = regexp(g.Sequence,seqrcomplement(primer(4).seq)); %index of 5'
                primer(4).index_5p = primer(4).index_3p + primer(4).length_binding - 1;                                
                break
            else
            end      
        end
        primer(4) = correctPrimerBindingSite(g,h,primer(4));
        clear k i pl p Tm s

        
        % ---------- Calculate length of PCR product H1_KO ----------        
        H1_length = primer(2).index_5p + 1 + primer(2).length_overhang - primer(1).index_5p;
        display(['Length of PCR product H1_KO is ',num2str(H1_length),' bp.'])
                
         % ---------- Calculate length of PCR product H2_KO ----------        
        H2_length = primer(4).index_5p - primer(3).index_5p + 1 + primer(3).length_overhang;
        display(['Length of PCR product H2_KO is ',num2str(H2_length),' bp.'])
                
        %write PCR product lenght of H1_KO and H2_KO to a text file
        filename_out3 = [datestr(date,'yyyymmdd'),'_length_of_H1_KO_and_H2_KO_PCR_products_for_',g.CDS(h).gene,'.txt']; display(filename_out3);
        fid = fopen(filename_out3, 'w','n','US-ASCII'); %write mode
        fprintf(fid, '%s\t%s\t%s\t\n','Gene','Name','Length (bp)');        
        fprintf(fid, '%s\t%s\t%d\t\n',g.CDS(h).gene,'H1_KO', H1_length);
        fprintf(fid, '%s\t%s\t%d\t\n',g.CDS(h).gene,'H2_KO', H2_length);
        fclose (fid);
        
        
        % ---------- Fifth primer: forward primer for colony PCR (KO Fc primer) ----------
        display('----- Designing KO Fc primer -----')           
        k = strfind(upstream,primer(1).seq);
        upstream2 = upstream(1:k-1); %region upstream of the F1 primer;
        k = regexp(upstream2, '[c|g]{2}'); %find all GC clamps in the region upstream of the KO F1 primer
        k = sort(k,'descend'); %start with the primer closest to the KO F1 primer and move further upstream from there

        for i = 1:size(k,2) %loop over all potential primers    
            pl = 20; %minimal primer length
            p = upstream2(k(i)-pl+2 : k(i)+1); %size(p,2) = 20   
            Tm = calculate_Tm_like_PrimerSelect(p);
            
            while Tm < 58.0 %Tm of primer needs to be above 58.0 C
                pl = pl +1; %increase primer length by 1       
                p = upstream2(k(i)-pl+2 : k(i)+1);
                Tm = calculate_Tm_like_PrimerSelect(p);
            end

            s = oligoprop(p); %calculate properties of primer
            if isempty(s.Hairpins)                        
                primer(5).name = [g.CDS(h).gene,'_KO_Fc'];
                primer(5).seq = p;
                primer(5).Tm = Tm;
                primer(5).overhang = '';
                primer(5).length_binding = length(p);    
                primer(5).length_overhang = length(primer(5).overhang);    
                primer(5).length_total = primer(5).length_overhang + primer(5).length_binding;                
                %primer(5).index_3p = regexp(g.Sequence,p); %20151108
                primer(5).index_5p = regexp(g.Sequence,p);
                %primer(5).index_5p = primer(5).index_3p + primer(5).length_binding - 1; %20151108                               
                primer(5).index_3p = primer(5).index_5p + primer(5).length_binding - 1;                                
                break
            else
            end      
        end
        primer(5) = correctPrimerBindingSite(g,h,primer(5));
        clear k i pl p Tm s

                                                                                                                                                                              
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % --- DNA antisense strand (gene has reverse orientation) ---
    elseif g.CDS(h).indices(end) < g.CDS(h).indices(1)         
        display([g.CDS(h).gene, ' gene has reverse orientation.'])    

        ORF = seqrcomplement(g.Sequence(g.CDS(h).indices(end)+3 : g.CDS(h).indices(1)));    
        start_codon = seqrcomplement(g.Sequence(g.CDS(h).indices(1)-2 : g.CDS(h).indices(1))); 
        display(['Start codon: ',upper(start_codon)])
        stop_codon = seqrcomplement(g.Sequence(g.CDS(h).indices(end) : g.CDS(h).indices(end)+3-1));
        display(['Stop codon: ',upper(stop_codon)])

        try
        %500-bp search window for F1 primer; search window starts 280 bp upstream of stop codon                
        upstream = seqrcomplement(g.Sequence(g.CDS(h).indices(1)+281 : g.CDS(h).indices(1)+780)); display(upstream);
        %500-bp search window for R2 primer; search window starts 280 bp upstream of stop codon
        downstream = seqrcomplement(g.Sequence(g.CDS(h).indices(end)-780 : g.CDS(h).indices(end)-281)); display(downstream);
        catch ME
            display(['Upstream or downstream sequence(s) are not long enough. Skipped primer design for gene ',g.CDS(h).gene,'.'])
            primer = addEmptyPrimerEntries(g.CDS(h).gene);
            return
        end
        
        % ---------- First primer (KO F1 primer) ----------
        display('----- Designing KO F1 primer -----')    
        k = regexp(upstream, '[c|g]{2}'); %find all GC clamps
        k = sort(k,'descend'); %start with the primer closest to stop codon and move further upstream from there
        
        for i = 1:size(k,2) %loop over all potential primers    
            pl = 20; %minimal primer length
            p = upstream(k(i)-pl+2 : k(i)+1); %size(p,2) = 20       
            Tm = calculate_Tm_like_PrimerSelect(p);
        
            while Tm < 58.0 %Tm of primer needs to be above 58.0 C
                pl = pl +1; %increase primer length by 1       
                p = upstream(k(i)-pl+2 : k(i)+1);
                Tm = calculate_Tm_like_PrimerSelect(p);
            end

            s = oligoprop(p); %calculate properties of primer
            if isempty(s.Hairpins)                        
                primer(1).name = [g.CDS(h).gene,'_KO_F1'];
                primer(1).seq = p;
                primer(1).Tm = Tm;
                primer(1).overhang = '';
                primer(1).length_binding = length(p);                        
                primer(1).length_overhang = 0;
                primer(1).length_total = primer(1).length_overhang + primer(1).length_binding;
                primer(1).index_3p = regexp(g.Sequence,seqrcomplement(primer(1).seq)); 
                primer(1).index_5p = primer(1).index_3p + primer(1).length_binding - 1;
                orderfields(primer(1),{'name','seq','Tm','overhang','length_binding','length_overhang','length_total','index_5p','index_3p'});                
                break
            else
            end      
        end
        primer(1) = correctPrimerBindingSite(g,h,primer(1));
        clear k i pl p Tm s
        %Note, not necessary to check whether KO F1 primer is larger than 60 nt since primer with 35x As or 35x Ts has a Tm of 58.1436 C  

        % ---------- Second primer (KO R1 primer) ----------
        display('----- Designing KO R1 primer -----')                    
        pl = 20; %minimal primer length            
        p = g.Sequence(g.CDS(h).indices(1)+1 : g.CDS(h).indices(1)+pl);
        Tm = calculate_Tm_like_PrimerSelect(p);
        
        while Tm < 58.0
            pl = pl + 1; %increase primer lenght by 1                    
            p = g.Sequence(g.CDS(h).indices(1)+1 : g.CDS(h).indices(1)+pl);
            Tm = calculate_Tm_like_PrimerSelect(p);
        end

        %Try to include GC clamp if possible; Dirk 20140513    
        display(['Original KO R1 primer: ',p])
        %no need to take the reverse complement here since the gene is already in reverse orientation

        %Case #1: KO R1 primer has GC clamp
        if regexp(p(end-1:end), '^[gc]{2}')
            display('KO R1 primer has GC clamp.') % --> do nothing

        %Case #2: KO R1 primer has G/C at first position but A/T at second position    
        elseif regexp(p(end-1:end), '^[at][gc]')
            display('KO R1 primer has G/C at first position but A/T at second position.')                                    
            p_new = g.Sequence(g.CDS(h).indices(1)+1 : g.CDS(h).indices(1)+(pl+1));%increase primer lenght by 1
            display(['KO R1 primer (+1): ',p_new])
            
            if regexp(p_new(end), '^[gc]') % new KO R1 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %use old seqeunce; KO R1 primer has G/C at first position
            end

        %Case #3: KO R1 primer has A/T at first position but G/C at second last position    
        elseif regexp(p(end-1:end), '^[gc][at]')
            display('KO R1 primer has A/T at first position but G/C at second position.')                                       
            p_new = g.Sequence(g.CDS(h).indices(1)+1 : g.CDS(h).indices(1)+(pl+2)); %increase primer lenght by 2
            display(['KO R1 primer (+2): ',p_new])
            
            if regexp(p_new(end-1:end), '^[gc]{2}') % new KO R1 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[at][gc]') % new KO R1 primer has G/C at first position and A/T at second first position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[gc][at]') % new KO R1 primer has A/T at first position and G/C at second first position
                %increase lenght of OLD KO R1 primer by 1 (this is same as decreasing length of NEW KO R1 primer by 1)                
                p_new = g.Sequence(g.CDS(h).indices(1)+1 : g.CDS(h).indices(1)+(pl+1)); 
                display(['KO R1 primer (+1): ',p_new])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else % new R1 primer has A/T at first and second position --> use old KO R1 primer
            end 

        %Case #4: KO R1 primer has A/T at first and second position
        elseif regexp(p(end-1:end), '^[at]{2}')        
            display('KO R1 primer has A/T at first and second position.')                                    
            p_new = g.Sequence(g.CDS(h).indices(1)+1 : g.CDS(h).indices(1)+(pl+2)); %increase primer lenght by 2
            display(['KO R1 primer (+2): ',p_new])
            
            if regexp(p_new(end-1:end), '^[gc]{2}') % new KO R1 primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[at][gc]') % new KO R1 primer has G/C at first position and A/T at second position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[gc][at]') % new KO R1 primer has A/T at first position and G/C at second  position
                %increase lenght of OLD KO R1 primer by 1 (this is same as decreasing length of NEW KO R1 primer by 1)                
                p_new = g.Sequence(g.CDS(h).indices(1)+1 : g.CDS(h).indices(1)+(pl+1));
                display(['KO R1 primer (+1): ',p_new])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %new KO R1 primer has A/T at first and second positions
                p = p_new; %NEW, Dirk 20151107; this extends the primer ending with e.g. 5'-AA to 5'-AAAA 
                Tm = calculate_Tm_like_PrimerSelect(p); %NEW, Dirk 20151107
            end                    
        end

        %check length of KO R1 primer
        %primer cannot be longer than 40 nucleotide since the overhang is already 20 nucleotide and the maximal primer length is 60 nucleotide
        if size(p,2) > (60 - 23) 
            display('Need to truncate KO R1 primer to 37 nucleotides (--> total length is 60 nucleotides).')
            p = p(1:37); % --> truncate to 37 nucleotides           
            Tm = calculate_Tm_like_PrimerSelect(p);
        else
        end    
        
        primer(2).name = [g.CDS(h).gene,'_KO_R1'];
        %no need to take the reverse complement here since the gene is already in reverse orientation        
        primer(2).overhang = 'GAGGGTATTCTGGGCCTCCATGT';  %this is the reverse complement of the longer F primer
        primer(2).seq = [primer(2).overhang,p]; %add 5' appendix to R1 primer
        primer(2).Tm = Tm;        
        primer(2).length_binding = length(p);
        primer(2).length_overhang = length(primer(2).overhang);        
        primer(2).length_total = primer(2).length_binding + primer(2).length_overhang; %add appendix to primer length             
        primer(2).index_5p = regexp(g.Sequence,p); %index of 5'
        primer(2).index_3p = primer(2).index_5p + primer(2).length_binding - 1; %index of 3'
        primer(2) = correctPrimerBindingSite(g,h,primer(2));
        clear k i pl p p_new Tm s
                               
        
        % ---------- Third primer (KO F2 primer) ----------
        display('----- Designing KO F2 primer -----')
        pl = 20;
        p = seqrcomplement(g.Sequence(g.CDS(h).indices(end)-pl : g.CDS(h).indices(end)-1)); %this is now a reverse primer
        Tm = calculate_Tm_like_PrimerSelect(p);
        
        while Tm < 58.0
            pl = pl +1; %increase primer lenght by 1
            p = seqrcomplement(g.Sequence(g.CDS(h).indices(end)-pl : g.CDS(h).indices(end)-1));
            Tm = calculate_Tm_like_PrimerSelect(p);
        end

        %Try to include GC clamp in KO F2 primer if possible; Dirk 20140513    
        display(['Original KO F2 primer: ',p])

        %Case #1: KO F2 primer has already GC clamp
        if regexp(p(end-1:end), '^[gc]{2}')
            display('KO F2 primer has GC clamp.') % --> do nothing

        %Case #2: KO F2 primer has G/C at first position but A/T at second position            
        elseif regexp(p(end-1:end), '^[at][gc]')
            display('KO F2 primer has G/C at first position but A/T at second position.')                   
            p_new = seqrcomplement(g.Sequence(g.CDS(h).indices(end)-(pl+1) : g.CDS(h).indices(end)-1)); %increase primer lenght by 1
            display(['KO F2 primer (+1): ',p_new])
            
            if regexp(p_new(end), '^[gc]') % new primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %use old primer seqeunce; KO F2 primer has G/C at first position
            end

        %Case #3: KO F2 primer has A/T at first position but G/C at second last position            
        elseif regexp(p(end-1:end), '^[gc][at]')
            display('KO F2 primer has A/T at first position but G/C at second position.')                               
            p_new = seqrcomplement(g.Sequence(g.CDS(h).indices(end)-(pl+2) : g.CDS(h).indices(end)-1)); %increase primer lenght by 2
            display(['KO F2 primer (+2): ',p_new])
            
            if regexp(p_new(end-1:end), '^[gc]{2}') % new primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[at][gc]') % new primer has G/C at first position and A/T at second first position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[gc][at]') % new primer has A/T at first position and G/C at second first position            
                %increase lenght of OLD KO F2 primer by 1 (this is same as decreasing length of NEW KO F2 primer by 1)
                p_new = seqrcomplement(g.Sequence(g.CDS(h).indices(end)-(pl+1) : g.CDS(h).indices(end)-1)); 
                display(['KO F2 primer (+1): ',p_new])
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            else %new KO F2 primer has A/T at first and second position --> use old F2 primer
            end         

        %Case #4: KO F2 primer has A/T at first and second position
        elseif regexp(p(end-1:end), '^[at]{2}')        
            display('KO F2 primer has A/T at first and second position.')                                
            p_new = seqrcomplement(g.Sequence(g.CDS(h).indices(end)-(pl+2) : g.CDS(h).indices(end)-1)); %increase primer lenght by 2
            display(['KO F2 primer (+2): ',p_new])
            
            if regexp(p_new(end-1:end), '^[gc]{2}') %new primer has now GC clamp
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[at][gc]') %new primer has G/C at first position and A/T at second position
                p = p_new;
                Tm = calculate_Tm_like_PrimerSelect(p);
            elseif regexp(p_new(end-1:end), '^[gc][at]') %new primer has A/T at first position and G/C at second  position            
            %increase lenght of OLD KO F2 primer by 1 (this is same as decreasing length of NEW KO F2 primer by 1) 
            p_new = seqrcomplement(g.Sequence(g.CDS(h).indices(end)-(pl+1) : g.CDS(h).indices(end)-1));            
            display(['KO F2 primer (+1): ',p_new])
            p = p_new;
            Tm = calculate_Tm_like_PrimerSelect(p);
            else %new KO F2 primer has A/T at first and second positions
                p = p_new; %NEW, Dirk 20151107; this extends the primer ending with e.g. 5'-AA to 5'-AAAA 
                Tm = calculate_Tm_like_PrimerSelect(p); %NEW, Dirk 20151107
            end                    
        end

        primer(3).name = [g.CDS(h).gene,'_KO_F2'];        
        primer(3).overhang = 'CGTATGTGAATGCTGGTCGCTATACTG'; %this is the reverse complememnt of the longer R primer
        %check length of KO F2 primer
        %primer cannot be longer than 33 nucleotide since the overhang is already 27 nucleotide and the maximal primer length is 60 nucleotide
        if size(p,2) > 60-length(primer(3).overhang) 
            display(['Need to truncate KO F2 primer to ',num2str(60-length(primer(3).overhang)),' nucleotides (--> total length is 60 nucleotides).'])
            %p = p(1:33); % --> truncate to 33 nucleotides           
            p = p(1:60-length(primer(3).overhang));
            Tm = calculate_Tm_like_PrimerSelect(p);
        else
        end                     
        primer(3).seq = [primer(3).overhang,p]; %add 5' appendix
        primer(3).Tm = Tm;
        primer(3).length_binding = length(p); %add appendix to primer length                
        primer(3).length_overhang = length(primer(3).overhang); %add appendix to primer length                
        primer(3).length_total = primer(3).length_overhang + primer(3).length_binding;
        primer(3).index_3p = regexp(g.Sequence,seqrcomplement(p));
        primer(3).index_5p = primer(3).index_3p + primer(3).length_binding - 1;
        primer(3) = correctPrimerBindingSite(g,h,primer(3));
        clear k i pl p p_new Tm s

        
        % ---------- Fourth primer (KO R2 primer) ----------
        display('----- Designing KO R2 primer -----')    
        k = regexp(downstream, '[c|g]{2}'); %find all GC clamps
        k = sort(k,'ascend'); %start with the primer closest to stop codon and move further downstream from there    
        
        for i = 1:size(k,2)   
            pl = 20;
            p = downstream(k(i) : k(i)+pl-1);    
            Tm = calculate_Tm_like_PrimerSelect(p);    
        
            while Tm < 58.0
                pl = pl+1;                   
                p = downstream(k(i) : k(i)+pl-1);
                Tm = calculate_Tm_like_PrimerSelect(p);                        
            end    

            s = oligoprop(p);
            if isempty(s.Hairpins)                        
                primer(4).name = [g.CDS(h).gene,'_KO_R2'];
                primer(4).seq = seqrcomplement(p);
                primer(4).Tm = Tm;
                primer(4).overhang = '';
                primer(4).length_binding = length(p);            
                primer(4).length_overhang = length(primer(4).overhang);
                primer(4).length_total = primer(4).length_overhang + primer(4).length_binding;            
                primer(4).index_5p = regexp(g.Sequence,primer(4).seq); %index of 5'
                primer(4).index_3p = primer(4).index_5p + primer(4).length_binding - 1;
                break
            else
            end      
        end
        primer(4) = correctPrimerBindingSite(g,h,primer(4));
        clear k i pl p Tm s
        
        
        
         % ---------- Calculate length of PCR product H1_KO ----------        
        H1_length = primer(1).index_5p - primer(2).index_5p + 1 + primer(2).length_overhang;
        display(['Length of PCR product H1_KO is ',num2str(H1_length),' bp.'])
                
         % ---------- Calculate length of PCR product H2_KO ----------        
        H2_length = primer(3).index_5p - primer(4).index_5p + 1 + primer(3).length_overhang;
        display(['Length of PCR product H2_KO is ',num2str(H2_length),' bp.'])
                        
        %write PCR product lenght of H1_KO and H2_KO to a text file
        filename_out3 = [datestr(date,'yyyymmdd'),'_length_of_H1_KO_and_H2_KO_PCR_products_for_',g.CDS(h).gene,'.txt']; display(filename_out3);
        fid = fopen(filename_out3, 'w','n','US-ASCII'); %write mode
        fprintf(fid, '%s\t%s\t%s\t\n','Gene','Name','Length (bp)');        
        fprintf(fid, '%s\t%s\t%d\t\n',g.CDS(h).gene,'H1_KO', H1_length);
        fprintf(fid, '%s\t%s\t%d\t\n',g.CDS(h).gene,'H2_KO', H2_length);
        fclose (fid);
        

        % ---------- Fifth primer: forward primer for colony PCR (KO Fc primer) ----------
        display('----- Designing KO Fc primer -----')           
        k = strfind(upstream,primer(1).seq);
        upstream2 = upstream(1:k-1); %region upstream of the KO F1 primer; k(i) comes from the calculation of the KO F1 primer        
        k = regexp(upstream2, '[c|g]{2}'); %find all GC clamps in the region upstream of the KO F1 primer
        k = sort(k,'descend'); %start with the primer closest to the KO F1 primer and move further upstream from there
        
        for i = 1:size(k,2) %loop over all potential primers    
            pl = 20; %minimal primer length
            p = upstream2(k(i)-pl+2 : k(i)+1); %size(p,2) = 20   
            Tm = calculate_Tm_like_PrimerSelect(p);
        
            while Tm < 58.0 %Tm of primer needs to be above 58.0 C
                pl = pl +1; %increase primer length by 1       
                p = upstream2(k(i)-pl+2 : k(i)+1);
                Tm = calculate_Tm_like_PrimerSelect(p);
            end

            s = oligoprop(p); %calculate properties of primer
            if isempty(s.Hairpins)                        
                primer(5).name = [g.CDS(h).gene,'_KO_Fc'];
                primer(5).seq = p;
                primer(5).Tm = Tm;
                primer(5).overhang = '';
                primer(5).length_binding = length(p);    
                primer(5).length_overhang = length(primer(5).overhang);    
                primer(5).length_total = primer(5).length_overhang + primer(5).length_binding;                
                primer(5).index_3p = regexp(g.Sequence,seqrcomplement(p));
                primer(5).index_5p = primer(5).index_3p + primer(5).length_binding -1;                                
                break
            else
            end      
        end
        primer(5) = correctPrimerBindingSite(g,h,primer(5));
        clear k i pl p Tm s
                      
    else
    end
    %}    
                   
    
    
    %{
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------ METHOD #4. Internal fusions ----- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'internal'
        
    %I haven't started working on this    
        
    %First determine whether gene is on the upper (sense) or lower (anti-sense) strand
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % --- DNA sense strand (gene has normal orientation) ---
    if g.CDS(h).indices(end) > g.CDS(h).indices(1)         
        display([g.CDS(h).gene, ' gene has normal orientation.'])    
        ORF = g.Sequence(g.CDS(h).indices(1) : g.CDS(h).indices(end)-3); %exclude stop codon
        stop_codon = g.Sequence(g.CDS(h).indices(end)-2 : g.CDS(h).indices(end)); display(['Stop codon: ',upper(stop_codon)])
        
        try
        %500-bp search window for F1 primer; search window starts 280 bp upstream of stop codon
        upstream = g.Sequence(g.CDS(h).indices(end)-782 : g.CDS(h).indices(end)-283); %Dirk 20140510        
        %500-bp search window for R2 primer; search window starts 280 bp upstream of stop codon
        downstream = g.Sequence(g.CDS(h).indices(end)+281 : g.CDS(h).indices(end)+780); %Dirk 20140510
        catch ME
            display(['Upstream or downstream sequence(s) are not long enough. Skipped primer design for gene ',g.CDS(h).gene,'.'])
            primer = addEmptyPrimerEntries(g.CDS(h).gene);
            return
        end
        
        %...

    else
        
        %...
    end
    %}
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ------ METHOD #... ----- %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    otherwise
    display('Method is not implemented yet.')        

    
    
    
    
    
    
end


%Display the designed primers
display(' ')
display(primer(1))
display(primer(2))
display(primer(3))
display(primer(4))
display(primer(5))


