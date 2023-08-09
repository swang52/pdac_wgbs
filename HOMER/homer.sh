# Description
#HOMER known transcription factor motif analyses for all DMRs

source ~/.bash_profile
cd HOMER

# Mouse DMRs
findMotifsGenome.pl MPN_DMRs_hyper.bed mm9 MPN_hyper_output -bg MPN_background.bed -size given -cpg
findMotifsGenome.pl MPN_DMRs_hypo.bed mm9 MPN_hypo_output -bg MPN_background.bed -size given -cpg
findMotifsGenome.pl TM_DMRs_hyper.bed mm9 TM_hyper_output -bg TM_background.bed -size given -cpg
findMotifsGenome.pl TM_DMRs_hypo.bed mm9 TM_hypo_output -bg TM_background.bed -size given -cpg

# Human DMRs
findMotifsGenome.pl NT_DMRs_hyper.bed hg38 MPN_hyper_output -bg NT_background.bed -size given -cpg
findMotifsGenome.pl NT_DMRs_hypo.bed hg38 MPN_hypo_output -bg NT_background.bed -size given -cpg
findMotifsGenome.pl EL_DMRs_hyper.bed hg38 TM_hyper_output -bg EL_background.bed -size given -cpg
findMotifsGenome.pl EL_DMRs_hypo.bed hg38 TM_hypo_output -bg EL_background.bed -size given -cpg
findMotifsGenome.pl sub_DMRs_hyper.bed hg38 sub_hyper_output -bg sub_background.bed -size given -cpg
findMotifsGenome.pl sub_DMRs_hypo.bed hg38 sub_hypo_output -bg sub_background.bed -size given -cpg