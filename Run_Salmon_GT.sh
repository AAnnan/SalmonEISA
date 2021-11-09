#! /bin/bash

## Resource Allocation
#SBATCH --time=1-00:00:00
#SBATCH --partition=local
#SBATCH --mem=64G
#SBATCH --cpus-per-task=30

#SBATCH --mail-user=ahrmad.annan@students.unibe.ch
#SBATCH --mail-type=end,fail
#SBATCH --job-name="Salmon"

source ${CONDA_ACTIVATE} salmon

cat c_elegans.PRJNA13758.WS279.mRNA_transcripts.fa CE_genes_seq.fa c_elegans.PRJNA13758.WS279.genomic.fa > CE_gentrome.fa

salmon index -p 30 -t CE_gentrome.fa -i EISA_index --decoys decoys.txt -k 19
#the k size selected here will act as the minimum acceptable length for a valid match.

########################################################################################################################
########################################################################################################################
data=/scratch/aannan/rnaseq/data

#366 b 2020
fa366b1_L1_R1=${data}/366b1/366B1_L1_R1_001_F6TKUgD2sWF3.fastq.gz
fa366b1_L2_R1=${data}/366b1/366B1_L2_R1_001_Aekf4XsT6VYC.fastq.gz
fa366b2_L1_R1=${data}/366b2/366B2_L1_R1_001_VukPzfB7UznO.fastq.gz
fa366b2_L2_R1=${data}/366b2/366B2_L2_R1_001_6eYI0PJPBw8Q.fastq.gz
fa366b3_L1_R1=${data}/366b3/366B3_L1_R1_001_9YNWmzpY47EG.fastq.gz
fa366b3_L2_R1=${data}/366b3/366B3_L2_R1_001_v1BuQK0d8xsQ.fastq.gz
fa366b4_L1_R1=${data}/366b4/366B4_L1_R1_001_oOmUjBH4ogzp.fastq.gz
fa366b4_L2_R1=${data}/366b4/366B4_L2_R1_001_KluTH7GqLaqT.fastq.gz

#382 b 2020
fa382b1_L1_R1=${data}/382b1/382B1_L1_R1_001_bBV52jLO2VX2.fastq.gz
fa382b1_L2_R1=${data}/382b1/382B1_L2_R1_001_wMojKBI3VV5t.fastq.gz
fa382b2_L1_R1=${data}/382b2/382B2_L1_R1_001_tt73nYLLPlmX.fastq.gz
fa382b2_L2_R1=${data}/382b2/382B2_L2_R1_001_u3AK2wvq6gvq.fastq.gz
fa382b3_L1_R1=${data}/382b3/382B3_L1_R1_001_yMMdM9dtpDE9.fastq.gz
fa382b3_L2_R1=${data}/382b3/382B3_L2_R1_001_JZWqrZz0eUf3.fastq.gz
fa382b4_L1_R1=${data}/382b4/382B4_L1_R1_001_Y6Xklx9GlZkI.fastq.gz
fa382b4_L2_R1=${data}/382b4/382B4_L2_R1_001_EGVBFP42rXm4.fastq.gz

#775 b 2020
fa775b1_L1_R1=${data}/775b1/775B1_L1_R1_001_sqW7EiJxR1jQ.fastq.gz
fa775b1_L2_R1=${data}/775b1/775B1_L2_R1_001_aY1ESmg57cND.fastq.gz
fa775b2_L1_R1=${data}/775b2/775B2_L1_R1_001_7cNYJTOgtfMO.fastq.gz
fa775b2_L2_R1=${data}/775b2/775B2_L2_R1_001_KpQ3uyTDcRn7.fastq.gz
fa775b3_L1_R1=${data}/775b3/775B3_L1_R1_001_MhTCU4uNQ8pN.fastq.gz
fa775b3_L2_R1=${data}/775b3/775B3_L2_R1_001_UYIRFgF1LfSl.fastq.gz
fa775b4_L1_R1=${data}/775b4/775B4_L1_R1_001_VqnDEIFMH0zU.fastq.gz
fa775b4_L2_R1=${data}/775b4/775B4_L2_R1_001_Vai70bmYnzo1.fastq.gz

#784 b 2020
fa784b1_L1_R1=${data}/784b1/784B1_L1_R1_001_s2C4bcxqbHMV.fastq.gz
fa784b1_L2_R1=${data}/784b1/784B1_L2_R1_001_oQla0AFVtg0m.fastq.gz
fa784b2_L1_R1=${data}/784b2/784B2_L1_R1_001_zjDxs06VWdYb.fastq.gz
fa784b2_L2_R1=${data}/784b2/784B2_L2_R1_001_CNh0JBmkbTz5.fastq.gz
fa784b4_L1_R1=${data}/784b4/784B4_L1_R1_001_LIP6tzX8axBm.fastq.gz
fa784b4_L2_R1=${data}/784b4/784B4_L2_R1_001_Uy5Ue5wY2n5n.fastq.gz
fa784b5_L1_R1=${data}/784b5/784B5_L1_R1_001_ffvRNNk1nHyM.fastq.gz
fa784b5_L2_R1=${data}/784b5/784B5_L2_R1_001_OBLGwAYsTO12.fastq.gz

#366 C 2021
fa366C1_L1_R1=${data}/366C_1/366C1_L1_R1_001_yYNqwishzWTA.fastq.gz
fa366C1_L1_R2=${data}/366C_1/366C1_L1_R2_001_ERQ6rABXFI7g.fastq.gz
fa366C1_L2_R1=${data}/366C_1/366C1_L2_R1_001_fzQXUDsOOODQ.fastq.gz
fa366C1_L2_R2=${data}/366C_1/366C1_L2_R2_001_YuRKEPO8cjem.fastq.gz
fa366C2_L1_R1=${data}/366_C2/366C2_L1_R1_001_pl8yzAZgjsvK.fastq.gz
fa366C2_L1_R2=${data}/366_C2/366C2_L1_R2_001_zhd5zGlpY8GX.fastq.gz
fa366C2_L2_R1=${data}/366_C2/366C2_L2_R1_001_3k6qtwyWLYjh.fastq.gz
fa366C2_L2_R2=${data}/366_C2/366C2_L2_R2_001_zlZk43H3Ghbd.fastq.gz
fa366C3_L1_R1=${data}/366_C3/366C3_L1_R1_001_y5hWcd7HkCrr.fastq.gz
fa366C3_L1_R2=${data}/366_C3/366C3_L1_R2_001_Z1XX40dvM868.fastq.gz
fa366C3_L2_R1=${data}/366_C3/366C3_L2_R1_001_KQg5NmSZshKq.fastq.gz
fa366C3_L2_R2=${data}/366_C3/366C3_L2_R2_001_DdcIruIUUJss.fastq.gz

#366 D 2021
fa366D1_L1_R1=${data}/366D_1/366D1_L1_R1_001_UuSFWY4G2FHZ.fastq.gz
fa366D1_L1_R2=${data}/366D_1/366D1_L1_R2_001_EhUNkdmr5gT7.fastq.gz
fa366D1_L2_R1=${data}/366D_1/366D1_L2_R1_001_8LNvGvyniyI7.fastq.gz
fa366D1_L2_R2=${data}/366D_1/366D1_L2_R2_001_Eqtj8kyBfb67.fastq.gz
fa366D2_L1_R1=${data}/366_D2/366D2_L1_R1_001_rOaWPSFG1otu.fastq.gz
fa366D2_L1_R2=${data}/366_D2/366D2_L1_R2_001_WCBjpodAaEJC.fastq.gz
fa366D2_L2_R1=${data}/366_D2/366D2_L2_R1_001_wu1s11EhFqfW.fastq.gz
fa366D2_L2_R2=${data}/366_D2/366D2_L2_R2_001_E4ttDA9R76zx.fastq.gz
fa366D3_L1_R1=${data}/366_D3/366D3_L1_R1_001_1E8qNYEyhaL4.fastq.gz
fa366D3_L1_R2=${data}/366_D3/366D3_L1_R2_001_4B6ETpDziqhX.fastq.gz
fa366D3_L2_R1=${data}/366_D3/366D3_L2_R1_001_kTZRIF0AqgwB.fastq.gz
fa366D3_L2_R2=${data}/366_D3/366D3_L2_R2_001_trzJtwANJTZ5.fastq.gz

#821 D 2021
fa821D1_L1_R1=${data}/821D_1/821D1_L1_R1_001_1O6TNic16ZsY.fastq.gz
fa821D1_L1_R2=${data}/821D_1/821D1_L1_R2_001_sdKot4uUzwt1.fastq.gz
fa821D1_L2_R1=${data}/821D_1/821D1_L2_R1_001_JRa1yULZNytF.fastq.gz
fa821D1_L2_R2=${data}/821D_1/821D1_L2_R2_001_GuQ16vl51mht.fastq.gz
fa821D2_L1_R1=${data}/821_D2/821D2_L1_R1_001_wUmamVGIE6v7.fastq.gz
fa821D2_L1_R2=${data}/821_D2/821D2_L1_R2_001_I5s0LYX72nG9.fastq.gz
fa821D2_L2_R1=${data}/821_D2/821D2_L2_R1_001_FEq0gMefVBAe.fastq.gz
fa821D2_L2_R2=${data}/821_D2/821D2_L2_R2_001_7UJ2UtoZjjyl.fastq.gz
fa821D3_L1_R1=${data}/821_D3/821D3_L1_R1_001_s3W6FhuNoZDl.fastq.gz
fa821D3_L1_R2=${data}/821_D3/821D3_L1_R2_001_u3N20gW4yBO2.fastq.gz
fa821D3_L2_R1=${data}/821_D3/821D3_L2_R1_001_p0UnXXbQs8tR.fastq.gz
fa821D3_L2_R2=${data}/821_D3/821D3_L2_R2_001_G3CSmZLXult2.fastq.gz

#822 A 2021
fa822A1_L1_R1=${data}/822_A1/822A1_L1_R1_001_7XaIcsVhP3pF.fastq.gz
fa822A1_L1_R2=${data}/822_A1/822A1_L1_R2_001_CQ6x14edcBBJ.fastq.gz
fa822A1_L2_R1=${data}/822_A1/822A1_L2_R1_001_2P63JiFjpFGk.fastq.gz
fa822A1_L2_R2=${data}/822_A1/822A1_L2_R2_001_aLzy7TuURwXW.fastq.gz
fa822A2_L1_R1=${data}/822_A2/822A2_L1_R1_001_hqTt43owtXv4.fastq.gz
fa822A2_L1_R2=${data}/822_A2/822A2_L1_R2_001_Q8nwCRQXdAFN.fastq.gz
fa822A2_L2_R1=${data}/822_A2/822A2_L2_R1_001_5i7lnLNFgIQo.fastq.gz
fa822A2_L2_R2=${data}/822_A2/822A2_L2_R2_001_08bLuGPD09DT.fastq.gz
fa822A3_L1_R1=${data}/822_A3/822A3_L1_R1_001_HbTjF2Obg3YH.fastq.gz
fa822A3_L1_R2=${data}/822_A3/822A3_L1_R2_001_KVGuzXrTs8A3.fastq.gz
fa822A3_L2_R1=${data}/822_A3/822A3_L2_R1_001_rL8hxUDHZ3jN.fastq.gz
fa822A3_L2_R2=${data}/822_A3/822A3_L2_R2_001_U5JA6jEr3494.fastq.gz

#822 C 2021
fa822C1_L1_R1=${data}/822_C1/822C1_L1_R1_001_ep2dS2QpAM4I.fastq.gz
fa822C1_L1_R2=${data}/822_C1/822C1_L1_R2_001_gubA3bvyzHJ0.fastq.gz
fa822C1_L2_R1=${data}/822_C1/822C1_L2_R1_001_TrH9CpfpMac4.fastq.gz
fa822C1_L2_R2=${data}/822_C1/822C1_L2_R2_001_XkjapjBVa91V.fastq.gz
fa822C2_L1_R1=${data}/822_C2/822C2_L1_R1_001_UGkMHBUbMBIS.fastq.gz
fa822C2_L1_R2=${data}/822_C2/822C2_L1_R2_001_pJjh5PaXssU2.fastq.gz
fa822C2_L2_R1=${data}/822_C2/822C2_L2_R1_001_EoVOaKVUeIRU.fastq.gz
fa822C2_L2_R2=${data}/822_C2/822C2_L2_R2_001_eSnjplSeoQRT.fastq.gz
fa822C3_L1_R1=${data}/822_C3/822C3_L1_R1_001_fzTzPctAU5eI.fastq.gz
fa822C3_L1_R2=${data}/822_C3/822C3_L1_R2_001_DEUVR6BgZWik.fastq.gz
fa822C3_L2_R1=${data}/822_C3/822C3_L2_R1_001_Gs7MOOxAROrQ.fastq.gz
fa822C3_L2_R2=${data}/822_C3/822C3_L2_R2_001_pUPlp5nNS9mi.fastq.gz

#828 C 2021
fa828C1_L1_R1=${data}/828C_1/828C1_L1_R1_001_92jLZEdNYE27.fastq.gz
fa828C1_L1_R2=${data}/828C_1/828C1_L1_R2_001_b6malwUIKmOR.fastq.gz
fa828C1_L2_R1=${data}/828C_1/828C1_L2_R1_001_eVO8OaPO6rD5.fastq.gz
fa828C1_L2_R2=${data}/828C_1/828C1_L2_R2_001_548VC9Dclsso.fastq.gz
fa828C2_L1_R1=${data}/828_C2/828C2_L1_R1_001_2qn9FVAXuJdS.fastq.gz
fa828C2_L1_R2=${data}/828_C2/828C2_L1_R2_001_2HC6cN9XiPwB.fastq.gz
fa828C2_L2_R1=${data}/828_C2/828C2_L2_R1_001_mEHkebaaTPLm.fastq.gz
fa828C2_L2_R2=${data}/828_C2/828C2_L2_R2_001_W9LyK6d9KkdU.fastq.gz
fa828C3_L1_R1=${data}/828_C3/828C3_L1_R1_001_2sObxWvXwbj6.fastq.gz
fa828C3_L1_R2=${data}/828_C3/828C3_L1_R2_001_tiIxD5saDFNJ.fastq.gz
fa828C3_L2_R1=${data}/828_C3/828C3_L2_R1_001_ZtAljFbnP2gp.fastq.gz
fa828C3_L2_R2=${data}/828_C3/828C3_L2_R2_001_MZdYIuLGuk0Z.fastq.gz

#844 C 2021
fa844C1_L1_R1=${data}/844C_1/844C1_L1_R1_001_vaUmNu514Bzh.fastq.gz
fa844C1_L1_R2=${data}/844C_1/844C1_L1_R2_001_krUKwMrER9Dh.fastq.gz
fa844C1_L2_R1=${data}/844C_1/844C1_L2_R1_001_7baBq6QheLQp.fastq.gz
fa844C1_L2_R2=${data}/844C_1/844C1_L2_R2_001_YrAKLUzH70Hv.fastq.gz
fa844C2_L1_R1=${data}/844_C2/844C2_L1_R1_001_cxuvCgbSzle5.fastq.gz
fa844C2_L1_R2=${data}/844_C2/844C2_L1_R2_001_nFrRLngNXk69.fastq.gz
fa844C2_L2_R1=${data}/844_C2/844C2_L2_R1_001_eGVziXCUU5nc.fastq.gz
fa844C2_L2_R2=${data}/844_C2/844C2_L2_R2_001_lsVc74C5eDx0.fastq.gz
fa844C3_L1_R1=${data}/844_C3/844C3_L1_R1_001_QHx5JcuIrij0.fastq.gz
fa844C3_L1_R2=${data}/844_C3/844C3_L1_R2_001_wfrU3LJntHwu.fastq.gz
fa844C3_L2_R1=${data}/844_C3/844C3_L2_R1_001_yPPs3NZg1wFi.fastq.gz
fa844C3_L2_R2=${data}/844_C3/844C3_L2_R2_001_Ct6I2d89VMUc.fastq.gz


#823 D 2021
fa823D1_L1_R1_1=${data}/823D_1/823D1_L1_R1_001_Pz2hOidIcWnC.fastq.gz
fa823D1_L1_R1_2=${data}/823D_1/823D1_L1_R1_001_ZUMePO9FW6GF.fastq.gz
fa823D1_L1_R2_1=${data}/823D_1/823D1_L1_R2_001_yLFumohAT1mu.fastq.gz
fa823D1_L1_R2_2=${data}/823D_1/823D1_L1_R2_001_ysEtwczGazic.fastq.gz
fa823D1_L2_R1_1=${data}/823D_1/823D1_L2_R1_001_2lDwa3LxJwta.fastq.gz
fa823D1_L2_R1_2=${data}/823D_1/823D1_L2_R1_001_59mj4sYN5lm4.fastq.gz
fa823D1_L2_R2_1=${data}/823D_1/823D1_L2_R2_001_WI6pOqzECjuT.fastq.gz
fa823D1_L2_R2_2=${data}/823D_1/823D1_L2_R2_001_lR3LWX9Kt0m5.fastq.gz
fa823D2_L1_R1=${data}/823_D2/843D2_L1_R1_001_gJUmYIZYM9uP.fastq.gz
fa823D2_L1_R2=${data}/823_D2/843D2_L1_R2_001_hMaTTx4H4k3m.fastq.gz
fa823D2_L2_R1=${data}/823_D2/843D2_L2_R1_001_HqG9wf4Nh6Cj.fastq.gz
fa823D2_L2_R2=${data}/823_D2/843D2_L2_R2_001_0hsSsrXdCPo1.fastq.gz
fa823D3_L1_R1=${data}/823_D3/823D3_L1_R1_001_Zu9VpYWbqfxm.fastq.gz
fa823D3_L1_R2=${data}/823_D3/823D3_L1_R2_001_Vi9kOxv8mQlv.fastq.gz
fa823D3_L2_R1=${data}/823_D3/823D3_L2_R1_001_osmdRWLTR1TH.fastq.gz
fa823D3_L2_R2=${data}/823_D3/823D3_L2_R2_001_MQXHzBBy60vy.fastq.gz

########################################################################################################################
########################################################################################################################


#366 b 2020
salmon quant -i EISA_index -l A -r ${fa366b1_L1_R1} ${fa366b1_L2_R1} --validateMappings -p 30 -o salmon_out/366b1 --seqBias --gcBias 
salmon quant -i EISA_index -l A -r ${fa366b2_L1_R1} ${fa366b2_L2_R1} --validateMappings -p 30 -o salmon_out/366b2 --seqBias --gcBias 
salmon quant -i EISA_index -l A -r ${fa366b3_L1_R1} ${fa366b3_L2_R1} --validateMappings -p 30 -o salmon_out/366b3 --seqBias --gcBias 
salmon quant -i EISA_index -l A -r ${fa366b4_L1_R1} ${fa366b4_L2_R1} --validateMappings -p 30 -o salmon_out/366b4 --seqBias --gcBias 

#382 b 2020
salmon quant -i EISA_index -l A -r ${fa382b1_L1_R1} ${fa382b1_L2_R1} --validateMappings -p 30 -o salmon_out/382b1 --seqBias --gcBias 
salmon quant -i EISA_index -l A -r ${fa382b2_L1_R1} ${fa382b2_L2_R1} --validateMappings -p 30 -o salmon_out/382b2 --seqBias --gcBias 
salmon quant -i EISA_index -l A -r ${fa382b3_L1_R1} ${fa382b3_L2_R1} --validateMappings -p 30 -o salmon_out/382b3 --seqBias --gcBias 
salmon quant -i EISA_index -l A -r ${fa382b4_L1_R1} ${fa382b4_L2_R1} --validateMappings -p 30 -o salmon_out/382b4 --seqBias --gcBias 

#775 b 2020
salmon quant -i EISA_index -l A -r ${fa775b1_L1_R1} ${fa775b1_L2_R1} --validateMappings -p 30 -o salmon_out/775b1 --seqBias --gcBias 
salmon quant -i EISA_index -l A -r ${fa775b2_L1_R1} ${fa775b2_L2_R1} --validateMappings -p 30 -o salmon_out/775b2 --seqBias --gcBias 
salmon quant -i EISA_index -l A -r ${fa775b3_L1_R1} ${fa775b3_L2_R1} --validateMappings -p 30 -o salmon_out/775b3 --seqBias --gcBias 
salmon quant -i EISA_index -l A -r ${fa775b4_L1_R1} ${fa775b4_L2_R1} --validateMappings -p 30 -o salmon_out/775b4 --seqBias --gcBias 

#784 b 2020
salmon quant -i EISA_index -l A -r ${fa784b1_L1_R1} ${fa784b1_L2_R1} --validateMappings -p 30 -o salmon_out/784b1 --seqBias --gcBias 
salmon quant -i EISA_index -l A -r ${fa784b2_L1_R1} ${fa784b2_L2_R1} --validateMappings -p 30 -o salmon_out/784b2 --seqBias --gcBias 
salmon quant -i EISA_index -l A -r ${fa784b4_L1_R1} ${fa784b4_L2_R1} --validateMappings -p 30 -o salmon_out/784b4 --seqBias --gcBias 
salmon quant -i EISA_index -l A -r ${fa784b5_L1_R1} ${fa784b5_L2_R1} --validateMappings -p 30 -o salmon_out/784b5 --seqBias --gcBias 

#366 C 2021
salmon quant -i EISA_index -l A -1 ${fa366C1_L1_R1} ${fa366C1_L2_R1} -2 ${fa366C1_L1_R2} ${fa366C1_L2_R2} --validateMappings -p 30 -o salmon_out/366C1 --seqBias --gcBias 
salmon quant -i EISA_index -l A -1 ${fa366C2_L1_R1} ${fa366C2_L2_R1} -2 ${fa366C2_L1_R2} ${fa366C2_L2_R2} --validateMappings -p 30 -o salmon_out/366C2 --seqBias --gcBias 
salmon quant -i EISA_index -l A -1 ${fa366C3_L1_R1} ${fa366C3_L2_R1} -2 ${fa366C3_L1_R2} ${fa366C3_L2_R2} --validateMappings -p 30 -o salmon_out/366C3 --seqBias --gcBias 

#366 D 2021
#salmon quant -i EISA_index -l A -1 ${fa366D1_L1_R1} ${fa366D1_L2_R1} -2 ${fa366D1_L1_R2} ${fa366D1_L2_R2} --validateMappings -p 30 -o salmon_out/366D1 --seqBias --gcBias 
#salmon quant -i EISA_index -l A -1 ${fa366D2_L1_R1} ${fa366D2_L2_R1} -2 ${fa366D2_L1_R2} ${fa366D2_L2_R2} --validateMappings -p 30 -o salmon_out/366D2 --seqBias --gcBias 
#salmon quant -i EISA_index -l A -1 ${fa366D3_L1_R1} ${fa366D3_L2_R1} -2 ${fa366D3_L1_R2} ${fa366D3_L2_R2} --validateMappings -p 30 -o salmon_out/366D3 --seqBias --gcBias 

#822 C 2021
#salmon quant -i EISA_index -l A -1 ${fa822C1_L1_R1} ${fa822C1_L2_R1} -2 ${fa822C1_L1_R2} ${fa822C1_L2_R2} --validateMappings -p 30 -o salmon_out/822C1 --seqBias --gcBias 
#salmon quant -i EISA_index -l A -1 ${fa822C2_L1_R1} ${fa822C2_L2_R1} -2 ${fa822C2_L1_R2} ${fa822C2_L2_R2} --validateMappings -p 30 -o salmon_out/822C2 --seqBias --gcBias 
#salmon quant -i EISA_index -l A -1 ${fa822C3_L1_R1} ${fa822C3_L2_R1} -2 ${fa822C3_L1_R2} ${fa822C3_L2_R2} --validateMappings -p 30 -o salmon_out/822C3 --seqBias --gcBias 

#828 C 2021
salmon quant -i EISA_index -l A -1 ${fa828C1_L1_R1} ${fa828C1_L2_R1} -2 ${fa828C1_L1_R2} ${fa828C1_L2_R2} --validateMappings -p 30 -o salmon_out/828C1 --seqBias --gcBias 
salmon quant -i EISA_index -l A -1 ${fa828C2_L1_R1} ${fa828C2_L2_R1} -2 ${fa828C2_L1_R2} ${fa828C2_L2_R2} --validateMappings -p 30 -o salmon_out/828C2 --seqBias --gcBias 
salmon quant -i EISA_index -l A -1 ${fa828C3_L1_R1} ${fa828C3_L2_R1} -2 ${fa828C3_L1_R2} ${fa828C3_L2_R2} --validateMappings -p 30 -o salmon_out/828C3 --seqBias --gcBias 

#844 C 2021
salmon quant -i EISA_index -l A -1 ${fa844C1_L1_R1} ${fa844C1_L2_R1} -2 ${fa844C1_L1_R2} ${fa844C1_L2_R2} --validateMappings -p 30 -o salmon_out/844C1 --seqBias --gcBias 
salmon quant -i EISA_index -l A -1 ${fa844C2_L1_R1} ${fa844C2_L2_R1} -2 ${fa844C2_L1_R2} ${fa844C2_L2_R2} --validateMappings -p 30 -o salmon_out/844C2 --seqBias --gcBias 
salmon quant -i EISA_index -l A -1 ${fa844C3_L1_R1} ${fa844C3_L2_R1} -2 ${fa844C3_L1_R2} ${fa844C3_L2_R2} --validateMappings -p 30 -o salmon_out/844C3 --seqBias --gcBias 

#822 A 2021
#salmon quant -i EISA_index -l A -1 ${fa822A1_L1_R1} ${fa822A1_L2_R1} -2 ${fa822A1_L1_R2} ${fa822A1_L2_R2} --validateMappings -p 30 -o salmon_out/822A1 --seqBias --gcBias 
#salmon quant -i EISA_index -l A -1 ${fa822A2_L1_R1} ${fa822A2_L2_R1} -2 ${fa822A2_L1_R2} ${fa822A2_L2_R2} --validateMappings -p 30 -o salmon_out/822A2 --seqBias --gcBias 
#salmon quant -i EISA_index -l A -1 ${fa822A3_L1_R1} ${fa822A3_L2_R1} -2 ${fa822A3_L1_R2} ${fa822A3_L2_R2} --validateMappings -p 30 -o salmon_out/822A3 --seqBias --gcBias 

#821 D 2021
#salmon quant -i EISA_index -l A -1 ${fa821D1_L1_R1} ${fa821D1_L2_R1} -2 ${fa821D1_L1_R2} ${fa821D1_L2_R2} --validateMappings -p 30 -o salmon_out/821D1 --seqBias --gcBias 
#salmon quant -i EISA_index -l A -1 ${fa821D2_L1_R1} ${fa821D2_L2_R1} -2 ${fa821D2_L1_R2} ${fa821D2_L2_R2} --validateMappings -p 30 -o salmon_out/821D2 --seqBias --gcBias 
#salmon quant -i EISA_index -l A -1 ${fa821D3_L1_R1} ${fa821D3_L2_R1} -2 ${fa821D3_L1_R2} ${fa821D3_L2_R2} --validateMappings -p 30 -o salmon_out/821D3 --seqBias --gcBias 

#823 D 2021
#salmon quant -i EISA_index -l A -1 ${fa823D1_L1_R1_1} ${fa823D1_L1_R1_2} ${fa823D1_L2_R1_1} ${fa823D1_L2_R1_2} -2 ${fa823D1_L1_R2_1} ${fa823D1_L1_R2_2} ${fa823D1_L2_R2_1} ${fa823D1_L2_R2_2} --validateMappings -p 30 -o salmon_out/823D1 --seqBias --gcBias 
#salmon quant -i EISA_index -l A -1 ${fa823D2_L1_R1} ${fa823D2_L2_R1} -2 ${fa823D2_L1_R2} ${fa823D2_L2_R2} --validateMappings -p 30 -o salmon_out/823D2 --seqBias --gcBias 
#salmon quant -i EISA_index -l A -1 ${fa823D3_L1_R1} ${fa823D3_L2_R1} -2 ${fa823D3_L1_R2} ${fa823D3_L2_R2} --validateMappings -p 30 -o salmon_out/823D3 --seqBias --gcBias 


conda deactivate




