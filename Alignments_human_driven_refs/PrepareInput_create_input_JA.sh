#!/bin/bash

module purge
module load gcc/4.9.3-gold
module load PYTHON/3.6.3

#We keep the species names for each one of the primates used for annotation from their path
INDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Alignments_refs/Input_clusters/
OUTDIR=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/human_driven_results/Alignments_refs/Input_fastas/
SRC=/scratch/devel/avalenzu/PhD_EvoGenom/GenomPhenom200primates/src/Alignments_human_driven_refs/
mkdir -p ${OUTDIR}

#create variable file
touch PrepareInput_create_input_JA.txt

#Loop for human genes
Genes="ADGRA3 APEX1 APEX2 APH1A APH1B API5 APIP APLN APLNR APLP2 APMAP APOA2 APOA4 APOB APOBEC1 APOBEC2 APOBEC3F ARL6IP5 ARRDC5 ARSA ARVCF ASB14 ASTL ATG9B ATIC ATP6AP1L ATP6AP2 BACH2 BAD BAG1 BAG2 BAG3 BPGM BPHL C1QTNF6 C1QTNF7 C1QTNF8 CAMKMT CAMKV CASP1 CASP2 CCDC102B CCDC103 CCDC177 CCDC178 CCDC8 CCDC88C CCDC89 CCDC90B CCDC91 CCND1 CCND2 CDC20 CDC23 CDC25C CDK5R1 CDK5R2 CDK5RAP2 CEACAM5 CHIC2 CLDND2 CLEC10A CNBP CNDP1 COX7A2 CPNE9 CS CXCL11 CXCL12 CXCL13 CXCL14 CXCL17 CYP2E1 CYP2J2 CYP2S1 CYP2U1 DIABLO DIAPH1 DIAPH2 DIAPH3 DICER1 DIDO1 DIMT1 DNAH5 DUSP2 DUSP28 EFR3A EFR3B ENTHD1 ENTPD1 ERGIC1 ERGIC2 ERGIC3 ETV2 ETV3 ETV3L EXOSC6 EXOSC7 EXOSC8 EXOSC9 EXT1 FAM217B FBXL2 FBXL22 FBXL3 FBXL4 FGF1 FUT7 FUT8 FUT9 GALNT9 GALNTL5 GALNTL6 GALP GET4 GFAP GJB7 GJC1 GJC2 GJC3 GOLGA3 GOLGA4 GOLGA5 GPR149 GPR150 GRIP1 GRIP2 GRIPAP1 HM13 HRH4 HRNR HS1BP3 HS2ST1 HTR1A HTR1B HTR1D IDNK IDO1 IGFBP7 IGSF9 IKBIP IMP3 IMP4 IMPA1 INTS9 INTU INVS IP6K1 JADE2 JADE3 KIRREL2 KLHL3 KLHL36 KLHL38 KRT3 KRT37 KRT38 KRT39 LGALS7 LGALS8 LLGL2 LRRC3 LRRC38 LRRC3B LRRC3C LYRM1 LYRM2 LYRM4 LYRM7 LYSMD1 MAD2L2 MARS2 MARVELD1 MFAP2 MFAP3L MGMT MIS18A MIS18BP1 MISP3 MORN3 MORN4 MTIF2 MTIF3 MTMR10 MTMR11 MYBPH MYBPHL MYCBP2 MYCBPAP MYOZ2 NDUFA1 NDUFA2 NDUFA3 NELL1 NELL2 NEMF NEMP1 NIPA2 NIPAL1 NOS1AP NRG1 NUDT11 OFD1 OGDHL OR4F6 OR4K15 OR6B2 PABPC4 PABPC4L PABPC5 PABPN1L PACRG PACRGL PCGF3 PCGF5 PCGF6 PCID2 PCIF1 PCYOX1 PCYT1A PHF11 PHF12 PIGV PKD2 PKD2L2 PLEC PLG PODNL1 PRAMEF12 PRKAA2 PRKAB2 PRR3 PRR35 PRR36 PRR5-ARHGAP8 PRR5L PTH PTK2 PTK2B PTK6 PTK7 RAB21 RAB22A RAB23 RAB24 RAB25 RAB26 RAB27A RARRES2 RBBP7 RORB RORC ROS1 RPL7 RPL7A RPL7L1 RPL8 RRP7A RRP8 RRP9 RRS1 SAA1 SAA2-SAA4 SCD SCFD1 SGMS2 SGO1 SGO2 SGPL1 SGPP1 SKI SKIV2L SLC34A2 SLC34A3 SMIM19 SMIM23 SMIM24 SNRNP70 SNRPA SNRPA1 SNRPB SNRPB2 SSUH2 STK3 STK38 STK39 TACC2 TACC3 TAF1A TAS1R2 TBXA2R TES TET1 TET3 TIMM23 TIMM29 TIMM44 TIMM50 TMEM150A TMEM68 TMEM69 TOM1 TOMM20 TOMM20L TOMM22 TOR2A TOR3A ZNF410 ZNF418 ZNF420 ZNF423 ZNF425 ZNF426 ZNF428 ZNF429 ZNF514 ZNF516 ZNF517 ZNF518A ZNF736"
echo $Genes | tr " " "\n" | while read gene;
do for filepath in $(ls ${INDIR}${gene}/*tar.gz);
do echo "$filepath">> PrepareInput_create_input_JA.txt
echo $filepath
done
done
