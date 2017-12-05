clear
clc
load inputdata.mat

p=0.25;
Option=2;

disease_sim02  = cosSim(md_adjmat02');
miRNA_sim02  = cosSim(md_adjmat02);
        

scores=reductionModel(disease_sim02, miRNA_sim02, md_adjmat02',Option,p);

[KMDR_rank,KMDR_rank_known] =Rank_miRNAs( scores',md_adjmat02,miRNA_list02,disease_list02 );

Write_file( KMDR_rank )
%clear w v p Option scores