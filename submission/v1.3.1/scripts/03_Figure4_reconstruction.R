suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
args <- commandArgs(trailingOnly=TRUE)
root <- if (length(args)) args[1] else getwd()
if (!dir.exists(file.path(root,'derived_tables'))) stop('Run from package root or pass package root as the first argument.')
out <- file.path(root,'figures'); dir.create(out,showWarnings=FALSE,recursive=TRUE)
formal <- fread(file.path(root,'derived_tables','HTN_MVP_SBP_formal_coloc_ALL_PRIORS.tsv'))
corr <- fread(file.path(root,'derived_tables','Figure4B_MVP_regional_correlation.tsv'))
ukb <- fread(file.path(root,'derived_tables','Figure4C_discovery_vs_UKB_M50_H4.tsv'))
formal[,Region:=ifelse(region=='chr10_region6','Chromosome 10','Chromosome 18')]
formal[,Prior:=factor(p12,levels=c(1e-6,1e-5,5e-5),labels=c('1e-6','1e-5','5e-5'))]
fl <- melt(formal,id.vars=c('Region','Prior','p12'),measure.vars=c('PP.H3','PP.H4'),variable.name='Hypothesis',value.name='Posterior')
fl[,Hypothesis:=ifelse(Hypothesis=='PP.H3','H3: distinct variants','H4: shared variant')]
ref <- fl[p12==1e-5]; ref[,label_y:=ifelse(Posterior>0.8,Posterior-0.08,Posterior+0.08)]
pA <- ggplot(fl,aes(Prior,Posterior,group=Hypothesis,color=Hypothesis,shape=Hypothesis))+geom_hline(yintercept=.5,linetype='dashed')+geom_vline(xintercept=2,linetype='dotted')+geom_line()+geom_point(size=2.7)+geom_text(data=ref,aes(y=label_y,label=sprintf('%.3f',Posterior)),show.legend=FALSE,size=3)+facet_wrap(~Region,nrow=1)+scale_y_continuous(limits=c(0,1))+labs(title='HTN-MVP-SBP regional colocalization',subtitle='H3 and H4 posterior probabilities across shared-variant priors',x='Shared-variant prior (p12)',y='Posterior probability',color=NULL,shape=NULL)+theme_classic(base_size=10)+theme(legend.position='top')
pB <- ggplot(corr,aes(Metric,Correlation,fill=Region))+geom_col(position=position_dodge(.75),width=.62)+geom_text(aes(label=sprintf('%.3f',Correlation)),position=position_dodge(.75),vjust=-.45,size=3)+scale_y_continuous(limits=c(0,.72))+labs(title='Regional concordance with MVP-SBP',subtitle='Aligned regional Z-score correlations',x=NULL,y='Correlation',fill=NULL)+theme_classic(base_size=10)+theme(legend.position='top')
ord <- c('NT5C2','RMC1','CNNM2','MARCKSL1P1'); ukb[,Gene:=factor(Gene,levels=rev(ord))]
ul <- rbind(ukb[,.(Gene,Dataset='Discovery CDD',H4=H4_discovery)],ukb[,.(Gene,Dataset='UKB M50',H4=H4_UKB)])
pC <- ggplot()+geom_vline(xintercept=.5,linetype='dashed')+geom_segment(data=ukb,aes(x=H4_UKB,xend=H4_discovery,y=Gene,yend=Gene),color='grey40')+geom_point(data=ul,aes(H4,Gene,shape=Dataset),size=2.8)+scale_x_continuous(limits=c(0,.82),breaks=c(0,.25,.5,.75))+labs(title='CDD-side colocalization in discovery and UKB M50',subtitle='Reference prior p12 = 1e-5',x='PP.H4',y=NULL,shape=NULL)+theme_classic(base_size=10)+theme(legend.position='top')
fig <- (pA/(pB|pC))+plot_layout(heights=c(1.05,1))+plot_annotation(tag_levels='A')
ggsave(file.path(out,'Figure4_external_genetic_evaluation_RECONSTRUCTED.pdf'),fig,width=10.5,height=8)
