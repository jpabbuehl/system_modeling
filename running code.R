eset<-initialize_eset()
eset<-rm_outlier_replicate(eset)
subpopulations<-knn_subpopulation(eset)
relevant_probeset<-betr_filter(eset,druglist)
signatures_E2<-model_fit(eset,relevant_probeset,choice='E2')
signatures_Marcus<-model_fit(eset,relevant_probeset,choice='Marcus')
signatures_M<-model_fit(eset,relevant_probeset,choice='M')
save(signatures_E2,file='signatures E2.rda')
save(signatures_Marcus,file='signatures Marcus.rda')
save(signatures_M,file='signatures M.rda')

Optimizatin signatures level
load(file='signatures E2.rda')
signatures_E2<-signature_optimization('Methotrexate',signatures_E2)
etc....

up<-vector()
down<-vector()
drug<-vector()
for(i in 1:length(signatures_E2)){
  stat<-table(signatures_E2[[i]][2])
  up<-append(up,stat[2])
  down<-append(down,stat[1])
  drug<-append(drug,names(signatures_E2)[i])
}
info<-data.frame(drug,up,down)

save(signatures_E2,file='RF optimized signatures.rda')
save(info,file='summary RF optimized signatures.rda')
save(subpopulations,file='knn subpopulation.rda')
save(targeting_matching,file='drug drug interaction matrix.rda')
# Here now
source('main.R')
source('data.R')
load('RF optimized signatures.rda')
load('knn subpopulation.rda')
load('drug drug interaction matrix.rda')
eset<-initialize_eset()
eset<-rm_outlier_replicate(eset)
targeting_matching<-population_targeting(eset,signatures_E2)
#current top
prediction(eset,subpopulations,targeting_matching,subpopulation=TRUE,spectrum=FALSE,scale=TRUE,top=20)
