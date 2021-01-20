
dir_inputGeneTreesMetadata='/data/home/mass/tmass/stylophora_yotam/data/stuff/gene-trees-biomineralization/gene-trees-tests-11/Results_Dec09-biomin-OGs_graphs/graphs1'
inputPrefix="tree.biomin"

#outidx = '/data/home/mass/tmass/stylophora_yotam/data/stuff/gene-trees-biomineralization/gene-trees-tests-11/Results_Dec09-biomin-OGs_graphs/graphs1/species-level-1000per'
outidx = '/data/home/mass/tmass/stylophora_yotam/data/stuff/gene-trees-biomineralization/gene-trees-tests-11/Results_Dec09-biomin-OGs_graphs/graphs1/group-level-1000per'

# use one of the following options: "reoccurrence_score_" or "reoccurrenceBySpeciesGroup_score_"
#scoreFunction="reoccurrence_score_"
scoreFunction="reoccurrenceBySpeciesGroup_score_"

col_id = 'Newick_label'
col_select = 'select_'
col_species= "species"

#col_cluster_ids = 'ids.node.children'
col_cluster_ids = 'ids'
col_cluster_reference = 'clades'
cluster_reference= 'complex corals; robust corals'
col_cluster_species = 'reference.species'
cluster_species = 'CniMergedstylophoraP'

species_for_test=c("CniAcroporaM", "CniAcroporaD" ,"CniMergedstylophoraP") #,"cniOculinaP") # the species that are included in the selected genes set
species_group = c("Complex","Complex","Robust")

###########################

library(ggplot2)

prefix1 = unique(sapply(list.files(path=dir_inputGeneTreesMetadata),function(x) regmatches(x,regexpr(paste0(inputPrefix,".OG\\d+"),x,perl=T))[1]))
prefix1 = prefix1[!is.na(prefix1)]

geneCluster2 = list()
clusters2 = list()
metadata2 = list()
selected_ids=c()
all_ids=c()
for(i in 1:length(prefix1)){
    cat(prefix1[i],"\n")
    f_metadata=paste0(dir_inputGeneTreesMetadata,"/",prefix1[i],".a.ggtree-metadata.txt")
    f_clusters=paste0(dir_inputGeneTreesMetadata,"/",prefix1[i],".a.orthologs-clusters.txt")
    if(file.exists(f_metadata) & file.exists(f_clusters)){
        metadata1 = read.csv(f_metadata,sep="\t",stringsAsFactors=F,header=T)
        metadata2[[prefix1[i]]] = metadata1[,c(col_id,col_select)]
        all_ids = c(all_ids,metadata1[,col_id])
        w = metadata1[metadata1[,col_species] %in% species_for_test,col_select]
        selected_ids = c(selected_ids,w[!is.na(w)])
        finfo = file.info(f_clusters)
        if(finfo$size > 3){ # file size > 3 bytes
            clusters1 = read.csv(f_clusters,sep="\t",stringsAsFactors=F,header=T)
            clusters2[[prefix1[i]]] = clusters1[  (clusters1[,col_cluster_reference] == cluster_reference)&(clusters1[,col_cluster_species] == cluster_species),]
            geneCluster1 = data.frame(id=c(),group=c())
            if(nrow(clusters2[[prefix1[i]]]) > 0){
                for(j in 1:nrow(clusters2[[prefix1[i]]])){
                    ids0 = strsplit(clusters2[[prefix1[i]]][j,col_cluster_ids],split="\\s*;\\s*",perl=T)[[1]]
                    ids1 = data.frame(id = ids0, group = rep(j,times=length(ids0)),stringsAsFactors=F)
                    geneCluster1 = rbind(geneCluster1,ids1)
                }
                geneCluster2[[prefix1[i]]] = geneCluster1
            }
        }
    }else{cat("file does not exist !!\n")}
}
length(geneCluster2)

reoccurrence_score_ = function(clusters,ids,minMatch = 2){
    groupSize = aggregate(clusters$id,by=list(clusters$group),length)
    names(groupSize) = c('group','size')
    c = clusters$id %in% ids
    score1=0
    matchesPerGroup2=NA
    if(sum(c) > 0){
        matchesPerGroup = aggregate(clusters[c,'id'],by=list(clusters[c,'group']),length)
        names(matchesPerGroup) =  c('group','match')
        matchesPerGroup1 = merge(matchesPerGroup,groupSize,by='group')
        matchesPerGroup1$density = apply(matchesPerGroup1[,c('match','size')],1,function(x) x[1]/x[2])
        z = matchesPerGroup1[matchesPerGroup1$match >= minMatch,]
        if(nrow(z) > 0){
            matchesPerGroup2 = z[z$density==max(z$density),][1,]
            score1 = matchesPerGroup2[1,'density']
        }else{
            matchesPerGroup2 = matchesPerGroup1[matchesPerGroup1$density==max(matchesPerGroup1$density),][1,]
        }
    }
    return(list(score1,matchesPerGroup2,unique(clusters$id[c])))
}

reoccurrenceBySpeciesGroup_score_ = function(clusters,ids,minMatch = 2,species_for_test1=species_for_test,species_group1=species_group){
    groupSize = aggregate(clusters$id,by=list(clusters$group),length)
    names(groupSize) = c('group','size')
    clusters$species=sapply(clusters$id,function(x) gsub("_.*","",x,perl=T))
    indices=match(clusters$species,species_for_test1)
    clusters$speciesGroup=species_group1[indices] # assigning species group name for each species, or NA if the species doesn't apper in species_for_test variables.
    c = clusters$id %in% ids
    score1=0
    matchesPerGroup2=NA
    if(sum(c) > 0){
        #matchesPerGroup = aggregate(clusters[c,'id'],by=list(clusters[c,'group']),length)
        clusters_ = clusters[c,]
        matchesPerGroup = aggregate(1:nrow(clusters_),by=list(clusters_[,'group']),
            function(x){u=unique(clusters_[x,"speciesGroup"]); u=u[!is.na(u)]; return(ifelse(length(u)==length(unique(species_group1)),length(u),0))})
        names(matchesPerGroup) =  c('group','match')
        matchesPerGroup1 = merge(matchesPerGroup,groupSize,by='group')
        matchesPerGroup1$density = apply(matchesPerGroup1[,c('match','size')],1,function(x) x[1]/x[2])
        z = matchesPerGroup1[matchesPerGroup1$match >= minMatch,]
        if(nrow(z) > 0){
            matchesPerGroup2 = z[z$density==max(z$density),][1,]
            score1 = matchesPerGroup2[1,'density']
        }else{
            matchesPerGroup2 = matchesPerGroup1[matchesPerGroup1$density==max(matchesPerGroup1$density),][1,]
        }
    }
    return(list(score1,matchesPerGroup2,unique(clusters$id[c])))
}

if(scoreFunction == "reoccurrenceBySpeciesGroup_score_"){
    reoccurrence_score = reoccurrenceBySpeciesGroup_score_
}else{
    reoccurrence_score = reoccurrence_score_
}

reoccurrence_score1 = list()
reoccurrence_matches=list()
reoccurrence_obs = data.frame(group=c(), match=c(), size=c(), density=c(), OG=c())
for(i in 1:length(geneCluster2)){ # visit each gene tree
    list1 = reoccurrence_score(geneCluster2[[i]],selected_ids,minMatch = 2)
    reoccurrence_score1[[i]] = list1[[1]]
    reoccurrence_matches[[i]] =  list1[[3]]
    if(is.data.frame(list1[[2]])){
        list1[[2]]$OG = names(geneCluster2)[i]
        reoccurrence_obs = rbind(reoccurrence_obs,list1[[2]])
    }
}
reoccurrence_scores_observed = unlist(reoccurrence_score1)
write.table(reoccurrence_obs,file=paste0(outidx,"-reoccurrence_observed.txt"),sep="\t",row.names=F)

reoccurrence_score_expected = list()
reoccurrence_expect = data.frame(group=c(), match=c(), size=c(), density=c(), OG=c(), permutation_num=c())
expected1 = c()
perm=1000
genes_for_test = all_ids[grepl(paste0(species_for_test,collapse="|"),all_ids,perl=T)]
for(j in 1:perm){
    cat(j,"\n")
    rand_ids = sample(genes_for_test,size=length(selected_ids),replace=F) # sample(all_ids,size=length(selected_ids),replace=F)
    reoccurrence_score_expected0 = list()
    for(i in 1:length(geneCluster2)){ # visit each gene tree
        list1 = reoccurrence_score(geneCluster2[[i]],rand_ids,minMatch = 2)
        reoccurrence_score_expected0[[i]] = list1[[1]]
        if(is.data.frame(list1[[2]])){
            list1[[2]]$OG = names(geneCluster2)[i]
            list1[[2]]$permutation_num = paste0(j," (",perm,")")
            reoccurrence_expect = rbind(reoccurrence_expect,list1[[2]])
        }
    }
    reoccurrence_score_expected[[j]] = unlist(reoccurrence_score_expected0) 
    expected1  = c(expected1,sum(reoccurrence_score_expected[[j]]))
}
obs1 = sum(reoccurrence_scores_observed)
p = sum(expected1 >= obs1)/perm
write.table(reoccurrence_expect,file=paste0(outidx,"-reoccurrence_expected.txt"),sep="\t")

totalClusters1=c()
for(z in 1:length(geneCluster2)){totalClusters1=c(totalClusters1,length(unique(geneCluster2[[z]]$group)))}
sum(totalClusters1)

save.image(paste0(outidx,".RData"))

# graphs:

reoccurrence_obs2 = reoccurrence_obs[reoccurrence_obs$match >= 0,]
reoccurrence_obs2$permutation_num='0'
reoccurrence_obs2$experiment='observed'
reoccurrence_expect2 = reoccurrence_expect[reoccurrence_expect$match >= 0,]
reoccurrence_expect2$experiment='permutation'
reoccurrence2=rbind(reoccurrence_obs2,reoccurrence_expect2)
reoccurrence2$density_ = sapply(reoccurrence2$density,function(x) ifelse(x>0.3,0.3,x))
q = quantile(reoccurrence_obs2[reoccurrence_obs2$match > 1,'density'],0.1)
reoccurrence2$enrichedClade = apply(reoccurrence2[,c('match','density')],1,function(x) ifelse((x[1]>1) & (x[2]>q),1,0))

reoccurrence2[reoccurrence2$permutation_num==0 & reoccurrence2$match>0,]

gg1 = ggplot(reoccurrence2,aes(x=as.factor(experiment),y=match,color=density_,alpha=0.05))
gg1 = gg1 + geom_jitter(height=0.05,width=0.25) + scale_y_continuous(breaks =0:9) + scale_color_gradient2(low='green',mid='black',high='red',midpoint=0.1)
pdf(paste0(outidx,"_barplot.pdf"))
gg1 + theme_bw()
dev.off()

gg2 = ggplot(reoccurrence2,aes(x=enrichedClade,group=as.factor(experiment),color=experiment))
gg2 = gg2 + geom_histogram(aes(y=..density..),alpha=0.5, position="identity",bins=2,fill=NA)
#gg2 = gg2 + xlab(paste0("enriched clade (Complex and Robust known genes in clades, high known genes density)",' p=',p))
gg2 = gg2 + xlab(paste0("enriched clade (high known genes density)",' p=',p))
pdf(paste0(outidx,"_hist.pdf"))
gg2 + theme_bw()
dev.off()

reoccurrence_obs3 = aggregate(reoccurrence_obs2[,c('match')],by=list(reoccurrence_obs2$match),length)
names(reoccurrence_obs3) = c('sekeltal_genes_count_in_cluster','clusters')
reoccurrence_obs3$experiment = 'observed'
reoccurrence_obs3$clusters_count = reoccurrence_obs3$clusters

reoccurrence_expect3 = aggregate(reoccurrence_expect2[,c('match')],by=list(reoccurrence_expect2$match),length)
names(reoccurrence_expect3) = c('sekeltal_genes_count_in_cluster','clusters')
reoccurrence_expect3$experiment = 'permutation-test'
reoccurrence_expect3$clusters_count = reoccurrence_expect3$clusters/100

reoccurrence3 = rbind(reoccurrence_obs3,reoccurrence_expect3)
reoccurrence3 = reoccurrence3[reoccurrence3$sekeltal_genes_count_in_cluster  > 1,]

save.image(paste0(outidx,".RData"))