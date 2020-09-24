import numpy as np
import scanpy as sc
import anndata as ad
import pandas as pd
import xlsxwriter
import pickle
from rpy2 import robjects

def mast_de_groups(adata, groupby, save):
    '''Compute differential expression with the MAST package by treatment covariate within clusters provided as "groupby" and export as excel file'''

    robjects.r('''
        mast_de_r <- function(adata, obs, var, clusters, groupby){
            #Prepare data sets for SingleCellExperiment data structure conversion
            #obs['wellKey'] = row.names(obs)
            #var['primerid'] = row.names(var)
            print('Deploying to R...')
            #Convert to SingleCellExperiment type
            #sca <- FromMatrix(exprsArray=data_mat, fData=var)
            sca <- SceToSingleCellAssay(adata, class = "SingleCellAssay")
            #Compute Gene detection rate
            colData(sca)$n_genes = scale(colData(sca)$n_genes)

            #Create a vector that will hold all the DE results
            output <- vector("list", length(clusters))

            count <- 0
            print('Begin computation...')
            #Loop over all louvain clusters
            for (i in clusters){
                count <- count+1
                print(i)
                #Create data subsets which should be used for testing
                if (groupby=='louvain_final') {
                    sca_sub <- subset(sca, with(colData(sca), louvain_final==i))
                } else if (groupby=='louvain_r1') {
                    sca_sub <- subset(sca, with(colData(sca), louvain_r1==i))
                } else {
                    stop()
                }

                #Filter out non-expressed genes in the subset
                sca_sub <- sca_sub[rowSums(assay(sca_sub)) != 0, ]

                #Define & run hurdle model
                zlmCond <- zlm(formula = ~condition + n_genes, sca=sca_sub)
                summaryCond <- summary(zlmCond, doLRT='conditionStress')

                summaryDt <- summaryCond$datatable

                result <- merge(summaryDt[contrast=='conditionStress' & component=='H',.(primerid, `Pr(>Chisq)`)], #p-vals
                                 summaryDt[contrast=='conditionStress' & component=='logFC', .(primerid, coef)], #logFC coefficients
                                 by='primerid') 

                #Correct for multiple testing (FDR correction) and filtering
                result[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]
                result[,coef:=result[,coef]/log(2)]
                names(result) <- c("gene", "pval", "log2FC", "qval")
                result = result[order(result$qval),]

                output[[count]] <- result

            }
            return(output)
        }
    ''')
    
    mast_de = robjects.globalenv['mast_de_r']
    

    
    #Create new Anndata object for use in MAST with non-batch corrected data as before
    adata_test = adata.copy()
    adata_test.X = adata.raw.X
    adata_test.obs['n_genes'] = (adata_test.X > 0).sum(1) 
    
    obs = adata_test.obs
    var = adata_test.var
    clusters = list(adata_test.obs[groupby].cat.categories)
    
    expr_dict = {
        adata_test.var.index[i]:{} for i in range(adata_test.shape[1])
    }
    expr_dict_stress = {
        adata_test.var.index[i]:{} for i in range(adata_test.shape[1])
    }
    expr_dict_ctrl = {
        adata_test.var.index[i]:{} for i in range(adata_test.shape[1])
    }
    
    for clust in adata_test.obs[groupby].cat.categories:
        expr = np.mean(adata_test[adata_test.obs[groupby] == clust].X, axis=0)
        
        expr_stress = np.mean(adata_test[
            (adata_test.obs['condition']=='Stress')
            & (adata_test.obs[groupby] == clust)
        ].X, axis=0)
        
        expr_ctrl = np.mean(adata_test[
            (adata_test.obs['condition']=='Control')
            & (adata_test.obs[groupby] == clust)
        ].X, axis=0)
        
        for i, gene in enumerate(adata_test.var.index):
            expr_dict[gene][clust] = expr[i]
            expr_dict_stress[gene][clust] = expr_stress[i]
            expr_dict_ctrl[gene][clust] = expr_ctrl[i]
    
    result = mast_de(adata_test, obs, var, clusters, groupby)
    result = {clusters[i]:datframe for i, datframe in enumerate(result)}
    
    writer = pd.ExcelWriter(save, engine='xlsxwriter')
    print('Number of significant DE genes:')    
    for clust in clusters:
        
        result[clust]['meanExpr'] = [
            expr_dict[gene][clust] for gene in result[clust]['gene'].values
        ]
        
        result[clust]['meanExprStress'] = [
            expr_dict_stress[gene][clust] for gene in result[clust]['gene'].values
        ]
        
        result[clust]['meanExprCtrl'] = [
            expr_dict_ctrl[gene][clust] for gene in result[clust]['gene'].values
        ]
        
        result[clust].to_excel(writer,sheet_name=str(clust))
        print(clust+':', np.sum([result[clust]['qval']<0.05]))

    writer.save()

    return result

def mast_de_bulk(adata, save):
    '''Compute differential expression with the MAST package by treatment covariate on the whole dataset'''
    
    robjects.r('''
        mast_de_r <- function(adata, obs, var){
            #Prepare data sets for SingleCellExperiment data structure conversion
            #obs['wellKey'] = row.names(obs)
            #var['primerid'] = row.names(var)
            print('Deploying to R...')
            #Convert to SingleCellExperiment type
            #sca <- FromMatrix(exprsArray=data_mat, fData=var)
            sca <- SceToSingleCellAssay(adata, class = "SingleCellAssay")
            #Compute Gene detection rate
            colData(sca)$n_genes = scale(colData(sca)$n_genes)

            #Create a vector that will hold all the DE results

            count <- 0
            print('Begin computation...')
            #Define & run hurdle model
            zlmCond <- zlm(formula = ~condition + n_genes, sca=sca)
            summaryCond <- summary(zlmCond, doLRT='conditionStress')
            summaryDt <- summaryCond$datatable
            
            result <- merge(summaryDt[contrast=='conditionStress' & component=='H',.(primerid, `Pr(>Chisq)`)], #p-vals
                             summaryDt[contrast=='conditionStress' & component=='logFC', .(primerid, coef)], #logFC coefficients
                             by='primerid') 
                             
            #Correct for multiple testing (FDR correction) and filtering
            result[,FDR:=p.adjust(`Pr(>Chisq)`, 'fdr')]
            result[,coef:=result[,coef]/log(2)]
            names(result) <- c("gene", "pval", "log2FC", "qval")
            result = result[order(result$qval),]
            output <- result
                
            return(output)
        }
    ''')
    
    mast_de = robjects.globalenv['mast_de_r']
    

    
    #Create new Anndata object for use in MAST with non-batch corrected data as before
    adata_test = adata.copy()
    adata_test.X = adata.raw.X
    adata_test.obs['n_genes'] = (adata_test.X > 0).sum(1) 
    
    obs = adata_test.obs
    var = adata_test.var
    
    expr_dict = {
        adata_test.var.index[i]:{} for i in range(adata_test.shape[1])
    }
    
    expr_dict_stress = {
        adata_test.var.index[i]:{} for i in range(adata_test.shape[1])
    }
    
    expr_dict_ctrl = {
        adata_test.var.index[i]:{} for i in range(adata_test.shape[1])
    }
    
    expr = np.mean(adata_test.X, axis=0)
    
    expr_stress = np.mean(
        adata_test[(adata_test.obs['condition']=='Stress')].X, axis=0
    )
    
    expr_ctrl = np.mean(
        adata_test[(adata_test.obs['condition']=='Control')].X, axis=0
    )
    
    for i, gene in enumerate(adata_test.var.index):
        expr_dict[gene] = expr[i]
        expr_dict_stress[gene] = expr_stress[i]
        expr_dict_ctrl[gene] = expr_ctrl[i]
    
    result = mast_de(adata_test, obs, var)
    result
    writer = pd.ExcelWriter(save, engine='xlsxwriter')  
    
    result['meanExpr'] = [
        expr_dict[gene] for gene in result['gene'].values
    ]
    
    result['meanExprStress'] = [
        expr_dict_stress[gene] for gene in result['gene'].values
    ]
    
    result['meanExprCtrl'] = [
        expr_dict_ctrl[gene] for gene in result['gene'].values
    ]
    result.to_excel(writer)

    writer.save()

    return result