#Overrepresentation analyses (ORA)

library("clusterProfiler")
library("org.Hs.eg.db")
library("ReactomePA")
library("DOSE")

#Generating ENTREZID gene lists of upregulated and downregulated genes from DGE analysis
upreg <- rownames(filtered_res05_up)
genes_up <- bitr(upreg, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
entrez_up <- genes_up$ENTREZID
downreg <- rownames(filtered_res05_down)
genes_down <- bitr(downreg, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb="org.Hs.eg.db")
entrez_down <- genes_down$ENTREZID

#KEGG ORA upregulated genes
kegg_up <- enrichKEGG(entrez_up)
results_kegg_up <- as.data.frame(kegg_up)

#KEGG ORA downregulated genes
kegg_down <- enrichKEGG(entrez_down)
results_kegg_down <- as.data.frame(kegg_down)

#Reactome ORA upregulated genes
e_up <- enrichPathway(gene=entrez_up, readable=TRUE)
results_react_up <- as.data.frame(e_up)

#Reactome ORA downregulated genes
e_down <- enrichPathway(gene=entrez_down, readable=TRUE)
results_react_down <- as.data.frame(e_down)

#GO-MF ORA upregulated genes
ego_up_mf <- enrichGO(gene = entrez_up, OrgDb = org.Hs.eg.db, ont = "MF", readable = TRUE)
ego_up_mf_simplify <- simplify(ego_up_mf)
results_go_up_mf <- as.data.frame(ego_up_mf_simplify)

#GO-MF ORA downregulated genes
ego_down_mf <- enrichGO(gene = entrez_down, OrgDb = org.Hs.eg.db, ont = "MF", readable = TRUE)
ego_down_mf_simplify <- simplify(ego_down_mf)
results_go_down_mf <- as.data.frame(ego_down_mf_simplify)

#GO-BP ORA upregulated genes
ego_up_bp <- enrichGO(gene = entrez_up, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
ego_up_bp_simplify <- simplify(ego_up_bp)
results_go_up_bp <- as.data.frame(ego_up_bp_simplify)

#GO-BP ORA downregulated genes
ego_down_bp <- enrichGO(gene = entrez_down, OrgDb = org.Hs.eg.db, ont = "BP", readable = TRUE)
ego_down_bp_simplify <- simplify(ego_down_bp)
results_go_down_bp <- as.data.frame(ego_down_bp_simplify)

#GO-CC ORA upregulated genes
ego_up_cc <- enrichGO(gene = entrez_up, OrgDb = org.Hs.eg.db, ont = "CC", readable = TRUE)
ego_up_cc_simplify <- simplify(ego_up_cc)
results_go_up_cc <- as.data.frame(ego_up_cc_simplify)
    
#GO-CC ORA downregulated genes
ego_down_cc <- enrichGO(gene = entrez_down, OrgDb = org.Hs.eg.db, ont = "CC", readable = TRUE)
ego_down_cc_simplify <- simplify(ego_down_cc)
results_go_down_cc <- as.data.frame(ego_down_cc_simplify)

#DO ORA upregulated genes
do_up <- enrichDO(entrez_up)
results_do_up <- as.data.frame(do_up)

#DO ORA downregulated genes
do_down <- enrichDO(entrez_down)
results_do_down <- as.data.frame(do_down)