if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("pathview")
BiocManager::install("GEOquery")

library(pathview)
library(GEOquery)

gse <- getGEO("GSE16873")
data_gse16873 <- exprs(gse[[1]])

# Define the KEGG pathway ID for the urea cycle
urea_cycle_id <- "map00220"

# Load gene expression data
data(gse16873)

data(gse16873.d)
data(demo.paths)

#KEGG view: gene data only
i <- 1
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id =
                     demo.paths$sel.paths[i], species = "hsa", out.suffix = "gse16873",
                   kegg.native = TRUE)
str(pv.out)
head(pv.out$plot.data.gene)
#result PNG file in current directory

#Graphviz view: gene data only
pv.out <- pathview(gene.data = gse16873.d[, 1], pathway.id =
                     demo.paths$sel.paths[i], species = "hsa", out.suffix = "gse16873",
                   kegg.native = FALSE, sign.pos = demo.paths$spos[i])
#result PDF file in current directory

#KEGG view: both gene and compound data
sim.cpd.data=sim.mol.data(mol.type="cpd", nmol=3000)
i <- 3
print(demo.paths$sel.paths[i])
pv.out <- pathview(gene.data = gse16873.d[, 1], cpd.data = sim.cpd.data,
                   pathway.id = demo.paths$sel.paths[i], species = "hsa", out.suffix =
                     "gse16873.cpd", keys.align = "y", kegg.native = TRUE, key.pos = demo.paths$kpos1[i])
str(pv.out)
head(pv.out$plot.data.cpd)

#multiple states in one graph
set.seed(10)
sim.cpd.data2 = matrix(sample(sim.cpd.data, 18000, 
                              replace = TRUE), ncol = 6)
pv.out <- pathview(gene.data = gse16873.d[, 1:3], 
                   cpd.data = sim.cpd.data2[, 1:2], pathway.id = demo.paths$sel.paths[i], 
                   species = "hsa", out.suffix = "gse16873.cpd.3-2s", keys.align = "y", 
                   kegg.native = TRUE, match.data = FALSE, multi.state = TRUE, same.layer = TRUE)
str(pv.out)
head(pv.out$plot.data.cpd)




#?? Visualize the urea cycle pathway using gene expression data
pv.out <- pathview(gene.data = gse16873, pathway.id = urea_cycle_id, species = "hsa")