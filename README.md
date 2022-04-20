# Microanatomy of the human atherosclerotic plaque by single-cell transcriptomics

_Authors_: Marie A.C. Depuydt, Koen H.M. Prange, Lotte Slenders, Tiit Örd, Danny Elbersen, Arjan Boltjes, Saskia C.A. de Jager, Folkert W. Asselbergs, Gert Jan de Borst, Einari Aavik, Tapio Lönnberg, Esther Lutgens, Christopher K. Glass, Hester M. den Ruijter, Minna U Kaikkonen, Ilze Bot, Bram Slütter, Sander W. van der Laan, Seppo Yla-Herttuala, Michal Mokry, Johan Kuiper, Menno P.J. de Winther, Gerard Pasterkamp.


This repository holds all scripts used to analyze data and create figures presented in ["Microanatomy of the human atherosclerotic plaque by single-cell transcriptomics"]().


> # Abstract
> 
> **Rationale:** Atherosclerotic lesions are known for their cellular heterogeneity, yet the molecular complexity within the cells of human plaques have not been fully assessed. 

> **Objective:** Using single-cell transcriptomics and chromatin accessibility we gained a better understanding of the pathophysiology underlying human atherosclerosis.

> **Methods and Results:** We performed single-cell RNA and single-cell ATAC sequencing on human carotid atherosclerotic plaques to define the cells at play and determine their transcriptomic and epigenomic characteristics. We identified 14 distinct cell populations including endothelial cells, smooth muscle cells, mast cells, B cells, myeloid cells, and T cells and identified multiple cellular activation states and suggested cellular interconversions. Within the endothelial cell population we defined subsets with angiogenic capacity plus clear signs of endothelial to mesenchymal transition. CD4+ and CD8+ T cells showed activation-based subclasses, each with a gradual decline from a cytotoxic to a more quiescent phenotype. Myeloid cells included two populations of pro-inflammatory macrophages showing IL1B or TNF expression as well as a foam cell-like population expressing TREM2 and displaying a fibrosis-promoting phenotype. ATACseq data identified specific transcription factors associated with the myeloid subpopulation and T cell cytokine profiles underlying mutual activation between both cell types. Finally, cardiovascular disease susceptibility genes identified using public GWAS data were particularly enriched in lesional macrophages, endothelial and smooth muscle cells. 

> **Conclusion:** This study provides a transcriptome-based cellular landscape of human atherosclerotic plaques and highlights cellular plasticity and intercellular communication at the site of disease. This detailed definition of cell communities at play in atherosclerosis will facilitate cell-based mapping of novel interventional targets with direct functional relevance for the treatment of human disease.


**Figure 1: Overall expression of some putative cardiovascular 'target' genes in carotid plaques from the Athero-Express Biobank Study**
![Overall expression of some putative cardiovascular 'target' genes in carotid plaques from the Athero-Express Biobank Study](FIGURES/20200625.TargetExpression_vs_1000genes.png)



# Scripts
We shared the scripts we used for the analyses of the data. These are provided as-is: they were not tested on different systems, _et cetera_.

* `SCRIPTS` : this folder contains several scripts used for the analyses.


# Data access
The data, _i.e._ the single-cell RNAseq and associated experimental/clinical data from the [Athero-Express Biobank Study](http://www.atheroexpress.nl), used for this article are submitted to [DataverseNL and available through here](https://doi.org/10.34894/RWAHLS). 


--------------

#### The MIT License (MIT)
##### Copyright (c) 1979-2021

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:   

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Reference: http://opensource.org.



