We  are writing several functions that will aletmiately be used for a pacakge , the package is dealing with locus boundary detection in genomic data.

We have FUMA identification, and loci reported in previous studies.
The tool will take as input the
GWAS summary statistics file, and the FUMA identified loci, and the previous reported loci.
Additionally, the plink reference panel file path is also provided as input. Will be used to calculate LD structure and visualization of LD structure.


If the input FUMA identified locus is crossed with three loci in the previuosly identified loci, we will break the locus in the three loci , where using the start and end of the biggest so the first new locus will have the first start, and the last new locus will have the bigger end (if the FUMA locus is biggere it will take it if the prevoiuos it will take it )

So The idea is to split the FUMA locus to the previuosl know independt loci and and to improve the start and end (not losing the start and end of the LOCUS taking the best)

We can have multiple FUMA we will algin them , so we can merge the intrecrossed loci , the ouput and can then be used with the previuolsy identified to be splitted

All the SNPs higher than 5E-8, will be extracted along with the 1000 around each siginficant SNP then usinque these SNPs will be extracted and the LD will be caluclated to idnetify how many idnentpent locus within each have a SNP with significant threshold