# EW_dmGWAS-Java
## Input file:
1. GWAS summary statistic file
2. Human protein-protein interaction network file

## 1. Infer node weight
p value to z score
## 2. Infer edge weight
normalizing probes
duplicate (aggregate)
first match with PPI
PCC

## 3. Run Java version of EW_dmGWAS on Windows
### 3.1 Command line
```
> c:
> mkdir C:\Users\fyan\EW_dmGWAS
> cd C:\Users\fyan\EW_dmGWAS
> java EW_dmGWAS_v3_1
usage: java EW_dmGWAS_v3_1 model={lambda|EW_dmGWAS}
        <node_weight_file=XX>
        <edge_weight_file=edge_weight_file>
        lambda=lambda
        <r=0.2>
        <output=output.modules.txt>
        <permutation={true|false}>
Example java EW_dmGWAS_v3_1 model=lambda node_weight_file=XX edge_weight_file=X
```
Example 1
```
java EW_dmGWAS_v3_1  model=lambda node_weight_file=XX edge_weight_file=XX 
```
The output looks like this:
```
input genes = XX
PPIs = XX edges, XX nodes
5.0     # modules=XX       avg_nodes=XX  avg_edges=XX estimated_lambda = XX
6.0     # modules=XX       avg_nodes=XX  avg_edges=XX estimated_lambda = XX
7.0     # modules=XX       avg_nodes=XX  avg_edges=XX estimated_lambda = XX
8.0     # modules=XX       avg_nodes=XX  avg_edges=XX estimated_lambda = XX
9.0     # modules=XX       avg_nodes=XX  avg_edges=XX estimated_lambda = XX
10.0     # modules=XX       avg_nodes=XX  avg_edges=XX estimated_lambda = XX
```

Example 2
```
java EW_dmGWAS_v3_1  model=EW_dmGWAS node_weight_file=XX edge_weight_file=XX lambda=1 r=0.2 output=XX.txt permutation=true
```
The output looks like this:
```
Start at: Fri 2019.04.05 at 11:22:45 AM CDT
Set lambda_1 = 1.0
input genes = XX
input genes = XX
PPIs = XX edges, XX nodes
10.0%...20.0%...30.0%...40.0%...9000.
50.0%...60.0%...70.0%...80.0%...90.0%...# find modules = XX
Start permutation for 1000 times. This may take some time.
0.1.2.3...997.998.999.Permutation finished !
Finished: Fri 2019.04.05 at 12:36:49 PM CDT
Module file: XX.txt
```
### 3.2 Output file
Output format

| module_genes | seed | module_score | edge_score| n_edges | z_perm |
| ---------- | ---------- |  :----:  |  :----:  |  :----:  |  :----:  |
| gene1 gene2 | gene1 | XX | XX | 1 | XX |
| gene2 gene3 gene4 | gene2 | XX | XX | 2 | XX |
## 4. Export the result to Cytoscape and Build Network


