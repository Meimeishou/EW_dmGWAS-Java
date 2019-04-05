# EW_dmGWAS-Java
## Infer node weight
p value to z score
## Infer edge weight
normalizing probes
duplicate (aggregate)
first match with PPI
PCC

## Java version of EW_dmGWAS on Windows
### Command line
c:
cd C:\Users\fyan\EW_dmGWAS
C:\Users\fyan\EW_dmGWAS> java EW_dmGWAS_v3_1
usage: java EW_dmGWAS_v3_1 model={lambda|EW_dmGWAS}
        <node_weight_file=XX>
        <edge_weight_file=edge_weight_file>
        lambda=lambda
        <r=0.2>
        <output=output.modules.txt>
        <permutation={true|false}>
Example java EW_dmGWAS_v3_1 model=lambda node_weight_file=XX edge_weight_file=XX
java EW_dmGWAS_v3_1  model=0.5 node_weight_file=XX edge_weight_file=XX r=0.2 output=XX  permutation=true
java EW_dmGWAS_v3_1  model=EW_dmGWAS node_weight_file=Asian_node_weight_scale.txt edge_weight_file=edge_weight_scale.txt lambda=1 r=0.2 output=Asian.lambda_1_permutated.txt permutation=true
java EW_dmGWAS_v3_1  model=lambda node_weight_file=Asian_node_weight_scale.txt edge_weight_file=edge_weight_scale.txt r=0.2 output=test.txt 
### Output file
## Export the result to Cytoscape and Build Network