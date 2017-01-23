Triangular Alignment (TAME) is a new formulation for the network alignment problem that uses network motifs as the driving factors to align a pair of given graphs. In order to install TAME, simply download the latest version and run:

```
make
```

To run TAME, the bare minimum configuraion is as follows:

```
tri-match -G <1st_net> -H <2nd_net>
```

where the network files are in tab separated format by default. In this sense, each file contains lines in the form: "src_node dst_node." Additional parameters can be provided to config different aspects of TAME. To get a comprehensive list of parameters, type:

```
./tri-match --help
```


There are three example scripts to help you get a sense of how to start using TAME:

1. *mini_run.sh*: runs a toy example as a proof of concept that code is compiled and working correctly. Visualization of input graphs is available as *.eps files under input/ folder.
2. *run_Family1.sh*: runs TAME on an example network from the NAPABench
3. *sparse_run.sh*: runs sparse formulation of TAME on the yeast and human interactome downloaded from the BioGRID, without any postprocessing, which is the fastest configuration for large networks. This can serve as a good initial solution for other alignment/postprocessing methods, or as a final alignment itself. It is notable here that postprocessing has a --significant-- effect on improving the quality of alignments.

# Cite
Mohammadi, S., Gleich, D. F., Kolda, T. G., & Grama, A. (2016). Triangular Alignment (TAME): A Tensor-based Approach for Higher-order Network Alignment. IEEE/ACM Transactions on Computational Biology and Bioinformatics (http://doi.org/10.1109/TCBB.2016.2595583)
