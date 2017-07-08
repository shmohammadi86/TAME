Triangular Alignment (TAME) is a new formulation for the network alignment problem that uses network motifs as the driving factors to align a pair of given graphs. In order to install TAME, simply download the latest version and run:

```
make
```

To run TAME, the bare minimum configuraion is as follows:

```
tri-match -t <type> -G <1st_net> -H <2nd_net>
```

where the network files are in tab separated format by default. In this sense, each file contains lines in the form: "src_node dst_node." Additional parameters can be provided to config different aspects of TAME. To get a comprehensive list of parameters, type:

```
./tri-match --help
```


There are three example scripts to help you get a sense of how to start using TAME:

1. *sparse_run.sh*: aligns yeast vs human interactome with sequence similarity in sparse mode.
2. *dense_run.sh*: aligns yeast vs human interactome with sequence similarity in dense mode.
3. *dense_run_noPrior.sh*: aligns yeast vs human interactome without sequence similarity in dense mode.

# Cite
Mohammadi, S., Gleich, D. F., Kolda, T. G., & Grama, A. (2016). Triangular Alignment (TAME): A Tensor-based Approach for Higher-order Network Alignment. IEEE/ACM Transactions on Computational Biology and Bioinformatics (http://doi.org/10.1109/TCBB.2016.2595583)
