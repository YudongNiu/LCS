# Discover Your Characteristic Community from Graph Hierarchical Communities

Build HIMOR-Index
-------
To build the HIMOR-Index for our COD problem, execute the following command on linux:

```sh
python3 himor.py dataset
```

The parameter "dataset" denotes the name of your network.

Running code
-------
To run the code for COD problem, execute the following command on linux:

```sh
python3 COD.py dataset labeled K cpr loc indexed
```

**There are 6 parameters**
* dataset: the name of your network
* labeled: 
* K: the 
* cpr: if the compressed COD framework is used
* loc: if the lore algorithm is used
* indexed: if the HIMOR-index is used

For example, the following command

```sh
python3 COD.py cora 1 5 1 1 1
```
performs 

Input Files
-----------
The COD.py requires three input files.