TC
===========================

## Preprocess from SNAP data

### Data Formats

#### SNAP
A SNAP dataset http://snap.stanford.edu/ is a text file that each line is an edge with source and destination vertex IDs, such as 
  ```11  25```

#### CSR

### converter

The codes are in the folder below.
```
TC/graph_converter/undirected_csr/
```

```
./text_to_bin input
```
It automatically generates CSR files ```adjacent.bin```, ```head.bin``` and ```begin.bin``` to current path. Each line in SNAP file (e.g., ```0,1```) generates two edges (0,1) and (1,0). The converted files are bit files. The begin position file ```begin.bin``` represents the offsets of each neighbor lists in both ```head.bin``` and ```adjacent.bin```. The neighbor lists are not sorted.

### cleaner

The codes are in the folder below.
```
TC/graph_cleaner/
```


```
./cleaner input_path
```
The cleaner sorts all neighbor lists, removes self-loop and multiple edges. After this step, the graph is an undirected graph in CSR format (both directions of each edge is stored in the CSR). The output files are written to current folder, so to copy the executable to target folder is suggested.

### rank-by-degree

```
./rank-by-degree input_path
```
The rank-by-degree orientation tool removes half of all edges. It only keeps the edges from lower degree vertex to higher degree vertex. After this step, the graph is an undirected graph in CSR format. The output files are written to current folder, so to copy the executable to target folder is suggested.


## Partition

The codes are in the folder below.
```
TC/TC-EXT/partitioner-x/
```


```
./partition-x input_path
```
The partitioner tool partition the graph in a 2-D fashion.
The rule is described in paper 
It only keeps the edges from lower degree vertex to higher degree vertex. After this step, the graph is an undirected graph in CSR format. The output files are written to current folder, so to copy the executable to target folder is suggested.


## Triangle counting

### GPU in-memory
The path of in-memory GPU triangle counting is ```TC/TC-GPU/work-steal/```.
Run the code with command ```./tc <input_path>```.

### CPU in-memory
The path of in-memory CPU triangle counting is ```TC/TC-CPU/tc-ne-cpu/```.
Run the code with command ```./tc <input_path>```.

### GPU on 2d partitioned data

### CPU on 2d partitioned data

## Graph Challenge Example Dataset and Toolkit

The dataset path is ```TC/example_dataset/amazon0302/```, the _\_adj.mmio_ is the original data format downloaded from https://graphchallenge.mit.edu/data-sets. Converter tool source code is in ```TC/gConv/```, copy the executable to dataset path, line 8 of the bash script ```converter.sh``` needs to be modified (set the number to (#vertex+1), the vertex count can be found by ```less <input>_adj.mmio```). Run the bash script with command ```./converter.sh <input>_adj.mmio``` (replace the input name, for example, amazon0302).

