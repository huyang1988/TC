# TC


## Preprocess from SNAP data

### converter
```
TC/graph_converter/undirected_csr/
```

```
./text_to_bin input
```
It automatically generates CSR files ```adjacent.bin```, ```head.bin``` and ```begin.bin``` to current path. Each line in SNAP file (e.g., ```0,1```) generates two edges (0,1) and (1,0). The converted files are bit files. The begin position file ```begin.bin``` represents the offsets of each neighbor lists in both ```head.bin``` and ```adjacent.bin```. The neighbor lists are not sorted.

### cleaner
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

## Triangle counting

### GPU

### CPU

### GPU on 2d partitioned data

### CPU on 2d partitioned data


