# DendriteArborAnalyzer
DendriteArborAnalyzer is an ImageJ plugin for the collection of data from images of neuronal dendritic arbors. Data collected include branch lengths, path distances, and Strahler analyses of branching structure. 

This plugin extends the functionality of the (Strahler Analysis plugin)[https://imagej.net/plugins/strahler-analysis].

## Usage
The input to the program should be a skeletonization (see (Skeletonize3D)[https://imagej.net/plugins/skeletonize3d]) of the neuronal dendritic arbor with a rectangular ROI around where the soma is on the skeleton. 
1. Open the image, subject of analysis.
2. Run plugin. Plugin can be found through ImageJ command menu: Analyze > Skeleton > DendriteArborAnalyzer2.
3. Select parameters of choice (described below).
4. Click OK.

### Data/Images Produced
* Strahler Table containing statistics on the number, ramification ratios, and average branch length of all branches with each Strahler Order.
* Branch Table containing data from each branch of the dendrite including Strahler order, length, path distance to soma, euclidean distance to soma, and parent order.
* Strahler Mask, a copy of input image with branches colored by Strahler order.
* BranchID Mask, Strahler Mask with branches labeled with their branch # in the Branch Table.
* Frequency Tables binning branches of a given order based on the measurements: path distance, parent order, euclidean distance, branch length.

### Parameters
Yes/No Parameters:
* Display Tables - if yes, generates Strahler Table and Branch Table.
* Display Strahler Mask - if yes, generates Strahler Mask.
* Display Branch ID - if yes, generates BranchID Mask.
* Reverse Branch Order - Numbers branch orders in reverse, with larger order branches being closer to the tip.
* Display Frequency Table - if yes, generates Frequency Tables.

Variable Parameters:
* Branch Length Bin Size - Bin size used in producing the branch length frequency table.
* Path Distance Bin Size - Bin size used in producing the path distance frequency table.
* Maximum Length to be Higher Order - Branches above this length will not be classified as higher order (see below).
* \# Orders to be higher Order (0 for no threshold) - Specifies the number of orders considered higher order
Note: the last two parameters allow for the functionality to distinguish long terminal branches from short terminal branches in assigning order. When the # orders > 0, branches above the specified length will be assigned to the maximum of their Strahler order and the # of higher order branches, preventing long branches from having order below the specified number. This may be preferable for comparing phenotypes with pruning of short terminal branches.

## Installation
1. Download the plugin found [here](target/DendriteArborAnalyzer__-1.0.13-SNAPSHOT.jar).
2. Move the plugin to your local Imagej `Plugins` folder as described [here](https://imagej.net/plugins/#:~:text=Advanced%20topics-,Installing%20plugins%20manually,-If%20the%20plugin).
3. Restart ImageJ.