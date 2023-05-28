# Time-Series Similarity-Search (Nearest Neighbor Search) and Clustering
This project implements two fundamental algorithmic aspects for time-series handling.
* approximation algorithms for similarity search on time-series.
* clustering algorithms on time-series.

## Part 1 : Similarity-Search Algorithms:
**Locality Sensitive Hashing (LSH) on time-series**\
The time-series are represented as:
* euclidean polygonal curves, and similarity is measured by implementing the Discrete Frechet Metric.
* one dimensional polygonal curves, after using dimensionality and complexity reduction techniques (filtering/ Grid Snapping). Similarity in this case is measured using the Continuous Frechet Metric.


## Part 2 : Clustering Algorithms:
The implemented clustering algorithms follow the general principle of centroid-based clustering. The clusters are defined by a set of centroids, each time-series is assigned to its (exact or approximate) closest centroid. So the clustering problem is reduced to the initialization and iterative improvement of the centroids, with respect to the input time-series. Various approaches are implemented:
### Centroid Initialization:
* k-means++ initialization technique
### Time-Series to Centroid Assignment:
* Lloyd's exact nearest centroid algorithm (simplest approach)
* Range search for approximate nearest centroid (using LSH or Hypercube random projection algorithms on time-series)

### Update:
* Mean Curve method

The Silhouette metric is used for the evaluation of the quality of the clustering.

## Part 3: Unit Tests are implemented

**The project was developed as part of the "Software Development for Hard Algorithmic Problems" course of the Computer Science Department of the university of Athens (NKUA)**
