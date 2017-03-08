Noparama
========

This software is about mathematics, in particular nonparametric Bayesian models.

Installation
------------

On Ubuntu:

.. code-block:: bash

	sudo apt install git cmake gcc libeigen3-dev
	echo "Go to your preferred installation directory"
	git clone git@github.com:mrquincle/noparama.git
	cd noparama
	make
	echo "The compiled binary you can find in build/bin/

A dataset can be found in a companion repository: https://github.com/mrquincle/noparama-datasets.git.

Use
---

After compiling the program it can be run as:

.. code-block:: bash

	./noparama

--data file                  file with data to process, e.g. datasets/twogaussians.data  
--config file                file with configuration values

Example:

.. code-block:: bash

	cd build/bin
	./noparama -d ../datasets/twogaussians.data

The datasets are not included with the repository, but in https://github.com/mrquincle/noparama-datasets.git.

Structure
---------

The elements in this software package or often vectors and matrices, hence Eigen is used.

Density functions
-----------------

Density functions are represented by their sufficient statistics. A normal probability distribution is represented by 
mean and covariance matrix.

Data
----

The data is clustered in sets. In MCMC we often loop through all data items, first removing a data item from the 
cluster it is currently assigned to, then assigning it to an existing or new cluster.

The clusters themselves are also updated where we iterate over all clusters and update a cluster given all data items 
assigned to it.

A hashmap for the data items allows us to uniquely identify each data item. A cluster refers to its data items through 
a vector of keys. Removing a data item from a cluster is by removing a key from this vector.

A cluster can also refer to its data items through a pointer to the data item including its identifier. Search is 
still easy by the use of keys. Removal is by deleting a pointer from a vector or set.

To remove a data item from a cluster, we will have to iterate through all clusters, or also store the membership 
information per data item. To also store a reference to a cluster per data item is cumbersome and leads to double 
updates.

Requirements:

* update membership information in O(1) and in one atomic operation
* iterate over data items in cluster X in O(D), with D number of data items that are a member of cluster X
* iterate over all data points in O(N), with N number of data items 

This data structure is a matrix with columns with only a single 1.

===  ===  ===  ===  ===  ===  ===
\    D0   D1   D2   D3   ...  DN
===  ===  ===  ===  ===  ===  ===
C0    1    0    0    0
C1    0    1    1    0
C2    0    0    0    0
...
...
CK
===  ===  ===  ===  ===  ===  ===

A data item that is not assigned is represented by a zero column-vector, a cluster without data points by a zero 
row-vector. C0, C1, etc. refers to an object with cluster parameters. D0, D1, etc. refer to an object with data values.

Hence, we need a matrix plus vectors to store references to the mentioned objects.

Then a function like matrix::assign(cluster, data) to be fast must be with indices into those vectors. If this is not 
the case we need to search through those vectors.

The indices can be stored in an encapsulating object, e.g. as a `std::pair<index, object>`, but that might not be 
necessary. It is also preferably to return data without capsulating structures. 

If we want to return a container with data items, we might also store the data items directly.

===  ===  ===  ===  ===
C0   D0   0    0    0
C1   0    D1   D2   0
C2   0    0    0    0
...
...
CK
===  ===  ===  ===  ===

This however, would still require us to create a set out of something like [0 D1 D2 0] @C1.

Hence, what we can do is to maintain two data structures. A matrix structure:

===  ===  ===  ===  ===  ===  ===
\    D0   D1   D2   D3   ...  DN
===  ===  ===  ===  ===  ===  ===
C0    1    0    0    0
C1    0    1    1    0
C2    0    0    0    0
...
...
CK
===  ===  ===  ===  ===  ===  ===

Plus a set structure:

===  ===  ===  ===  ===  ===  ===
C0   D0
C1   D1   D2
C2
...
...
CK
===  ===  ===  ===  ===  ===  ===

Here we do not have the property anymore that the update is atomic! Setting something to 1 or 0 in the assignment 
matrix, needs also an update in the set structure.

In Progress
-----------
   
This is a work-in-progress. This means it normally can be compiled and should be functionality complete.
However, don't use it yet.

:Authors:
    Anne van Rossum

:Version: 0.1.53
