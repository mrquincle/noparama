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

A dataset can be found in a companion repository: dataset_

.. _dataset: https://github.com/mrquincle/noparama-datasets.git.

This is a Monte Carlo approach, so you will have to take additional care that there are good random numbers available
on your system. The ``entropy_avail`` file has typically a value higher than 3000 on my system. 

.. code-block:: bash

	cat /proc/sys/kernel/random/entropy_avail

You can remediate this by installing ``rng-tools``.

This simulation has been run with g++ version 6.2.0, 7.0.1, 7.1.0 without trouble. However, running it with 5.4.0 
lead to significantly slower convergence. There might be some errors in older g++ compilers that reduce the diversity
in random numbers. I'm too lazy to figure that out, just use a more modern compiler if you encounter this.

Use
---

After compiling the program it can be run as:

.. code-block:: bash

	./noparama

--data file                  file with data to process, e.g. datasets/twogaussians.data  
--config file                file with configuration values

Example:

.. code-block:: bash

	./build/bin/noparama -d ../datasets/twogaussians.data -a triadic -T 1000 -c clustering

Typically, the purity should be almost 1, the rand index a bit lower and the adjusted rand index a bit lower again.

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

Literature
----------

Currently implemented are algorithms 2 and 8 by Neal, the split-merge sampler with the simple random split procedure
by Jain and Neal, and the nonconjugate SAMS sampler by Dahl. The Smart-Dumb/Dumb-Smart sampler just as the Triadic
sampler makes the acceptance between splits/merges more balanced. 

The Informed Sampler uses a mixture of a global Metropolized independence sampler and a local (normal) Metropolis
sampler. Brilliant! We'll definitely implement a variant on this later.

1. Markov chain sampling methods for Dirichlet process mixture models (`Neal, 2000`_).
2. A split-merge Markov chain Monte Carlo procedure for the Dirichlet process mixture model (`Jain, Neal, 2004`_).
3. An improved merge-split sampler for conjugate Dirichlet process mixture models (`Dahl, 2003`_).
4. Sequentially-allocated merge-split sampler for conjugate and nonconjugate Dirichlet process mixture models (`Dahl, 2005`_).
5. A Smart-Dumb/Dumb-Smart Algorithm for Efficient Split-Merge MCMC (`Wang, Russell, 2015`_).
6. The Informed Sampler: A Discriminative Approach to Bayesian Inference in Computer Vision (`Jampani, Nowozin, Loper, Gehler, 2015`_)

.. _Neal, 2000: https://pdfs.semanticscholar.org/de79/8ab2f2e7ca312c12ba34a0d9c05cff9fbf3c.pdf
.. _Jain, Neal, 2004: https://pdfs.semanticscholar.org/6305/dcc03c8378e371e73b0a68ff29f1167a65f0.pdf
.. _Dahl, 2003: https://pdfs.semanticscholar.org/cfbe/dd12e5040c76e7c9981b19e6e333d6111656.pdf
.. _Dahl, 2005: https://pdfs.semanticscholar.org/f49c/620fa006d2e1e07c71092d9692ba5d71f14f.pdf
.. _Wang, Russell, 2015: https://pdfs.semanticscholar.org/c444/ee208269dcfe6c96ede88525893549c39add.pdf
.. _Jampani, Nowozin, Loper, Gehler, 2015: https://arxiv.org/pdf/1402.0859.pdf

In Progress
-----------
   
This is a work-in-progress. This means it normally can be compiled and should be functionality complete.

Current state:

* The split-merge algorithms etc have been tested extensively.
* The simpler algorithms were first implemented using matlab/octave and are now ported to C++. This is nontrivial because they assume conjugacy. 
* That means that now I'm implementing update and downdate of sufficient statistics. For example the NIW distribution
  should have an update(data_t data) or even update(dataset_t dataset) function that at once adjusts its values.


:Authors:
    Anne van Rossum

:Version: 0.1.78
