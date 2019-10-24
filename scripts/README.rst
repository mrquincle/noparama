Different scripts
=================

There are many scripts that help with running a complete trial or visualize the results.

Example
-------

We can run the binary itself in the following manner... 

.. code-block:: bash

    build/bin/noparama -d ../data/lines/dpmline4.pnts.Sat_Aug_29_19:28:36_2015.811530.data.txt -c regression -T 10000 -a algorithm8

This will however only run one particular item. 

Run entire batch
----------------

The batch with all data frames can be run at once: 

.. code-block:: bash

    scripts/run_all.sh ../data/lines algorithm8

This will write down the results in output/algorithm8. Each run creates a separate directory. If this directory already exists, the particular run will be aborted, so no results will be accidentally overwritten.

Visualize the results
---------------------

Collect all the results

.. code-block:: bash

    cd scripts
    ./collect.sh ../output/algorithm8

And now create a directory and visualize the results:

.. code-block:: bash

    mkdir -p local
    Rscript visualize.R algorithm8 ../output/algorithm8 local

Visualize an individual data frame
----------------------------------

The script is called `visualize_input.R` and it only works on lines or line segments:

.. code-block:: bash

    mkdir -p temp
    Rscript visualize_input.R ../../data/lines/dpmline4.pnts.Sat_Aug_29_19:28:33_2015.811501.data.txt temp

Visualize a particular result 
-----------------------------

The script is called `visualize_single_fit.R` and it only works on lines or line segments:

.. code-block:: bash

    mkdir -p temp
    Rscript visualize_single_fit.R ../output/algorithm8/dpmline4.pnts.Sat_Aug_29_19:28:36_2015.811530.data.txt/LATEST temp

Analysis of results
-------------------

Analyze the results, mean, median, etc.:

.. code-block:: bash

    ./analyze.m ../output/algorithm8

Note that the violin plot shows the median values very clearly.

Other scripts
-------------

Check also the scripts at the repository: emd_suite_

.. _emd_suite: https://github.com/mrquincle/emd-suite/tree/master/analyse

It contains R scripts to visualize cubes for example.
