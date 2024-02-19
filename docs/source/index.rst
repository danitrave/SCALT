Welcome to SCALT's documentation!
=================================

**SCALT** (Single Cell Annotation Likelihood Tool) is an innovative tool which classifies cells from single cell RNA seq experiments at single cell resolution level taking the advantage of a maximum likelihood-based approach and a collection of pre-compiled cell type specific lists of genes constructed by extensive re-analysis of comprehensive and expert curated catalogues. 

The application does not require any **clustering**, **dimensionality reduction** or **manual annotation**.

For further details, we advise to read the paper available at **DOI**

As mentioned in the paper, SCALT is composed of other parallel utilities that allow to perform more complex analysis such as:

1. build cell type specific **lists of genes**  in a **deterministic** fashion starting from a **counts matrix** and the correspoding **annotation** for each cell;
2. build cell type specific **lists of genes** starting from a **count matrix** and a gathering of **user-defined** cell type specific lists of genes making the use of an **hypergeometric** test.

If you want to learn more about, please see the manual for point to point for instructions and tips.

.. toctree::
   :maxdepth: 2
   :caption: SCALT: Basic usage
   
   intro.rst

.. toctree::
   :maxdepth: 2
   :caption: Prerequisites
   
   prerequisites.rst

.. toctree::
   :maxdepth: 2
   :caption: Databases
   
   databases.rst

.. toctree::
   :maxdepth: 2
   :caption: SCALT: Classification
   
   classification.rst

.. toctree::
   :maxdepth: 2
   :caption: SCALT: Advanced
   
   advanced.rst

.. toctree::
   :maxdepth: 2
   :caption: SCALT: From Annotation
   
   annoListsBuilder.rst

.. toctree::
   :maxdepth: 2
   :caption: SCALT: From User-defiend Lists
   
   naiveListsBuilder.rst

.. toctree::
   :maxdepth: 2
   :caption: SCALT: For Brave People
   
   bravePeople.rst

.. toctree::
   :maxdepth: 2
   :caption: Auxiliary programs
   
   auxiliary.rst


