Welcome to SCALT's documentation!
=================================

**SCALT** (Single Cell Annotation Likelihood Tool) is an innovative tool which classifies cells from single cell RNA seq experiments at single cell resolution level combining a likelihood-based approach and a collection of pre-compiled cell type specific lists of genes constructed through an extensive re-analysis of comprehensive and expert curated catalogues. 

The tool does not require any **clustering**, **dimensionality reduction** or **manual annotation** and cells are annotated at single cell level.

Beyond clasification, SCALT arranges additional utilities for eventually handling problems related to the analysis of single cell RNA sequencing data found on a daily:

1. build estensive cell type specific **lists of genes** in a **deterministic** fashion starting from a **counts matrix** and the correspoding **annotation** for each cell;
2. ascertain the **similarity** and **reliability** of the classification of the cell types used to generate the lists and eventually suggesting an **operative definition** in case of highly similar annotations.

For further details, we advise to read the paper available at **DOI**

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
   :caption: SCALT: Similarity through Mutual Information
   
   mutual_info.rst

.. toctree::
   :maxdepth: 2
   :caption: SCALT: For Brave People
   
   bravePeople.rst

.. toctree::
   :maxdepth: 2
   :caption: Auxiliary programs
   
   auxiliary.rst


