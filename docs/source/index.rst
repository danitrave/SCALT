Welcome to SCARLET's documentation!
=================================

**SCARLET** (Single Cell Annotation Likelihood Tool) is an innovative tool which classifies cells from single cell RNA seq experiments at single cell resolution level combining a likelihood-based approach and a collection of pre-compiled cell type specific lists of genes constructed through an extensive re-analysis of comprehensive and expert curated catalogues. 

The tool does not require any **clustering**, **dimensionality reduction** or **manual annotation** and cells are annotated at single cell level.

Beyond clasification, SCARLET arranges additional utilities for eventually handling problems related to the analysis of single cell RNA sequencing data found on a daily basis:

1. build estensive cell type specific **lists of genes** in a **deterministic** fashion starting from a **counts matrix** and the correspoding **annotation** for each cell;
2. ascertain the **similarity** and **reliability** of the classification of the cell types used to generate the lists and eventually suggesting an **operative definition** in case of highly similar annotations.

For further details, we advise to read the paper available at **DOI**

If you want to learn more about, please see the manual for point to point for instructions and tips.

.. toctree::
   :maxdepth: 2
   :caption: SCARLET: Basic usage
   
   intro.rst

.. toctree::
   :maxdepth: 2
   :caption: Prerequisites
   
   prerequisites.rst

.. toctree::
   :maxdepth: 2
   :caption: SCARLET: Classification
   
   classification.rst

.. toctree::
   :maxdepth: 2
   :caption: SCARLET: Case Study
   
   tutorial_pbmc.rst

.. toctree::
   :maxdepth: 2
   :caption: SCARLET: Advanced
   
   advanced.rst

.. toctree::
   :maxdepth: 2
   :caption: SCARLET: From Annotation
   
   annoListsBuilder.rst

.. toctree::
   :maxdepth: 2
   :caption: SCARLET: Mutual Information
   
   mutual_info.rst

.. toctree::
   :maxdepth: 2
   :caption: SCARLET: For Brave People
   
   bravePeople.rst

.. toctree::
   :maxdepth: 2
   :caption: Auxiliary programs
   
   auxiliary.rst


