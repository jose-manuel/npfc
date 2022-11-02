======================
Command Line Interface
======================

The annotation of PNPs require the processing of three different datasets: fragments, natural and synthetic.
To better organize and orchestrate the work, snakemake workflows were developped.

Each workflow is based on a set of Python or Bash scripts, that can be used individually as well in fragments, natural and synthetic workflows.

Current scripts are regrouped in different categories and described in the table below:


.. csv-table:: List of commands
   :header: "Category","Name","Description"

   "Core","fc_classify","Classify fragment combinations by pairwise comparison of fragment hits for each molecule."
   "Core","fcg_annotate_pnp","Annotate FCGs with PNP information."
   "Core","fcg_generate","Generate fragment combination graphs using all fragment combinations found for each molecules."
   "Core","frags_annotate_fcp","Annotate the fragments with Fragment Combination Points (FCP)."
   "Core","mols_dedupl","Filter out duplicate molecules and (if provided) update a reference file for processing in parallel."
   "Core","mols_depict","Generate 2D-coordinates for molecules."
   "Core","mols_draw","Draw molecules into individual files. Useful for exporting fragment images for FCG representation."
   "Core","mols_extract_murcko","Extract the Murcko scaffold from molecules."
   "Core","mols_load","Load molecules into RDKit format."
   "Core","mols_standardize","Apply a serie of standardization operations on molecules."
   "Core","mols_subset","Filter out molecules from a molecular file using another molecular file as a reference."
   "Miscelleneous","fhits_filter","Filter out one or several fragment hits from fragment- search and combination results using their fragment ID."
   "Report","report/report_fcg_chunk","Preprocess the report of the FCG information for a given chunk  (optimization for cluster)."
   "Report","report/report_fcg_concat","Concatenate all chunks of FCG informations produced by report_fcg_chunk and generate the report."
   "Report","report/report_fcg2","Report the number of fragment hits and fragment combinations and fragment combination graphs."
   "Report","report/report_fcp","Report the number of symmetry centers found in fragments." 
   "Report","report/report_mols_count","Report the number of molecules found for a given chunk in all steps of a workflow."
   "Report","report/report_pnp","Report the number of PNP and NP-like molecules."
   "Report","report/report_prep","Report the number of passed filtered and error molecules during preparation."
   "Report","report/report_subset","Report the number of molecules filtered out of the input dataset during subsetting."
   "Report","report/report_time","Report  the computation duration of a given chunk in all steps of a workflow."
   "Workflow management","chunk_check","Check if there are missing chunks in the output directory."
   "Workflow management","chunk_sdf","Split a SDF into chunks without parsing molecules."
   "Workflow management","concat_sdf","Concatenate SDF chunks into a single SDF without parsing molecules."
   "Workflow management","concat_synonyms","Concatenate all synonym files (references for identifying duplicate molecules) into a single file."
   "Workflow management","mols_concat","Concatenate two molecular files."
   "Workflow management","mols_count","Count the number of molecules in a molecular file."
   "Workflow management","report_protocol","DEPRECATED. Reports are now generated automatically during workflow executions."
   "Workflow management","run_protocol_fc","Run a snakemake workflow (either fragments or natural or synthetic)."

Some commands are used in the workflows (Core and Report categories), whereas some others were used only once but were kept in the tool for convenience.
One example would be the the fhits_filter script, which was used to filter out benzene fhits from the fragment search and fragment combination results,
without having to restart these computations.