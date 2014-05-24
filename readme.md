#FSpw
A tool to perform statistical power analyses for FreeSurfer data
--

See Liem et al. (submitted). Reliability and power analysis of cortical and subcortical FreeSurfer metrics.

In this repository you will find
- 1) surface data of the reliability analysis.
- 2) the matlab tool to perform different kinds of statistical power analyses. 


1) Reliabilty
----
The file reliabilty.zip contains ICC overlay files in fsaverage space

2) Power
----
**Dependencies**
- FreeSurfer's save_mgh (in FreeSurfer's matlab directory)


**Input**

To perform statistical power analyses a configuration file has to be created.
In this file a *cfg* structure is provided with all the necessary information to compute power analyses. This *cfg*
structure is then passed on to the FSpwr.m function.

For an example see **FSpwr_config.m**

**Output**

Data is written to the directory given in *cfg.outputPath*. Surfaces are in fsaverage space.
