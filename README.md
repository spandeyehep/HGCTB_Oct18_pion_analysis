# HGCAL_analysis_code

## Pion analysis code for OctTB2018


Any CMSSW version environment should be fine (do `cmsenv`).

### Downlaod scripts <br/>
`git clone git@github.com:spandeyehep/HGCTB_Oct18_pion_analysis.git .` <br/>

### How to run the script for data and simulation: <br/>

`cd scripts`<br/>
`make`<br/>

`./analyzeHGCOctTB <file_list> outFileName.root <dataset> <energy>`<br/>


`<file_list.txt>` : contains files to be analyzed (can be found in file_list folder)<br/>
`<dataset>` : `data` or `FTFP_BERT_EMN` or `QGSP_FTFP_BERT_EMN` <br/>
`<energy>` : beam_energy<br/>
Example: <br/>
`./analyzeHGCOctTB file_list/data/pion_120.txt out.root data 120`

### Descrition of scripts: <br/>
`AnalyzeHGCOctTB.cc` : Main analysis code <br/>
`AnalyzeHGCOctTB.h` : Initialize histos here <br/>
`HGCNtupleVariables.h`: Tree variable initialization <br/>
`txt_maps/`: Dependencies  <br/>

### Important flags: <br/>
`AnalyzeHGCOctTB.cc` has following three flags in `main(...)` function:<br/> <br/>
`shower`: Set `true` for shower energy histograms <br/>
`longi` : Set `true` for longitudinal shower profiles <br/>
`trans` : Set `true` for transverse shower profiles <br/>

<br/>
It is recommeded that only flag is set to `true` for minimum resource consumption.
