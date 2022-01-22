# HGCAL_analysis_code

## Neutral Pion analysis code for OctTB2018 (33m bealine)


Any CMSSW version environment should be fine (do `cmsenv`).

### Downlaod scripts <br/>
`git clone -b  check_pi0_analysis git@github.com:spandeyehep/HGCTB_Oct18_pion_analysis.git .` <br/>

### How to run the script for data and simulation: <br/>

`cd scripts`<br/>
`make`<br/>

`./analyzeHGCOctTB <file_list> outFileName.root <dataset> <energy>`<br/>


`<file_list.txt>` : contains files to be analyzed (can be found in file_list folder)<br/>
`<dataset>` : `FTFP_BERT_EMN` or `QGSP_FTFP_BERT_EMN` <br/>
`<energy>` : beam_energy<br/>
Example: <br/>
`./analyzeHGCOctTB file_list/33m_1060pre4X/FTFP_BERT_EMN/pion_100.txt out.root FTFP_BERT_EMN 100`

### Descrition of scripts: <br/>
`AnalyzeHGCOctTB.cc` : Main analysis code <br/>
`AnalyzeHGCOctTB.h` : Initialize histos here <br/>
`HGCNtupleVariables.h`: Tree variable initialization <br/>
`txt_maps/`: Dependencies  <br/>
