#ifndef AnalyzeHGCOctTB_H
#define AnalyzeHGCOctTB_H

//the following are UBUNTU/LINUX, and MacOS ONLY terminal color codes.
#define RESET   "\033[0m"
#define BLACK   "\033[30m"      /* Black */
#define RED     "\033[31m"      /* Red */
#define GREEN   "\033[32m"      /* Green */
#define YELLOW  "\033[33m"      /* Yellow */
#define BLUE    "\033[34m"      /* Blue */
#define MAGENTA "\033[35m"      /* Magenta */
#define CYAN    "\033[36m"      /* Cyan */
#define WHITE   "\033[37m"      /* White */
#define BOLDBLACK   "\033[1m\033[30m"      /* Bold Black */
#define BOLDRED     "\033[1m\033[31m"      /* Bold Red */
#define BOLDGREEN   "\033[1m\033[32m"      /* Bold Green */
#define BOLDYELLOW  "\033[1m\033[33m"      /* Bold Yellow */
#define BOLDBLUE    "\033[1m\033[34m"      /* Bold Blue */
#define BOLDMAGENTA "\033[1m\033[35m"      /* Bold Magenta */
#define BOLDCYAN    "\033[1m\033[36m"      /* Bold Cyan */
#define BOLDWHITE   "\033[1m\033[37m"      /* Bold White */


#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include "HGCNtupleVariables.h"
#include "RecHit.h"
#include "LateralRings.h"
#include "RecHit_AHCAL.h"
#include "LateralSquares.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TGraph.h"
#include "TF1.h"
#include "TFile.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TDirectory.h"

class AnalyzeHGCOctTB : public HGCNtupleVariables{

 public:
  AnalyzeHGCOctTB(const TString &inputFileList="foo.txt", const char *outFileName="histo.root",const char *dataset="data", const char *config="alpha", const char* energy = "-1", bool shower = true, bool longi = true, bool trans = true);
  ~AnalyzeHGCOctTB();
  // Bool_t   FillChain(TChain *chain, TChain *chain2, TChain *chain3, const TString &inputFileList);
  Bool_t   FillChain(TChain *chain, const TString &inputFileList);
  Long64_t LoadTree(Long64_t entry);
  //void     EventLoop(const char *,const char *);
  void     EventLoop(const char *, bool, int);
  void     BookHistogram(const char *, const char *, const char* energy);

  void Module_Map_Init(const char *);
  void Alignment_Map_Init();
  void Noise_Map_Init();
  void Layer_Position_Init();
  void Weight_Map_Init();
  void Chi2_Weight_Map_Init();
  void Scaling_Factor_Init();

  bool SHOWER_RECO_HIST;
  bool LONGI_PROFILE_HIST;
  bool TRANSVERSE_PROFILE_HIST;
  TFile *oFile;

  
  ////////////// shower profile ///////////

  TDirectory *d_showerProfile;
  TDirectory *d_showerProfile_SummedEE;

  TDirectory *d_SS;
  TProfile *h_longi_profile_raw[40];
  TProfile *h_longi_profile_gev[40];

  TDirectory *d_SS_fraction;
  TProfile *h_longi_profile_raw_fraction[40];
  TProfile *h_longi_profile_gev_fraction[40];

  TProfile *h_longi_profile_MipsInEE_withAH;
  TProfile *h_longi_profile_MipsInEE_withAH_energy_fraction;
  TProfile *h_longi_profile_MipsInEE_withAH_abs_energy_fraction;


  TProfile *h_longi_profile_MipsInEE;
  TProfile *h_longi_profile_ShowerInEE;
  TProfile *h_longi_profile_inclusive;
  TProfile *h_longi_profile_inclusive_frac;

  TProfile *h_longi_profile_MipsInEE_gev;
  TProfile *h_longi_profile_ShowerInEE_gev;

  
  TProfile *h_longi_profile_MipsInEE_SS_ref;
  TH2F *h_longi_2D_MipsInEE_SS_ref;
  
  TProfile *h_longi_profile_MipsInEE_noWeight;
  TProfile *h_longi_profile_ShowerInEE_noWeight;
  TProfile *h_longi_profile_inclusive_noWeight;
  TProfile *h_longi_profile_MipsInEE_fractionalE;
  TProfile *h_longi_profile_ShowerInEE_fractionalE;
  TProfile *h_longi_profile_inclusive_fractionalE;
  TProfile *h_longi_profile_MipsInEE_abs_fractionalE;
  TProfile *h_longi_profile_ShowerInEE_abs_fractionalE;
  TProfile *h_longi_profile_inclusive_abs_fractionalE;
  TProfile *h_longi_profile_ShowerInEE_summed_up_EE;
  TProfile *h_longi_profile_ShowerInEE_summed_up_15EE;
  TProfile *h_longi_profile_ShowerInEE_summed_up_22EE;
  TProfile *h_longi_profile_ShowerInEE_summed_up_[40];

  
  TDirectory *d_showerProfile_layer;
  TProfile *h_longi_profile_raw_ShowerInEE_layer;
  TProfile *h_longi_profile_gev_ShowerInEE_layer;
  TProfile *h_longi_profile_raw_MipsInEE_layer;
  TProfile *h_longi_profile_gev_MipsInEE_layer;
  TProfile *h_longi_profile_raw_MipsInEE_SS_ref_layer;
  TProfile *h_longi_profile_raw_inclusive_layer;
  TProfile *h_longi_profile_raw_inclusive_frac_layer;
  TProfile *h_longi_profile_raw_ShowerInEE_noWeight_layer;
  TProfile *h_longi_profile_raw_MipsInEE_noWeight_layer;
  TProfile *h_longi_profile_raw_inclusive_noWeight_layer;
  TProfile *h_longi_profile_MipsInEE_fractionalE_layer;
  TProfile *h_longi_profile_ShowerInEE_fractionalE_layer;
  TProfile *h_longi_profile_inclusive_fractionalE_layer;
  TProfile *h_longi_profile_raw_SS10_check;
  TProfile *h_longi_profile_gev_SS10_check;

  TDirectory *d_SS_layer;
  TProfile *h_longi_profile_raw_layer[40];
  TProfile *h_longi_profile_gev_layer[40];

  TDirectory *d_SS_layer_fraction;
  TProfile *h_longi_profile_raw_layer_fraction[40];
  TProfile *h_longi_profile_gev_layer_fraction[40];

  
  /////////////////////////////////////////////////////////


  /////////////////////////////////////////////////////////////
  /////////// longitudinal, inclusive in SS for EE ////////////
  /////////////////////////////////////////////////////////////


  TDirectory *d_showerProfile_SS_inclusive;
  TDirectory *d_lambdaint;
  TProfile *h_longi_profile_raw_inc[16];
  TProfile *h_longi_profile_gev_inc[16];

  TDirectory *d_layer;
  TProfile *h_longi_profile_raw_layer_inc[16];
  TProfile *h_longi_profile_gev_layer_inc[16];
 


  /////////////////////////////////////////////////////////////

  /////////////////////////////////////////////////////  
  //////////////// transverse shower shape ////////////
  /////////////////////////////////////////////////////

  TDirectory *d_transverse;
  TDirectory *d_track_seed_diff;
  TH1F* h_track_seed_diff_x[40];
  TH1F* h_track_seed_diff_y[40];
  TH1F* h_track_seed_diff_dR[40];

  TDirectory *d_track_cog_diff;
  TH1F* h_track_cog_diff_x[40];
  TH1F* h_track_cog_diff_y[40];
  TH1F* h_track_cog_diff_dR[40];

  TDirectory *d_transverse_Edep;
  TDirectory *d_transverse_SS[40];
  TH1F* h_transverse_dR0p56[40][79];
  TH1F* h_transverse_dR1[40][79];
  TH1F* h_transverse_dR2[40][79];
  TH1F* h_transverse_dR3[40][79];
  TH1F* h_transverse_dR5[40][79];
  TH1F* h_transverse_dR8[40][79];
  TH1F* h_transverse_dR10[40][79];
  TH1F* h_transverse_dR12[40][79];
  TH1F* h_transverse_dR18[40][79];
  TH1F* h_transverse_dR20[40][79];

  TDirectory *d_transverse_prof;
  TDirectory *d_transverse_SS_prof[40];
  TProfile* h_transverse_prof_layer[40][79];
  // TProfile* h_transverse_dR5[40][79];
  // TProfile* h_transverse_dR8[40][79];
  // TProfile* h_transverse_dR10[40][79];
  TProfile* h_transverse_prof_EE[40];
  TProfile* h_transverse_prof_FH[40];
  TProfile* h_transverse_prof_AH[40];



  TDirectory *d_trans_distance_all;
  TDirectory *d_transverse_SS_dis_all[40];
  TH1F* h_transverse_distance_allSS[40][79];

  /////////////////////////////////////////////////////



  /////////////////////////////////////////////////////////////
  /////////// Transverse, inclusive in SS for EE ////////////
  /////////////////////////////////////////////////////////////


  TDirectory *d_transverse_inclusive;

  TDirectory *d_trans_mips;
  TDirectory *d_transverse_SS_inc_mips[16];
  TH1F* h_transverse_dR0p56_inc_mips[16][79];
  TH1F* h_transverse_dR1_inc_mips[16][79];
  TH1F* h_transverse_dR2_inc_mips[16][79];
  TH1F* h_transverse_dR3_inc_mips[16][79];
  TH1F* h_transverse_dR5_inc_mips[16][79];
  TH1F* h_transverse_dR8_inc_mips[16][79];
  TH1F* h_transverse_dR10_inc_mips[16][79];
  TH1F* h_transverse_dR12_inc_mips[16][79];
  TH1F* h_transverse_dR18_inc_mips[16][79];
  TH1F* h_transverse_dR20_inc_mips[16][79];


  TDirectory *d_trans_mips_summed;
  TDirectory *d_transverse_SS_inc_mips_summed[16];
  TH1F* h_transverse_dR0p56_inc_mips_summed[16][55];
  TH1F* h_transverse_dR1_inc_mips_summed[16][55];
  TH1F* h_transverse_dR2_inc_mips_summed[16][55];
  TH1F* h_transverse_dR3_inc_mips_summed[16][55];
  TH1F* h_transverse_dR5_inc_mips_summed[16][55];
  TH1F* h_transverse_dR8_inc_mips_summed[16][55];
  TH1F* h_transverse_dR10_inc_mips_summed[16][55];
  TH1F* h_transverse_dR12_inc_mips_summed[16][55];
  TH1F* h_transverse_dR18_inc_mips_summed[16][55];
  TH1F* h_transverse_dR20_inc_mips_summed[16][55];

    
  TDirectory *d_trans_gev;
  TDirectory *d_transverse_SS_inc_gev[16];
  TH1F* h_transverse_dR0p56_inc_gev[16][79];
  TH1F* h_transverse_dR1_inc_gev[16][79];
  TH1F* h_transverse_dR2_inc_gev[16][79];
  TH1F* h_transverse_dR3_inc_gev[16][79];
  TH1F* h_transverse_dR5_inc_gev[16][79];
  TH1F* h_transverse_dR8_inc_gev[16][79];
  TH1F* h_transverse_dR10_inc_gev[16][79];
  TH1F* h_transverse_dR12_inc_gev[16][79];
  TH1F* h_transverse_dR18_inc_gev[16][79];
  TH1F* h_transverse_dR20_inc_gev[16][79];

  TDirectory *d_trans_gev_summed;
  TDirectory *d_transverse_SS_inc_gev_summed[16];
  TH1F* h_transverse_dR0p56_inc_gev_summed[16][55];
  TH1F* h_transverse_dR1_inc_gev_summed[16][55];
  TH1F* h_transverse_dR2_inc_gev_summed[16][55];
  TH1F* h_transverse_dR3_inc_gev_summed[16][55];
  TH1F* h_transverse_dR5_inc_gev_summed[16][55];
  TH1F* h_transverse_dR8_inc_gev_summed[16][55];
  TH1F* h_transverse_dR10_inc_gev_summed[16][55];
  TH1F* h_transverse_dR12_inc_gev_summed[16][55];
  TH1F* h_transverse_dR18_inc_gev_summed[16][55];
  TH1F* h_transverse_dR20_inc_gev_summed[16][55];

  TDirectory *d_trans_gev_summed_ratio;
  TDirectory *d_transverse_SS_inc_gev_summed_ratio[16];
  TH1F* h_transverse_dR0p56_inc_gev_summed_ratio[16][55];
  TH1F* h_transverse_dR1_inc_gev_summed_ratio[16][55];
  TH1F* h_transverse_dR2_inc_gev_summed_ratio[16][55];
  TH1F* h_transverse_dR3_inc_gev_summed_ratio[16][55];
  TH1F* h_transverse_dR5_inc_gev_summed_ratio[16][55];
  TH1F* h_transverse_dR8_inc_gev_summed_ratio[16][55];
  TH1F* h_transverse_dR10_inc_gev_summed_ratio[16][55];
  TH1F* h_transverse_dR12_inc_gev_summed_ratio[16][55];
  TH1F* h_transverse_dR18_inc_gev_summed_ratio[16][55];
  TH1F* h_transverse_dR20_inc_gev_summed_ratio[16][55];


  
  TDirectory *d_trans_ratio;
  TDirectory *d_transverse_SS_inc_ratio[16];
  TH1F* h_transverse_dR0p56_dR5_inc_ratio[16][79];
  TH1F* h_transverse_dR1_dR5_inc_ratio[16][79];
  TH1F* h_transverse_dR2_dR5_inc_ratio[16][79];
  TH1F* h_transverse_dR3_dR5_inc_ratio[16][79];


  TDirectory *d_trans_ratio_summed;
  TDirectory *d_transverse_SS_inc_ratio_summed[16];
  TH1F* h_transverse_dR0p56_dR5_inc_ratio_summed[16][55];
  TH1F* h_transverse_dR1_dR5_inc_ratio_summed[16][55];
  TH1F* h_transverse_dR2_dR5_inc_ratio_summed[16][55];
  TH1F* h_transverse_dR3_dR5_inc_ratio_summed[16][55];

  TDirectory *d_trans_distance;
  TDirectory *d_transverse_SS_ratio[16];
  TDirectory *d_transverse_SS_dis[16];
  TH1F* h_transverse_distance[16][79];
  TH1F* h_transverse_distance_inclusive[16];
  TH1F* h_transverse_distance_inclusive_compart[16][3];

  TDirectory *d_trans_distance_zoom;
  TDirectory *d_transverse_SS_ratio_zoom[16];
  TDirectory *d_transverse_SS_dis_zoom[16];
  TH1F* h_transverse_distance_zoom[16][79];
  TH1F* h_transverse_distance_inclusive_zoom[16];
  TH1F* h_transverse_distance_inclusive_compart_zoom[16][3];

  TDirectory *d_trans_distance_prof;
  TDirectory *d_transverse_SS_ratio_prof[16];
  TDirectory *d_transverse_SS_dis_prof[16];
  TProfile* h_transverse_distance_prof[16][79];

  TProfile* h_transverse_distance_inclusive_prof[16];
  TProfile* h_transverse_distance_inclusive_compart_prof[16][3];

  // TProfile* p_htemp_check38;

  TDirectory *d_trans_prof_frac;
  TDirectory *d_trans_prof_frac_SS[3];
  TProfile* h_transverse_prof_frac[3][79];

  TProfile* h_transverse_prof_frac_Inc[3][3];

  
  /////////////////////////////////////////////////////////////


  // spandey debug 20 GeV /////

  TH1F* cog_x_db;
  TH1F* cog_y_db;
  TH2F* cog_xy_db;
  TH1F* run_db;
  TH1F* nrechit_db;
  TH1F* en_mip_db;
  TH1F* dR_db;
  TH1F* rechitx_db;
  TH1F* rechity_db;
  TH2F* rechitxy_db;

  // spandey debug, shape of trans dist
  TH1F* h_dR_SS1_7_L24_cog;
  TH1F* h_dR_SS1_7_L24_cog_weighted;
  TH1F* h_dR_SS1_7_L24_seed;
  TH1F* h_dR_SS1_7_L24_seed_weighted;
  TH1F* h_dR_SS1_7_L24_track;
  TH1F* h_dR_SS1_7_L24_track_weighted;

  TProfile* prof_energy_vs_dR_SS1_7_L24_cog;
  TProfile* prof_energy_vs_dR_SS1_7_L24_seed;
  TProfile* prof_energy_vs_dR_SS1_7_L24_track;

  



  /////////// D E B U G     :     dR(seed - track) > 2cm ////////
  // TH1F* h_debug_EE_layer;
  // TH1F* h_debug_FH_layer;

  // TH1F* h_debug_EE_energy;
  // TH1F* h_debug_FH_energy;

  
  //////////////////// weight scan /////////////////////

  // TDirectory *d_weightScan;
  // TDirectory *d_weightScan_EH;
  // TH1F* h_rechit_energy_FB_rel_weightScan[50];
  // TH1F* h_rechit_energy_EH_rel_weightScan[50];
  // TGraph* h_reso_weights;
  // TGraph* h_reso_weights_EH;


  /////////////////////////////////////////////////////


  //////////////////// track info /////////////////////

  TDirectory *d_trackInfo;
  TDirectory *d_trackLayer[40];


  
  TH1F* h_trackX[40];
  TH1F* h_trackY[40];
  TH2F* h_trackX_trackY[40];



  
  
  /////////////////////////////////////////////////////

  //////////////////// c o g ////// //////////////

  TDirectory *d_COG;
  TH1F* h_cogX[40];
  TH1F* h_cogY[40];
  TH2F* h_cogX_cogY[40];



  /////////////////////////////////////////////////////

  ///////////////// beam profile //////////////////////

  TDirectory *d_beamProfile;
  TDirectory *d_unWeighted;
  TDirectory *d_Weighted;


  TH2F* h_trackX_trackY_EE03;
  TH2F* h_trackX_trackY_EE03_corrected;
  TH2F* h_trackX_trackY_EE03_FULL_corrected;
  TH1F* h_trackX_EE03_track_corrected;
  TH1F* h_trackY_EE03_track_corrected;
  TH1F* h_rechitX_EE03;
  TH1F* h_rechitY_EE03;
  TH1F* h_trackX_EE03_FULL_corrected;
  TH1F* h_trackY_EE03_FULL_corrected;
  TH2F* h_rechitX_rechitY_EE03;
  TH2F* h_rechitX_rechitY_EE03_corrected;
  TH2F* h_trackX_rechitX_EE03;
  TH2F* h_trackX_rechitX_EE03_track_corrected;
  TH2F* h_trackX_rechitX_EE03_track_rechit_corrected;
  TH2F* h_trackY_rechitY_EE03;
  TH2F* h_trackY_rechitY_EE03_track_corrected;
  TH2F* h_trackY_rechitY_EE03_track_rechit_corrected;
  TH2F* h_dX_dY_EE03;
  TH2F* h_dX_dY_EE03_after_full_correction;

  TH2F* h_trackX_trackY_EE15;
  TH2F* h_trackX_trackY_EE15_corrected;
  TH2F* h_trackX_trackY_EE15_FULL_corrected;
  TH1F* h_trackX_EE15_track_corrected;
  TH1F* h_trackY_EE15_track_corrected;
  TH1F* h_rechitX_EE15;
  TH1F* h_rechitY_EE15;
  TH1F* h_trackX_EE15_FULL_corrected;
  TH1F* h_trackY_EE15_FULL_corrected;
  TH2F* h_rechitX_rechitY_EE15;
  TH2F* h_rechitX_rechitY_EE15_corrected;
  TH2F* h_trackX_rechitX_EE15;
  TH2F* h_trackX_rechitX_EE15_track_corrected;
  TH2F* h_trackX_rechitX_EE15_track_rechit_corrected;
  TH2F* h_trackY_rechitY_EE15;
  TH2F* h_trackY_rechitY_EE15_track_corrected;
  TH2F* h_trackY_rechitY_EE15_track_rechit_corrected;
  TH2F* h_dX_dY_EE15;
  TH2F* h_dX_dY_EE15_after_full_correction;

  TH2F* h_trackX_trackY_FH10;
  TH2F* h_trackX_trackY_FH10_corrected;
  TH2F* h_trackX_trackY_FH10_FULL_corrected;
  TH1F* h_trackX_FH10_track_corrected;
  TH1F* h_trackY_FH10_track_corrected;
  TH1F* h_rechitX_FH10;
  TH1F* h_rechitY_FH10;
  TH1F* h_trackX_FH10_FULL_corrected;
  TH1F* h_trackY_FH10_FULL_corrected;
  TH2F* h_rechitX_rechitY_FH10;
  TH2F* h_rechitX_rechitY_FH10_corrected;
  TH2F* h_trackX_rechitX_FH10;
  TH2F* h_trackX_rechitX_FH10_track_corrected;
  TH2F* h_trackX_rechitX_FH10_track_rechit_corrected;
  TH2F* h_trackY_rechitY_FH10;
  TH2F* h_trackY_rechitY_FH10_track_corrected;
  TH2F* h_trackY_rechitY_FH10_track_rechit_corrected;
  TH2F* h_dX_dY_FH10;
  TH2F* h_dX_dY_FH10_after_full_correction;


  /////////////////////////////////////////////////////

  ////////////// lateral distribution ///////////////

  TDirectory *d_lateral;
  TDirectory *d_debug;

  TH1F* h_E1_E7_EE;
  TH1F* h_E1_E19_EE;
  TH1F* h_E7_E19_EE;
  TH1F* h_E1_E7_FH;
  TH1F* h_E1_E19_FH;
  TH1F* h_E7_E19_FH;

  TH1F* h_E1_E7_HGCAL;
  TH1F* h_E1_E19_HGCAL;
  TH1F* h_E7_E19_HGCAL;

  TH1F* h_S1_S9_AH;
  TH1F* h_S1_S25_AH;
  TH1F* h_S9_S25_AH;

  TH1F* h_S1_S9_AH_debug;
  TH1F* h_S1_S25_AH_debug;
  TH1F* h_S9_S25_AH_debug;

  TH1F* h_S1_S9_AH_debug_inv;
  TH1F* h_S1_S25_AH_debug_inv;
  TH1F* h_S9_S25_AH_debug_inv;


  TH1F* h_E1_E7_SS_EE;
  TH1F* h_E1_E19_SS_EE;
  TH1F* h_E7_E19_SS_EE;

  TH1F* h_E1_E7_MIPs_in_EE;
  TH1F* h_E1_E19_MIPs_in_EE;
  TH1F* h_E7_E19_MIPs_in_EE;


  TH1F* h_E1_E7_SS_FH;
  TH1F* h_E1_E19_SS_FH;
  TH1F* h_E7_E19_SS_FH;

  TH1F* h_E1_E7_MIPs_in_FH;
  TH1F* h_E1_E19_MIPs_in_FH;
  TH1F* h_E7_E19_MIPs_in_FH;

  
  /////////////////////////////////////////////////////

  ////////////// energy & nrechit distribution /////////

  
  TDirectory *d_rechitEn_layer;
  TDirectory *d_rechitNrechit_layer;

  TDirectory *d_prof_layers[41];
  TDirectory *d_prof_layers_Nrec[41];

  TH1F *h_Rechits_En_SS[41][79];
  TH1F *h_Rechits_nrec_SS[41][79];
  
  TH1F *h_ShowerInEE_EE_fraction;
  TH1F *h_ShowerInEE_FH_fraction;
  TH1F *h_ShowerInEE_AH_fraction;
  TH1F *h_EE_fraction[41];
  TH1F *h_FH_fraction[41];
  TH1F *h_AH_fraction[41];
  /////////////////////////////////////////////////////

  /////////// shower energy reco //////////////////////

  TDirectory *d_shower_energy_reco;
  TDirectory *d_showering_in_EE;
  TDirectory *d_MIPs_in_EE;
  TDirectory *d_MIPs_in_FH;

  TH1F* h_EE_showerinEE;
  TH1F* h_FH_showerinEE;
  TH1F* h_AH_showerinEE;
  TH1F* h_all_showerinEE;
  TH1F* h_showerinEE_elecpi_scale;
  TH1F* h_showerinEE_elecpi_scale_withoutAH;

  TH1F* h_EE_MIPsinEE;
  TH1F* h_FH_MIPsinEE;
  TH1F* h_AH_MIPsinEE;
  TH1F* h_all_MIPsinEE;
  TH1F* h_MIPsinEE_elecpi_scale;
  TH1F* h_MIPsinEE_elecpi_scale_withoutAH;

  TH1F* h_EE_MIPsinFH;
  TH1F* h_FH_MIPsinFH;
  TH1F* h_AH_MIPsinFH;
  TH1F* h_all_MIPsinFH;
  TH1F* h_MIPsinFH_elecpi_scale;


  TH1F* h_EE_inclusive;
  TH1F* h_FH_inclusive;
  TH1F* h_AH_inclusive;
  TH1F* h_all_inclusive;
  TH1F* h_inclusive_elecpi_scale;

  
  /////////////////////////////////////////////////////

  //////////////////// FH debug ///////////////////////

  // TDirectory *d_FH_debug;
  // TH1F* h_nrechit_FH[12][7];
  // TH1F* h_energy_FH[12][7];


  /////////////////////////////////////////////////////
  
  ////////////// layerwise distribution ///////////////

  TDirectory *d_layerwise_distribution;
  TDirectory *d_ld_nrechit;
  TDirectory *d_ld_energy;
  
  TH1F* h_ld_nrechit[79];
  TH1F* h_ld_energy[79];

  /////////////////////////////////////////////////////

  /////////////////// selection cuts ////////////////

  TDirectory *d_selection_cut_check_showering_in_EE;
  TH1F* h_baseline;
  TH1F* h_SS1_reject;
  TH1F* h_SS2_reject;
  TH1F* h_trackwindow;
  TH1F* h_trackwindow_SS1_reject;
  TH1F* h_trackwindow_SS2_reject;


  TDirectory *d_selection_cut_check_MIPs_in_EE;
  TH1F* h_baseline_FH;
  TH1F* h_SS1_reject_FH;
  TH1F* h_SS2_reject_FH;
  TH1F* h_trackwindow_FH;
  TH1F* h_trackwindow_SS1_reject_FH;
  TH1F* h_trackwindow_SS2_reject_FH;

  /////////////////////////////////////////////////////

  /////////// selection cuts chi2 distributions ////////

  TDirectory *d_selection_cut_chi2;
  TDirectory *d_selcection_cut_EE;
  TDirectory *d_selcection_cut_FH;
  TH1F* h_baseline_chi2_EE;
  TH1F* h_SS1_reject_chi2_EE;
  TH1F* h_SS2_reject_chi2_EE;
  TH1F* h_trackwindow_chi2_EE;
  TH1F* h_trackwindow_SS1_reject_chi2_EE;
  TH1F* h_trackwindow_SS2_reject_chi2_EE;
  TH1F* h_trackwindow_SS2_reject_chi2_EE_recoE;
  
  TH1F* h_baseline_chi2_FH;
  TH1F* h_trackwindow_chi2_FH;
  TH1F* h_trackwindow_chi2_FH_recoE;
  TH1F* h_trackwindow_SS2_reject_chi2_FH;
  TH1F* h_trackwindow_SS2_reject_chi2_FH_recoE;

  TH1F* h_trackwindow_SS2_reject_chi2_all;
  TH1F* h_trackwindow_SS2_reject_chi2_all_recoE;

  
  /////////////////////////////////////////////////////

  
  
  TH2F* h_rechitX_rechitY_EE03_En_weighted;
  TH2F* h_rechitX_rechitY_EE15_En_weighted;


 
  TH1F* h_configuration;
  TH1F* h_particleID;
  TH1F* h_nTracks;
  TH1F* h_beamEnergy;
  TH1F* h_runNumber;
  TH1F* h_moduleID;
  TH1F* h_nRechits_layer[28];

  TH1F* h_rechit_energy_raw;
  TH1F* h_rechit_energy_raw_all;
  TH1F* h_rechit_energy_raw_EE;
  TH1F* h_rechit_energy_raw_FH;
  TH1F* h_rechit_energy_raw_AH;

  TH1F* h_rechit_energy_raw_AH_shower;
  TH1F* h_rechit_energy_raw_AH_shower_low_extended;



  TH1F* h_rechit_energy_raw_all_MipsInEEFH;
  TH1F* h_rechit_energy_raw_all_MipsInEE;
  TH1F* h_rechit_energy_raw_all_ShowerInEE;
  TH1F* h_rechit_energy_raw_all_ShowerInEE_modulo_EE1;


  TH1F* h_rechit_energy_raw_low_all;
  TH1F* h_rechit_energy_raw_low_EE;
  TH1F* h_rechit_energy_raw_low_FH;
  TH1F* h_rechit_energy_raw_low_AH;

  TH1F* h_rechit_energy_raw_verylow_all;
  TH1F* h_rechit_energy_raw_verylow_EE;
  TH1F* h_rechit_energy_raw_verylow_FH;
  TH1F* h_rechit_energy_raw_verylow_AH;



  TH1F* h_rechit_energy_raw_EE_MipsInEEFH;
  TH1F* h_rechit_energy_raw_FH_MipsInEEFH;
  TH1F* h_rechit_energy_raw_AH_MipsInEEFH;


  TH1F* h_cog_layer_EE;
  TH1F* h_cog_layer_FH;
  TH1F* h_cog_layer_AH;
  TH1F* h_cog_layer_HGCAL;



  TH2F* h_rechit_energy_raw_EE_vs_FH_MipsInEEFH;
  TH2F* h_rechit_energy_raw_FH_vs_AH_MipsInEEFH;
  TH2F* h_rechit_energy_raw_EE_vs_FH_MipsInEE;
  TH2F* h_rechit_energy_raw_FH_vs_AH_MipsInEE;
  TH2F* h_rechit_energy_raw_EE_vs_FH_ShowerInEE;
  TH2F* h_rechit_energy_raw_FH_vs_AH_ShowerInEE;


  TH2F* h_rechit_energy_raw_EE_vs_FH_all;
  TH1F* h_rechit_energy_raw_FH_MipsInEE_extended;
  TH1F* h_rechit_energy_raw_FH_MipsInEE_extended_v1;
  TH1F* h_rechit_energy_raw_AH_MipsInEE_extended;
  TH1F* h_rechit_energy_raw_AH_MipsInEE_extended_v1;

  TH1F* h_rechit_energy_raw_EE_MipsInFH_extended;
  TH1F* h_rechit_energy_raw_FH_MipsInFH_extended;
  TH1F* h_rechit_energy_raw_AH_MipsInFH_extended;
  TH1F* h_rechit_energy_raw_AH_MipsInFH_extended_v1;


  TH1F* h_rechit_energy_raw_EE_MipsInEEFH_extended;
  TH1F* h_rechit_energy_raw_FH_MipsInEEFH_extended;
  TH1F* h_rechit_energy_raw_AH_MipsInEEFH_extended;
  TH1F* h_rechit_energy_raw_AH_MipsInEEFH_extended_v1;

  TH1F* h_rechit_energy_raw_EE_MipsInEE;
  TH1F* h_rechit_energy_raw_EE_MipsInEE_extended;
  TH1F* h_rechit_energy_raw_FH_MipsInEE;
  TH1F* h_rechit_energy_raw_AH_MipsInEE;

  TH1F* h_rechit_energy_raw_EE_ShowerInEE;
  TH1F* h_rechit_energy_raw_FH_ShowerInEE;
  TH1F* h_rechit_energy_raw_AH_ShowerInEE;

  TH1F* h_SS_forLowerEwindow_ShowerInEE;

  TH1F* h_rechit_energy_raw_mid_all;
  TH1F* h_rechit_energy_raw_mid_EE;
  TH1F* h_rechit_energy_raw_mid_FH;

  TH1F* h_rechit_energy_raw_high_all;
  TH1F* h_rechit_energy_raw_high_EE;
  TH1F* h_rechit_energy_raw_high_FH;

 
  TH2F* h_trackX_trackY_EE1;
  TH2F* h_trackX_trackY_FH12;


  TH2F* h_Nrechit_EE_vs_FH;




  TH1F* h_Nrechit_EE;
  TH1F* h_Nrechit_FH;
  TH1F* h_Nrechit_AH;

  TH1F* h_Nrechit_low_EE;
  TH1F* h_Nrechit_low_FH;
  TH1F* h_Nrechit_low_AH;

  TH1F* h_shower_start;
  TH1F* h_shower_start_dN_dLambda;
  
  TH1F* h_shower_start_full_collapsed_EE;
  TH1F* h_shower_start_part_collapsed_EE;

  TH1F* h_shower_start_reg1;
  TH1F* h_shower_start_reg2;
  TH1F* h_shower_start_reg3;

  TH2F* h_Nrechit_EE_vs_FH_ShowerInEE; 
  TH2F* h_rechit_En_NRechits;


  int inEnergy_;
  const char *conf_;  

};
#endif

#ifdef AnalyzeHGCOctTB_cxx

void AnalyzeHGCOctTB::BookHistogram(const char *outFileName, const char* conf,  const char* energy) {
  
  char* hname = new char[1000];
  char* dir_name = new char[200];
  double y_max = -1.0;
  double x_2D_max = -1.0;
  double xbin = -1.0;
  double y_2D_max = -1.0;
  double mip_x_max = 3000;
  double mip_y_max = 1500;
  double r2_x_max = -1.0;
  double r2_x_bin = -1.0;
  
  double nrechit_xmax_EE = 0.0;
  double nrechit_xmax_FH = 0.0;
  double nrechit_xmax_AH = 0.0;


  if(!strcmp(energy, "20")) {
    y_max = 200.0;
    x_2D_max = 3000;
    xbin = 300;
    y_2D_max = 1500;
    inEnergy_ = 20;
    r2_x_max = inEnergy_*100;
    r2_x_bin = inEnergy_*10;
    nrechit_xmax_EE = 1000.0;
    nrechit_xmax_FH = 1000.0;
    nrechit_xmax_AH = 1000.0;
    
  }
  else if(!strcmp(energy, "50")) {
    y_max = 200.0;
    /* x_2D_max = 5000; */
    /* x_2D_max = 5000; */
    /* xbin = 500; */
    x_2D_max = 10000;
    xbin = 500;
    y_2D_max = 200;
    inEnergy_ = 50;
    r2_x_max = 0.5*inEnergy_*100;
    r2_x_bin = 0.5*inEnergy_*10;
    nrechit_xmax_EE = 1000.0;
    nrechit_xmax_FH = 1000.0;
    nrechit_xmax_AH = 1000.0;


  }
  else if(!strcmp(energy, "80")) {
    y_max = 200.0;
    /* x_2D_max = 8000; */
    /* xbin = 400; */
    x_2D_max = 12000;
    xbin = 300;
    y_2D_max = 3000;
    inEnergy_ = 80;
    r2_x_max = 0.5*inEnergy_*100;
    r2_x_bin = 0.5*inEnergy_*10;

    nrechit_xmax_EE = 1000.0;
    nrechit_xmax_FH = 1000.0;
    nrechit_xmax_AH = 1000.0;

  }
  else if(!strcmp(energy, "100")) {
    y_max = 300.0;
    /* x_2D_max = 10000; */
    /* xbin = 500; */
    x_2D_max = 18000;
    xbin = 360;
    y_2D_max = 3000;
    inEnergy_ = 100;
    r2_x_max = 0.5*inEnergy_*100;
    r2_x_bin = 0.5*inEnergy_*10;

    nrechit_xmax_EE = 1000.0;
    nrechit_xmax_FH = 1000.0;
    nrechit_xmax_AH = 1000.0;

  }
  else if(!strcmp(energy, "120")) {
    y_max = 400.0;
    /* x_2D_max = 12000; */
    /* xbin = 300; */
    x_2D_max = 20000;
    xbin = 400;
    y_2D_max = 3200;
    inEnergy_ = 120;
    r2_x_max = 0.5*inEnergy_*100;
    r2_x_bin = 0.5*inEnergy_*10;

    nrechit_xmax_EE = 1000.0;
    nrechit_xmax_FH = 1000.0;
    nrechit_xmax_AH = 2000.0;

  }
  else if(!strcmp(energy, "200")) {
    y_max = 600.0;
    /* x_2D_max = 18000; */
    /* xbin = 360; */
    x_2D_max = 40000;
    xbin = 800;
    y_2D_max = 4800;
    inEnergy_ = 200;
    r2_x_max = 0.5*inEnergy_*100;
    r2_x_bin = 0.5*inEnergy_*10;

    nrechit_xmax_EE = 2000.0;
    nrechit_xmax_FH = 1000.0;
    nrechit_xmax_AH = 2000.0;


  }
  else if(!strcmp(energy, "250")) {
    y_max = 700.0;
    /* x_2D_max = 25000; */
    /* xbin = 500; */
    x_2D_max = 50000;
    xbin = 1000;
    y_2D_max = 6000;
    inEnergy_ = 250;
    r2_x_max = 0.5*inEnergy_*100;
    r2_x_bin = 0.5*inEnergy_*10;

    nrechit_xmax_EE = 4000.0;
    nrechit_xmax_FH = 1000.0;
    nrechit_xmax_AH = 2000.0;


  }
  else if(!strcmp(energy, "300")) {
    y_max = 800.0;
    x_2D_max = 60000;
    /* x_2D_max = 20000; */
    xbin = 1000;
    y_2D_max = 6000;
    inEnergy_ = 300;
    r2_x_max = 0.5*inEnergy_*100;
    r2_x_bin = 0.5*inEnergy_*10;

    nrechit_xmax_EE = 4000.0;
    nrechit_xmax_FH = 1000.0;
    nrechit_xmax_AH = 2000.0;

  }
  else if(!strcmp(energy, "420")) {
    y_max = 20.0;
    x_2D_max = 200;
    /* x_2D_max = 20000; */
    xbin = 1000;
    y_2D_max = 6000;
    inEnergy_ = 420;
    r2_x_max = 0.5*200*100;
    r2_x_bin = 0.5*200*10;

    nrechit_xmax_EE = 1000.0;
    nrechit_xmax_FH = 1000.0;
    nrechit_xmax_AH = 1000.0;

  }

  else {
    cout<<"invalid energy!!!"<<endl; 
    return;
  }
  conf_ = conf;
  oFile = new TFile(outFileName, "recreate");


  if(TRANSVERSE_PROFILE_HIST) {
    cog_x_db = new TH1F("cog_x_db","cog_x_db",800,-40,40);
    cog_y_db = new TH1F("cog_y_db","cog_y_db",800,-40,40);
    cog_xy_db = new TH2F("cog_xy_db","cog_xy_db",800,-40,40,800,-40,40);
    run_db = new TH1F("run_db","run_db",100,650,700);
    nrechit_db = new TH1F("nrechit_db","nrechit_db",50,0,50);
    en_mip_db = new TH1F("en_mip_db","en_mip_db",800,0,400);
    dR_db = new TH1F("dR_db","dR_db",400,0,40);
    rechitx_db = new TH1F("rechitx_db","rechitx_db",800,-40,40);
    rechity_db = new TH1F("rechity_db","rechity_db",800,-40,40);
    rechitxy_db = new TH2F("rechitxy_db","rechitxy_db",800,-40,40,800,-40,40);

    h_dR_SS1_7_L24_cog = new TH1F("h_dR_SS1_7_L24_cog","h_dR_SS1_7_L24_cog",400,0,40);
    h_dR_SS1_7_L24_cog_weighted = new TH1F("h_dR_SS1_7_L24_cog_weighted","h_dR_SS1_7_L24_cog_weighted",400,0,40);
    h_dR_SS1_7_L24_seed = new TH1F("h_dR_SS1_7_L24_seed","h_dR_SS1_7_L24_seed",400,0,40);
    h_dR_SS1_7_L24_seed_weighted = new TH1F("h_dR_SS1_7_L24_seed_weighted","h_dR_SS1_7_L24_seed_weighted",400,0,40);
    h_dR_SS1_7_L24_track = new TH1F("h_dR_SS1_7_L24_track","h_dR_SS1_7_L24_track",400,0,40);
    h_dR_SS1_7_L24_track_weighted = new TH1F("h_dR_SS1_7_L24_track_weighted","h_dR_SS1_7_L24_track_weighted",400,0,40);

  
    prof_energy_vs_dR_SS1_7_L24_cog = new TProfile("prof_energy_vs_dR_SS1_7_L24_cog","prof_energy_vs_dR_SS1_7_L24_cog",400,0,40);
    prof_energy_vs_dR_SS1_7_L24_cog->GetXaxis()->SetTitle("dR_cog");
    prof_energy_vs_dR_SS1_7_L24_cog->GetYaxis()->SetTitle("<E>");
    prof_energy_vs_dR_SS1_7_L24_seed = new TProfile("prof_energy_vs_dR_SS1_7_L24_seed","prof_energy_vs_dR_SS1_7_L24_seed",400,0,40);
    prof_energy_vs_dR_SS1_7_L24_seed->GetXaxis()->SetTitle("dR_seed");
    prof_energy_vs_dR_SS1_7_L24_seed->GetYaxis()->SetTitle("<E>");
    prof_energy_vs_dR_SS1_7_L24_track = new TProfile("prof_energy_vs_dR_SS1_7_L24_track","prof_energy_vs_dR_SS1_7_L24_track",400,0,40);
    prof_energy_vs_dR_SS1_7_L24_track->GetXaxis()->SetTitle("dR_track");
    prof_energy_vs_dR_SS1_7_L24_track->GetYaxis()->SetTitle("<E>");

  }
  h_nTracks = new TH1F("h_nTracks","Tracks",12,-1,5);
  h_nTracks->GetXaxis()->SetTitle("Tracks");
  h_moduleID = new TH1F("h_moduleID","moduleID",4,-2,2);
  h_configuration = new TH1F("h_configuration","Configuration",60,-10,50);
  h_particleID = new TH1F("h_particleID","particleID",500,0,500);
  h_beamEnergy = new TH1F("h_beamEnergy","beamEnergy",500,0,500);
  h_runNumber = new TH1F("h_runNumber","RunNumber",2000,0,2000);
  h_trackX_trackY_EE1 = new TH2F("h_trackX_trackY_EE1","h_trackX_trackY_EE1",16,-8.0,8.0,16,-8.0,8.0);
  h_trackX_trackY_FH12 = new TH2F("h_trackX_trackY_FH12","h_trackX_trackY_FH12",16,-8.0,8.0,16,-8.0,8.0);

  


  h_cog_layer_EE = new TH1F("h_cog_layer_EE","h_cog_layer_EE",120,0,30);
  h_cog_layer_FH = new TH1F("h_cog_layer_FH","h_cog_layer_FH",120,0,15);
  h_cog_layer_HGCAL = new TH1F("h_cog_layer_HGCAL","h_cog_layer_HGCAL",168,0,42);
  h_cog_layer_AH = new TH1F("h_cog_layer_AH","h_cog_layer_AH",80,0,40);


  if(SHOWER_RECO_HIST) {
  
    h_Nrechit_EE = new TH1F("h_Nrechit_EE","NRechits in EE",nrechit_xmax_EE/2,0.0,nrechit_xmax_EE);
    h_Nrechit_FH = new TH1F("h_Nrechit_FH","NRechits in FH",nrechit_xmax_FH/2,0.0,nrechit_xmax_FH);
    h_Nrechit_AH = new TH1F("h_Nrechit_AH","NRechits in AH",nrechit_xmax_AH/2,0.0,nrechit_xmax_AH);


    h_Nrechit_low_EE = new TH1F("h_Nrechit_low_EE","NRechits in EE",100,0.0,100.0);
    h_Nrechit_low_FH = new TH1F("h_Nrechit_low_FH","NRechits in FH",100,0.0,100.0);
    h_Nrechit_low_AH = new TH1F("h_Nrechit_low_AH","NRechits in AH",100,0.0,100.0);



    h_rechit_energy_raw_FH_MipsInEE_extended = new TH1F("h_rechit_energy_raw_FH_MipsInEE_extended","h_rechit_energy_raw_FH_MipsInEE_extended",500.0,0,1000.0);
    h_rechit_energy_raw_FH_MipsInEE_extended_v1 = new TH1F("h_rechit_energy_raw_FH_MipsInEE_extended_v1","h_rechit_energy_raw_FH_MipsInEE_extended_v1",r2_x_bin,0.0,r2_x_max);
    h_rechit_energy_raw_AH_MipsInEE_extended = new TH1F("h_rechit_energy_raw_AH_MipsInEE_extended","h_rechit_energy_raw_AH_MipsInEE_extended",500.0,0,1000.0);
    h_rechit_energy_raw_AH_MipsInEE_extended_v1 = new TH1F("h_rechit_energy_raw_AH_MipsInEE_extended_v1","h_rechit_energy_raw_AH_MipsInEE_extended_v1",r2_x_bin,0.0,r2_x_max);


    h_rechit_energy_raw_EE_MipsInFH_extended = new TH1F("h_rechit_energy_raw_EE_MipsInFH_extended","h_rechit_energy_raw_EE_MipsInFH_extended",500.0,0,1000.0);
    h_rechit_energy_raw_FH_MipsInFH_extended = new TH1F("h_rechit_energy_raw_FH_MipsInFH_extended","h_rechit_energy_raw_FH_MipsInFH_extended",500.0,0,1000.0);
    h_rechit_energy_raw_AH_MipsInFH_extended = new TH1F("h_rechit_energy_raw_AH_MipsInFH_extended","h_rechit_energy_raw_AH_MipsInFH_extended",500.0,0,1000.0);
    h_rechit_energy_raw_AH_MipsInFH_extended_v1 = new TH1F("h_rechit_energy_raw_AH_MipsInFH_extended_v1","h_rechit_energy_raw_AH_MipsInFH_extended_v1",r2_x_bin,0.0,r2_x_max);



    h_rechit_energy_raw_EE_MipsInEEFH_extended = new TH1F("h_rechit_energy_raw_EE_MipsInEEFH_extended","h_rechit_energy_raw_EE_MipsInEEFH_extended",500.0,0,1000.0);
    h_rechit_energy_raw_EE_MipsInEEFH_extended->GetXaxis()->SetTitle("Energy in units of MIPs");
    h_rechit_energy_raw_FH_MipsInEEFH_extended = new TH1F("h_rechit_energy_raw_FH_MipsInEEFH_extended","h_rechit_energy_raw_FH_MipsInEEFH_extended",500.0,0,1000.0);
    h_rechit_energy_raw_FH_MipsInEEFH_extended->GetXaxis()->SetTitle("Energy in units of MIPs");
    h_rechit_energy_raw_AH_MipsInEEFH_extended = new TH1F("h_rechit_energy_raw_AH_MipsInEEFH_extended","h_rechit_energy_raw_AH_MipsInEEFH_extended",500.0,0,1000.0);
    h_rechit_energy_raw_AH_MipsInEEFH_extended->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_rechit_energy_raw_AH_MipsInEEFH_extended_v1 = new TH1F("h_rechit_energy_raw_AH_MipsInEEFH_extended_v1","h_rechit_energy_raw_AH_MipsInEEFH_extended_v1",r2_x_bin,0.0,r2_x_max);
    h_rechit_energy_raw_AH_MipsInEEFH_extended_v1->GetXaxis()->SetTitle("Energy in units of MIPs");





    h_rechit_energy_raw_EE_MipsInEEFH = new TH1F("h_rechit_energy_raw_EE_MipsInEEFH","h_rechit_energy_raw_EE_MipsInEEFH",250.0,0,500.0);
    h_rechit_energy_raw_EE_MipsInEEFH->GetXaxis()->SetTitle("Energy in units of MIPs");
    h_rechit_energy_raw_FH_MipsInEEFH = new TH1F("h_rechit_energy_raw_FH_MipsInEEFH","h_rechit_energy_raw_FH_MipsInEEFH",250.0,0,500.0);
    h_rechit_energy_raw_FH_MipsInEEFH->GetXaxis()->SetTitle("Energy in units of MIPs");
    h_rechit_energy_raw_AH_MipsInEEFH = new TH1F("h_rechit_energy_raw_AH_MipsInEEFH","h_rechit_energy_raw_AH_MipsInEEFH",250.0,0,500.0);
    h_rechit_energy_raw_AH_MipsInEEFH->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_rechit_energy_raw_EE_vs_FH_MipsInEEFH = new TH2F("h_rechit_energy_raw_EE_vs_FH_MipsInEEFH","h_rechit_energy_raw_EE_vs_FH_MipsInEEFH",250.0,0,500.0,250.0,0,500.0);
    h_rechit_energy_raw_EE_vs_FH_MipsInEEFH->GetXaxis()->SetTitle("Energy in EE units of MIPs");
    h_rechit_energy_raw_EE_vs_FH_MipsInEEFH->GetYaxis()->SetTitle("Energy in FH units of MIPs");

    h_rechit_energy_raw_FH_vs_AH_MipsInEEFH = new TH2F("h_rechit_energy_raw_FH_vs_AH_MipsInEEFH","h_rechit_energy_raw_FH_vs_AH_MipsInEEFH",250.0,0,500.0,250.0,0,500.0);
    h_rechit_energy_raw_FH_vs_AH_MipsInEEFH->GetXaxis()->SetTitle("Energy in FH units of MIPs");
    h_rechit_energy_raw_FH_vs_AH_MipsInEEFH->GetYaxis()->SetTitle("Energy in AH units of MIPs");
    h_rechit_energy_raw_FH_vs_AH_MipsInEEFH->GetYaxis()->SetTitleOffset(1.0);

    h_rechit_energy_raw_EE_MipsInEE = new TH1F("h_rechit_energy_raw_EE_MipsInEE","h_rechit_energy_raw_EE_MipsInEE",250.0,0,500.0);
    h_rechit_energy_raw_EE_MipsInEE_extended = new TH1F("h_rechit_energy_raw_EE_MipsInEE_extended","h_rechit_energy_raw_EE_MipsInEE_extended",500.0,0,1000.0);

    h_rechit_energy_raw_EE_MipsInEE->GetXaxis()->SetTitle("Energy in units of MIPs");
    h_rechit_energy_raw_FH_MipsInEE = new TH1F("h_rechit_energy_raw_FH_MipsInEE","h_rechit_energy_raw_FH_MipsInEE",r2_x_bin,0.0,r2_x_max);
    h_rechit_energy_raw_FH_MipsInEE->GetXaxis()->SetTitle("Energy in units of MIPs");
    h_rechit_energy_raw_AH_MipsInEE = new TH1F("h_rechit_energy_raw_AH_MipsInEE","h_rechit_energy_raw_AH_MipsInEE",r2_x_bin,0.0,r2_x_max);
    h_rechit_energy_raw_AH_MipsInEE->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_rechit_energy_raw_EE_vs_FH_MipsInEE = new TH2F("h_rechit_energy_raw_EE_vs_FH_MipsInEE","h_rechit_energy_raw_EE_vs_FH_MipsInEE",250.0,0,500.0,r2_x_bin,0.0,r2_x_max);
    h_rechit_energy_raw_EE_vs_FH_MipsInEE->GetXaxis()->SetTitle("Energy in EE units of MIPs");
    h_rechit_energy_raw_EE_vs_FH_MipsInEE->GetYaxis()->SetTitle("Energy in FH units of MIPs");
    h_rechit_energy_raw_EE_vs_FH_MipsInEE->GetYaxis()->SetTitleOffset(1.0);


    h_rechit_energy_raw_FH_vs_AH_MipsInEE = new TH2F("h_rechit_energy_raw_FH_vs_AH_MipsInEE","h_rechit_energy_raw_FH_vs_AH_MipsInEE",r2_x_bin,0.0,r2_x_max,r2_x_bin,0.0,r2_x_max);
    h_rechit_energy_raw_FH_vs_AH_MipsInEE->GetXaxis()->SetTitle("Energy in FH units of MIPs");
    h_rechit_energy_raw_FH_vs_AH_MipsInEE->GetYaxis()->SetTitle("Energy in AH units of MIPs");
    h_rechit_energy_raw_FH_vs_AH_MipsInEE->GetYaxis()->SetTitleOffset(1.0);



    h_rechit_energy_raw_EE_ShowerInEE = new TH1F("h_rechit_energy_raw_EE_ShowerInEE","h_rechit_energy_raw_EE_ShowerInEE",xbin,0.0,x_2D_max);
    h_rechit_energy_raw_EE_ShowerInEE->GetXaxis()->SetTitle("Energy in units of MIPs");
    h_rechit_energy_raw_FH_ShowerInEE = new TH1F("h_rechit_energy_raw_FH_ShowerInEE","h_rechit_energy_raw_FH_ShowerInEE",r2_x_bin,0.0,r2_x_max);
    h_rechit_energy_raw_FH_ShowerInEE->GetXaxis()->SetTitle("Energy in units of MIPs");
    h_rechit_energy_raw_AH_ShowerInEE = new TH1F("h_rechit_energy_raw_AH_ShowerInEE","h_rechit_energy_raw_AH_ShowerInEE",r2_x_bin,0.0,r2_x_max);
    h_rechit_energy_raw_AH_ShowerInEE->GetXaxis()->SetTitle("Energy in units of MIPs");


    h_SS_forLowerEwindow_ShowerInEE = new TH1F("h_SS_forLowerEwindow_ShowerInEE","SS_forLowerEwindow_ShowerInEE, 100<EE<200 mips",88,-2,42);


    h_rechit_energy_raw_EE_vs_FH_ShowerInEE = new TH2F("h_rechit_energy_raw_EE_vs_FH_ShowerInEE","h_rechit_energy_raw_EE_vs_FH_ShowerInEE",xbin,0.0,x_2D_max,r2_x_bin,0.0,r2_x_max);
    h_rechit_energy_raw_EE_vs_FH_ShowerInEE->GetXaxis()->SetTitle("Energy in EE units of MIPs");
    h_rechit_energy_raw_EE_vs_FH_ShowerInEE->GetYaxis()->SetTitle("Energy in FH units of MIPs");
    h_rechit_energy_raw_EE_vs_FH_ShowerInEE->GetYaxis()->SetTitleOffset(1.0);

    h_rechit_energy_raw_FH_vs_AH_ShowerInEE = new TH2F("h_rechit_energy_raw_FH_vs_AH_ShowerInEE","h_rechit_energy_raw_FH_vs_AH_ShowerInEE",xbin,0.0,x_2D_max,r2_x_bin,0.0,r2_x_max);
    h_rechit_energy_raw_FH_vs_AH_ShowerInEE->GetXaxis()->SetTitle("Energy in FH units of MIPs");
    h_rechit_energy_raw_FH_vs_AH_ShowerInEE->GetYaxis()->SetTitle("Energy in AH units of MIPs");
    h_rechit_energy_raw_FH_vs_AH_ShowerInEE->GetYaxis()->SetTitleOffset(1.0);


    h_rechit_energy_raw_EE_vs_FH_all = new TH2F("h_rechit_energy_raw_EE_vs_FH_all","h_rechit_energy_raw_EE_vs_FH_all",xbin,0.0,x_2D_max,xbin,0.0,x_2D_max);
    h_rechit_energy_raw_EE_vs_FH_all->GetXaxis()->SetTitle("Energy in EE units of MIPs");
    h_rechit_energy_raw_EE_vs_FH_all->GetYaxis()->SetTitle("Energy in FH units of MIPs");



    h_rechit_energy_raw = new TH1F("h_rechit_energy_raw","FH+AHCAL mips",xbin,0.0,x_2D_max);
    h_rechit_energy_raw->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_rechit_energy_raw_EE = new TH1F("h_rechit_energy_raw_EE","h_rechit_energy_raw_EE",xbin,0.0,x_2D_max);
    h_rechit_energy_raw_EE->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_rechit_energy_raw_FH = new TH1F("h_rechit_energy_raw_FH","h_rechit_energy_raw_FH",xbin,0.0,x_2D_max);
    h_rechit_energy_raw_FH->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_rechit_energy_raw_AH = new TH1F("h_rechit_energy_raw_AH","h_rechit_energy_raw_AH",xbin,0.0,x_2D_max);
    h_rechit_energy_raw_AH->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_rechit_energy_raw_AH_shower = new TH1F("h_rechit_energy_raw_AH_shower","h_rechit_energy_raw_AH_shower",xbin,0.0,x_2D_max/2);
    h_rechit_energy_raw_AH_shower->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_rechit_energy_raw_AH_shower_low_extended = new TH1F("h_rechit_energy_raw_AH_shower_low_extended","h_rechit_energy_raw_AH_shower_low_extended",500,0.0,1000.0);
    h_rechit_energy_raw_AH_shower_low_extended->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_rechit_energy_raw_all = new TH1F("h_rechit_energy_raw_all","h_rechit_energy_raw_all",xbin,0.0,x_2D_max);
    h_rechit_energy_raw_all->GetXaxis()->SetTitle("Energy in units of MIPs");


    h_rechit_energy_raw_all_MipsInEEFH = new TH1F("h_rechit_energy_raw_all_MipsInEEFH","h_rechit_energy_raw_all_MipsInEEFH",xbin,0.0,x_2D_max);
    h_rechit_energy_raw_all_MipsInEEFH->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_rechit_energy_raw_all_MipsInEE = new TH1F("h_rechit_energy_raw_all_MipsInEE","h_rechit_energy_raw_all_MipsInEE",xbin,0.0,x_2D_max);
    h_rechit_energy_raw_all_MipsInEE->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_rechit_energy_raw_all_ShowerInEE = new TH1F("h_rechit_energy_raw_all_ShowerInEE","h_rechit_energy_raw_all_ShowerInEE",xbin,0.0,x_2D_max);
    h_rechit_energy_raw_all_ShowerInEE->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_rechit_energy_raw_all_ShowerInEE_modulo_EE1 = new TH1F("h_rechit_energy_raw_all_ShowerInEE_modulo_EE1","h_rechit_energy_raw_all_ShowerInEE_modulo_EE1",xbin,0.0,x_2D_max);
    h_rechit_energy_raw_all_ShowerInEE_modulo_EE1->GetXaxis()->SetTitle("Energy in units of MIPs");


    h_rechit_energy_raw_low_EE = new TH1F("h_rechit_energy_raw_low_EE","h_rechit_energy_raw_low_EE",250.0,0.0,250.0);
    h_rechit_energy_raw_low_EE->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_rechit_energy_raw_low_FH = new TH1F("h_rechit_energy_raw_low_FH","h_rechit_energy_raw_low_FH",300,0.0,300.0);
    h_rechit_energy_raw_low_FH->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_rechit_energy_raw_low_AH = new TH1F("h_rechit_energy_raw_low_AH","h_rechit_energy_raw_low_AH",300,0.0,300.0);
    h_rechit_energy_raw_low_AH->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_rechit_energy_raw_low_all = new TH1F("h_rechit_energy_raw_low_all","h_rechit_energy_raw_low_all",500,0.0,500.0);
    h_rechit_energy_raw_low_all->GetXaxis()->SetTitle("Energy in units of MIPs");


    h_rechit_energy_raw_verylow_EE = new TH1F("h_rechit_energy_raw_verylow_EE","h_rechit_energy_raw_verylow_EE",500.0,0.0,50.0);
    h_rechit_energy_raw_verylow_EE->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_rechit_energy_raw_verylow_FH = new TH1F("h_rechit_energy_raw_verylow_FH","h_rechit_energy_raw_verylow_FH",500,0.0,50.0);
    h_rechit_energy_raw_verylow_FH->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_rechit_energy_raw_verylow_AH = new TH1F("h_rechit_energy_raw_verylow_AH","h_rechit_energy_raw_verylow_AH",500,0.0,50.0);
    h_rechit_energy_raw_verylow_AH->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_rechit_energy_raw_verylow_all = new TH1F("h_rechit_energy_raw_verylow_all","h_rechit_energy_raw_verylow_all",500,0.0,50.0);
    h_rechit_energy_raw_verylow_all->GetXaxis()->SetTitle("Energy in units of MIPs");


    float x_mid_range_min = x_2D_max*0.1;
    float x_mid_range_max = x_2D_max*0.66;
    float x_mid_bin = (x_mid_range_max - x_mid_range_min)/2;

    h_rechit_energy_raw_mid_EE = new TH1F("h_rechit_energy_raw_mid_EE","h_rechit_energy_raw_mid_EE",x_mid_bin,x_mid_range_min,x_mid_range_max);
    h_rechit_energy_raw_mid_EE->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_rechit_energy_raw_mid_FH = new TH1F("h_rechit_energy_raw_mid_FH","h_rechit_energy_raw_mid_FH",x_mid_bin,x_mid_range_min,x_mid_range_max);
    h_rechit_energy_raw_mid_FH->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_rechit_energy_raw_mid_all = new TH1F("h_rechit_energy_raw_mid_all","h_rechit_energy_raw_mid_all",x_mid_bin,x_mid_range_min,x_mid_range_max);
    h_rechit_energy_raw_mid_all->GetXaxis()->SetTitle("Energy in units of MIPs");

    float x_high_range_min = x_2D_max*0.25;
    float x_high_range_max = x_2D_max*0.7;
    float x_high_bin = (x_high_range_max - x_high_range_min)/4;

    h_rechit_energy_raw_high_EE = new TH1F("h_rechit_energy_raw_high_EE","h_rechit_energy_raw_high_EE",x_high_bin,x_high_range_min,x_high_range_max);
    h_rechit_energy_raw_high_EE->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_rechit_energy_raw_high_FH = new TH1F("h_rechit_energy_raw_high_FH","h_rechit_energy_raw_high_FH",x_high_bin,x_high_range_min,x_high_range_max);
    h_rechit_energy_raw_high_FH->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_rechit_energy_raw_high_all = new TH1F("h_rechit_energy_raw_high_all","h_rechit_energy_raw_high_all",x_high_bin,x_high_range_min,x_high_range_max);
    h_rechit_energy_raw_high_all->GetXaxis()->SetTitle("Energy in units of MIPs");



    h_Nrechit_EE_vs_FH = new TH2F("h_Nrechit_EE_vs_FH","NRechits EE vs FH",500,0.0,1000.0,500,0.0,1000.0);
    h_Nrechit_EE_vs_FH->GetXaxis()->SetTitle("NRechits in EE");
    h_Nrechit_EE_vs_FH->GetYaxis()->SetTitle("NRechits in FH");
    h_Nrechit_EE_vs_FH->GetYaxis()->SetTitleOffset(1.0);


  }


  
  h_shower_start = new TH1F("h_shower_start","h_shower_start",10000,-2.0,5.2);
  h_shower_start->GetXaxis()->SetTitle("#lambda_{0}");
  h_shower_start->Sumw2();

  h_shower_start_dN_dLambda = new TH1F("h_shower_start_dN_dLambda","h_shower_start_dN_dLambda",10000,-2.0,5.2);
  h_shower_start_dN_dLambda->GetXaxis()->SetTitle("#lambda_{0}");
  h_shower_start_dN_dLambda->Sumw2();

  h_shower_start_full_collapsed_EE = new TH1F("h_shower_start_full_collapsed_EE","h_shower_start_full_collapsed_EE",10000,-2.0,5.2);
  h_shower_start_full_collapsed_EE->GetXaxis()->SetTitle("#lambda_{0}");
  h_shower_start_full_collapsed_EE->Sumw2();

  h_shower_start_part_collapsed_EE = new TH1F("h_shower_start_part_collapsed_EE","h_shower_start_part_collapsed_EE",10000,-2.0,5.2);
  h_shower_start_part_collapsed_EE->GetXaxis()->SetTitle("#lambda_{0}");
  h_shower_start_part_collapsed_EE->Sumw2();

  
  h_shower_start_reg1 = new TH1F("h_shower_start_reg1","h_shower_start_reg1",10000,-2.0,5.2);
  h_shower_start_reg1->GetXaxis()->SetTitle("#lambda_{0}");
  h_shower_start_reg1->Sumw2();

  
  h_shower_start_reg2 = new TH1F("h_shower_start_reg2","h_shower_start_reg2",10000,-2.0,5.2);
  h_shower_start_reg2->GetXaxis()->SetTitle("#lambda_{0}");
  h_shower_start_reg2->Sumw2();

  
  h_shower_start_reg3 = new TH1F("h_shower_start_reg3","h_shower_start_reg3",10000,-2.0,5.2);
  h_shower_start_reg3->GetXaxis()->SetTitle("#lambda_{0}");
  h_shower_start_reg3->Sumw2();



 
  h_Nrechit_EE_vs_FH_ShowerInEE = new TH2F("h_Nrechit_EE_vs_FH_ShowerInEE","NRechits EE vs FH",500,0.0,1000.0,500,0.0,1000.0);
  h_Nrechit_EE_vs_FH_ShowerInEE->GetXaxis()->SetTitle("NRechits in EE");
  h_Nrechit_EE_vs_FH_ShowerInEE->GetYaxis()->SetTitle("NRechits in FH");
  h_Nrechit_EE_vs_FH_ShowerInEE->GetYaxis()->SetTitleOffset(1.0);

  char* h_name = new char[500];  

  if(LONGI_PROFILE_HIST) {


    d_showerProfile = oFile->mkdir("shower_profile_lambdaPi");
    d_showerProfile->cd();



    sprintf(h_name,"h_longi_profile_raw_EE_ShowerInEE");
    h_longi_profile_ShowerInEE = new TProfile(h_name,h_name,10000,-2.0,10.0);
    sprintf(h_name,"h_longi_profile_gev_EE_ShowerInEE");
    h_longi_profile_ShowerInEE_gev = new TProfile(h_name,h_name,10000,-2.0,10.0);
    sprintf(h_name,"h_longi_profile_raw_EE_MipsInEE");
    h_longi_profile_MipsInEE = new TProfile(h_name,h_name,10000,-2.0,10.0);
    sprintf(h_name,"h_longi_profile_gev_EE_MipsInEE");
    h_longi_profile_MipsInEE_gev = new TProfile(h_name,h_name,10000,-2.0,10.0);
    sprintf(h_name,"h_longi_profile_raw_EE_MipsInEE_SS_ref");
    h_longi_profile_MipsInEE_SS_ref = new TProfile(h_name,h_name,10000,-2.0,10.0);
    sprintf(h_name,"h_longi_2D_raw_EE_MipsInEE_SS_ref");
    h_longi_2D_MipsInEE_SS_ref = new TH2F(h_name,h_name,10000,-2.0,10.0, 1000, 0.0, 1000);
    sprintf(h_name,"h_longi_profile_raw_inclusive");
    h_longi_profile_inclusive = new TProfile(h_name,h_name,10000,-2.0,10.0);
    sprintf(h_name,"h_longi_profile_raw_inclusive_frac");
    h_longi_profile_inclusive_frac = new TProfile(h_name,h_name,10000,-2.0,10.0);



  
    d_SS = d_showerProfile->mkdir("profile_SS");
    d_SS->cd();
    for(int i = 0; i < 40; i++) {
      sprintf(h_name,"h_longi_profile_raw_SS_%02d",i+1);
      h_longi_profile_raw[i] = new TProfile(h_name,h_name,10000, -2.0, 10);
      h_longi_profile_raw[i]->GetXaxis()->SetTitle("#lambda_{int}");
      h_longi_profile_raw[i]->GetYaxis()->SetTitle("Mean energy deposited (MIPs)");

      sprintf(h_name,"h_longi_profile_raw_gev_%02d",i+1);
      h_longi_profile_gev[i] = new TProfile(h_name,h_name,10000, -2.0, 10);
      h_longi_profile_gev[i]->GetXaxis()->SetTitle("#lambda_{int}");
      h_longi_profile_gev[i]->GetYaxis()->SetTitle("Mean energy deposited (GeV)");


    }

    d_SS_fraction = d_showerProfile->mkdir("profile_SS_fraction");
    d_SS_fraction->cd();
    for(int i = 0; i < 40; i++) {
      sprintf(h_name,"h_longi_profile_raw_SS_%02d_fraction",i+1);
      h_longi_profile_raw_fraction[i] = new TProfile(h_name,h_name,10000, -2.0, 10);
      h_longi_profile_raw_fraction[i]->GetXaxis()->SetTitle("#lambda_{int}");
      h_longi_profile_raw_fraction[i]->GetYaxis()->SetTitle("Mean energy deposited (MIPs)");

      sprintf(h_name,"h_longi_profile_raw_gev_%02d_fraction",i+1);
      h_longi_profile_gev_fraction[i] = new TProfile(h_name,h_name,10000, -2.0, 10);
      h_longi_profile_gev_fraction[i]->GetXaxis()->SetTitle("#lambda_{int}");
      h_longi_profile_gev_fraction[i]->GetYaxis()->SetTitle("Mean energy deposited (GeV)");


    }

    d_showerProfile_SummedEE = d_showerProfile->mkdir("Summed_EE");
    d_showerProfile_SummedEE->cd();
  
    for(int i = 0; i < 40; i++) {
      sprintf(h_name,"h_longi_profile_ShowerInEE_summed_up_%02dEE",i+1);
    
      h_longi_profile_ShowerInEE_summed_up_[i] = new TProfile(h_name,h_name,10000,-2.0,10.0);

    }


    d_showerProfile_layer = oFile->mkdir("shower_profile_layer");
    d_showerProfile_layer->cd();
    sprintf(h_name,"h_longi_profile_raw_EE_ShowerInEE_layer");
    h_longi_profile_raw_ShowerInEE_layer = new TProfile(h_name,h_name,810,-1,80);
    sprintf(h_name,"h_longi_profile_gev_EE_ShowerInEE_layer");
    h_longi_profile_gev_ShowerInEE_layer = new TProfile(h_name,h_name,810,-1,80);
    sprintf(h_name,"h_longi_profile_raw_EE_MipsInEE_layer");
    h_longi_profile_raw_MipsInEE_layer = new TProfile(h_name,h_name,810,-1,80);
    sprintf(h_name,"h_longi_profile_gev_EE_MipsInEE_layer");
    h_longi_profile_gev_MipsInEE_layer = new TProfile(h_name,h_name,810,-1,80);
    sprintf(h_name,"h_longi_profile_raw_EE_MipsInEE_SS_ref_layer");
    h_longi_profile_raw_MipsInEE_SS_ref_layer = new TProfile(h_name,h_name,810,-1,80);
    sprintf(h_name,"h_longi_profile_raw_inclusive_layer");
    h_longi_profile_raw_inclusive_layer = new TProfile(h_name,h_name,810,-1,80);
    sprintf(h_name,"h_longi_profile_raw_inclusive_frac_layer");
    h_longi_profile_raw_inclusive_frac_layer = new TProfile(h_name,h_name,810,-1,80);
    sprintf(h_name,"h_longi_profile_raw_EE_ShowerInEE_noWeight_layer");
    h_longi_profile_raw_ShowerInEE_noWeight_layer = new TProfile(h_name,h_name,810,-1,80);
    sprintf(h_name,"h_longi_profile_raw_EE_MipsInEE_noWeight_layer");
    h_longi_profile_raw_MipsInEE_noWeight_layer = new TProfile(h_name,h_name,810,-1,80);
    sprintf(h_name,"h_longi_profile_raw_inclusive_noWeight_layer");
    h_longi_profile_raw_inclusive_noWeight_layer = new TProfile(h_name,h_name,810,-1,80);
  
    sprintf(h_name,"h_longi_profile_MipsInEE_fractionalE_layer");
    h_longi_profile_MipsInEE_fractionalE_layer = new TProfile(h_name,h_name,810,-1,80);
    sprintf(h_name,"h_longi_profile_ShowerInEE_fractionalE_layer");
    h_longi_profile_ShowerInEE_fractionalE_layer = new TProfile(h_name,h_name,810,-1,80);
    sprintf(h_name,"h_longi_profile_inclusive_fractionalE_layer");
    h_longi_profile_inclusive_fractionalE_layer = new TProfile(h_name,h_name,810,-1,80);

  
    h_longi_profile_raw_SS10_check = new TProfile("h_longi_profile_raw_SS10_check","h_longi_profile_raw_SS10_check",810,-1,80);
    h_longi_profile_raw_SS10_check->GetXaxis()->SetTitle("Layer #");
    h_longi_profile_raw_SS10_check->GetYaxis()->SetTitle("Mean energy deposited (MIPs)");

    h_longi_profile_gev_SS10_check = new TProfile("h_longi_profile_gev_SS10_check","h_longi_profile_gev_SS10_check",810,-1,80);
    h_longi_profile_gev_SS10_check->GetXaxis()->SetTitle("Layer #");
    h_longi_profile_gev_SS10_check->GetYaxis()->SetTitle("Mean energy deposited (MIPs)");


    d_SS_layer = d_showerProfile_layer->mkdir("profile_SS_layer");
    d_SS_layer->cd();
    for(int i = 0; i < 40; i++) {
      sprintf(h_name,"h_longi_profile_raw_SS_%02d_layer",i+1);
      h_longi_profile_raw_layer[i] = new TProfile(h_name,h_name,810,-1,80);
      h_longi_profile_raw_layer[i]->GetXaxis()->SetTitle("Layer #");
      h_longi_profile_raw_layer[i]->GetYaxis()->SetTitle("Mean energy deposited (MIPs)");
      sprintf(h_name,"h_longi_profile_gev_SS_%02d_layer",i+1);
      h_longi_profile_gev_layer[i] = new TProfile(h_name,h_name,810,-1,80);
      h_longi_profile_gev_layer[i]->GetXaxis()->SetTitle("Layer #");
      h_longi_profile_gev_layer[i]->GetYaxis()->SetTitle("Mean energy deposited (GeV)");

    
    }

    d_SS_layer_fraction = d_showerProfile_layer->mkdir("profile_SS_layer_fraction");
    d_SS_layer_fraction->cd();
    for(int i = 0; i < 40; i++) {
      sprintf(h_name,"h_longi_profile_raw_SS_%02d_layer_fraction",i+1);
      h_longi_profile_raw_layer_fraction[i] = new TProfile(h_name,h_name,810,-1,80);
      h_longi_profile_raw_layer_fraction[i]->GetXaxis()->SetTitle("Layer #");
      h_longi_profile_raw_layer_fraction[i]->GetYaxis()->SetTitle("Mean energy deposited (MIPs)");
      sprintf(h_name,"h_longi_profile_gev_SS_%02d_layer_fraction",i+1);
      h_longi_profile_gev_layer_fraction[i] = new TProfile(h_name,h_name,810,-1,80);
      h_longi_profile_gev_layer_fraction[i]->GetXaxis()->SetTitle("Layer #");
      h_longi_profile_gev_layer_fraction[i]->GetYaxis()->SetTitle("Mean energy deposited (GeV)");

    
    }

    

    ////////////////////////////////////
    // longi shower shape, inclusive  //
    ///////////////////////////////////
    
    d_showerProfile_SS_inclusive = oFile->mkdir("shower_profile_SS_inclusive");
    d_showerProfile_SS_inclusive->cd();
    d_lambdaint = d_showerProfile_SS_inclusive->mkdir("lambda_int");
    d_lambdaint->cd();

    for(int i = 0; i < 16; i++) {
      switch(i) {
      case 0: sprintf(h_name,"h_longi_profile_raw_SS_01_07_inc"); break;
      case 1: sprintf(h_name,"h_longi_profile_raw_SS_08_14_inc"); break;
      case 2: sprintf(h_name,"h_longi_profile_raw_SS_15_21_inc"); break;
      case 3: sprintf(h_name,"h_longi_profile_raw_SS_22_28_inc"); break;
      default: sprintf(h_name,"h_longi_profile_raw_SS_%02d_inc",i+1+24); break;
      }
      h_longi_profile_raw_inc[i] = new TProfile(h_name,h_name,10000, -2.0, 10);
      h_longi_profile_raw_inc[i]->GetXaxis()->SetTitle("#lambda_{int}");
      h_longi_profile_raw_inc[i]->GetYaxis()->SetTitle("Mean energy deposited (MIPs)");

      switch(i) {
      case 0: sprintf(h_name,"h_longi_profile_gev_SS_01_07_inc"); break;
      case 1: sprintf(h_name,"h_longi_profile_gev_SS_08_14_inc"); break;
      case 2: sprintf(h_name,"h_longi_profile_gev_SS_15_21_inc"); break;
      case 3: sprintf(h_name,"h_longi_profile_gev_SS_22_28_inc"); break;
      default: sprintf(h_name,"h_longi_profile_gev_SS_%02d_inc",i+1+24); break;
      }
      h_longi_profile_gev_inc[i] = new TProfile(h_name,h_name,10000, -2.0, 10);
      h_longi_profile_gev_inc[i]->GetXaxis()->SetTitle("#lambda_{int}");
      h_longi_profile_gev_inc[i]->GetYaxis()->SetTitle("Mean energy deposited (GeV)");
   

    }

    d_layer = d_showerProfile_SS_inclusive->mkdir("layer");
    d_layer->cd();
    for(int i = 0; i < 16; i++) {
      switch(i) {
      case 0: sprintf(h_name,"h_longi_profile_raw_SS_01_07_layer_inc"); break;
      case 1: sprintf(h_name,"h_longi_profile_raw_SS_08_14_layer_inc"); break;
      case 2: sprintf(h_name,"h_longi_profile_raw_SS_15_21_layer_inc"); break;
      case 3: sprintf(h_name,"h_longi_profile_raw_SS_22_28_layer_inc"); break;
      default: sprintf(h_name,"h_longi_profile_raw_SS_%02d_layer_inc",i+1+24); break;
      }
      h_longi_profile_raw_layer_inc[i] = new TProfile(h_name,h_name,810,-1,80);
      h_longi_profile_raw_layer_inc[i]->GetXaxis()->SetTitle("Layer");
      h_longi_profile_raw_layer_inc[i]->GetYaxis()->SetTitle("Mean energy deposited (MIPs)");

      switch(i) {
      case 0: sprintf(h_name,"h_longi_profile_gev_SS_01_07_layer_inc"); break;
      case 1: sprintf(h_name,"h_longi_profile_gev_SS_08_14_layer_inc"); break;
      case 2: sprintf(h_name,"h_longi_profile_gev_SS_15_21_layer_inc"); break;
      case 3: sprintf(h_name,"h_longi_profile_gev_SS_22_28_layer_inc"); break;
      default: sprintf(h_name,"h_longi_profile_gev_SS_%02d_layer_inc",i+1+24); break;
      }
      h_longi_profile_gev_layer_inc[i] = new TProfile(h_name,h_name,810,-1,80);
      h_longi_profile_gev_layer_inc[i]->GetXaxis()->SetTitle("Layer");
      h_longi_profile_gev_layer_inc[i]->GetYaxis()->SetTitle("Mean energy deposited (GeV)");
   

    }

    ////////////////////////////////////////////////////////

  }

  char* SS_text = new char[500];

  /////////////////////////////////
  ///// transerve shower shape ////
  /////////////////////////////////

  if(TRANSVERSE_PROFILE_HIST) {
    
    d_transverse = oFile->mkdir("transverse_shower_shape");
    d_transverse->cd();
    d_track_seed_diff = d_transverse->mkdir("track_seed_difference");
    for(int i = 0; i < 40; i++) {
      d_track_seed_diff->cd();
      sprintf(h_name,"h_track_seed_diff_x_L%02d",i+1);
      h_track_seed_diff_x[i] = new TH1F(h_name,h_name,400,-20.0,20.0);
      sprintf(h_name,"h_track_seed_diff_y_L%02d",i+1);
      h_track_seed_diff_y[i] = new TH1F(h_name,h_name,400,-20.0,20.0);
      sprintf(h_name,"h_track_seed_diff_dR_L%02d",i+1);
      h_track_seed_diff_dR[i] = new TH1F(h_name,h_name,620,-1.0,30.0);
   
    }

    d_track_cog_diff = d_transverse->mkdir("track_cog_difference");
    for(int i = 0; i < 40; i++) {
      d_track_cog_diff->cd();
      sprintf(h_name,"h_track_cog_diff_x_L%02d",i+1);
      h_track_cog_diff_x[i] = new TH1F(h_name,h_name,400,-20.0,20.0);
      sprintf(h_name,"h_track_cog_diff_y_L%02d",i+1);
      h_track_cog_diff_y[i] = new TH1F(h_name,h_name,400,-20.0,20.0);
      sprintf(h_name,"h_track_cog_diff_dR_L%02d",i+1);
      h_track_cog_diff_dR[i] = new TH1F(h_name,h_name,620,-1.0,30.0);
   
    }

    // d_transverse_Edep = d_transverse->mkdir("transverse_E_deposit");
    // d_transverse_Edep->cd();
    // for(int i = 0; i < 40; i++) {
    //   sprintf(h_name,"transerve_SS_%02d",i+1);
    //   d_transverse_SS[i] = d_transverse_Edep->mkdir(h_name);
    //   if(!d_transverse_SS[i]) { cout<<"Not found dir :: "<<h_name<<endl; }
    //   d_transverse_SS[i]->cd();
    //   for(int j = 0; j < 79; j++) {
    // 	sprintf(h_name,"h_transverse_dR0p56_SS%02d_L%02d",i+1,j+1);
    // 	h_transverse_dR0p56[i][j] = new TH1F(h_name,h_name,1000,0,4000);
    // 	sprintf(h_name,"h_transverse_dR1_SS%02d_L%02d",i+1,j+1);
    // 	h_transverse_dR1[i][j] = new TH1F(h_name,h_name,1000,0,4000);
    // 	sprintf(h_name,"h_transverse_dR2_SS%02d_L%02d",i+1,j+1);
    // 	h_transverse_dR2[i][j] = new TH1F(h_name,h_name,1000,0,4000);
    // 	sprintf(h_name,"h_transverse_dR3_SS%02d_L%02d",i+1,j+1);
    // 	h_transverse_dR3[i][j] = new TH1F(h_name,h_name,1000,0,4000);
    // 	sprintf(h_name,"h_transverse_dR5_SS%02d_L%02d",i+1,j+1);
    // 	h_transverse_dR5[i][j] = new TH1F(h_name,h_name,1000,0,4000);
    // 	sprintf(h_name,"h_transverse_dR8_SS%02d_L%02d",i+1,j+1);
    // 	h_transverse_dR8[i][j] = new TH1F(h_name,h_name,1000,0,4000);
    // 	sprintf(h_name,"h_transverse_dR10_SS%02d_L%02d",i+1,j+1);
    // 	h_transverse_dR10[i][j] = new TH1F(h_name,h_name,1000,0,4000);
    // 	sprintf(h_name,"h_transverse_dR12_SS%02d_L%02d",i+1,j+1);
    // 	h_transverse_dR12[i][j] = new TH1F(h_name,h_name,1000,0,4000);
    // 	sprintf(h_name,"h_transverse_dR18_SS%02d_L%02d",i+1,j+1);
    // 	h_transverse_dR18[i][j] = new TH1F(h_name,h_name,1000,0,4000);
    // 	sprintf(h_name,"h_transverse_dR20_SS%02d_L%02d",i+1,j+1);
    // 	h_transverse_dR20[i][j] = new TH1F(h_name,h_name,1000,0,4000);

      
    //   }
    // }


    // d_transverse_prof = d_transverse->mkdir("transverse_profile");
    // d_transverse_prof->cd();
    // for(int i = 0; i < 40; i++) {
    //   sprintf(h_name,"transerve_SS_%02d",i+1);
    //   d_transverse_SS_prof[i] = d_transverse_prof->mkdir(h_name);
    //   if(!d_transverse_SS_prof[i]) { cout<<"Not found dir :: "<<h_name<<endl; }
    //   sprintf(h_name,"h_transverse_prof_EE_SS%02d",i+1);
    //   h_transverse_prof_EE[i] = new TProfile(h_name,h_name,24,0.0,24.0);
    //   sprintf(h_name,"h_transverse_prof_FH_SS%02d",i+1);
    //   h_transverse_prof_FH[i] = new TProfile(h_name,h_name,24,0.0,24.0);
    
    //   d_transverse_SS_prof[i]->cd();
    //   for(int j = 0; j < 79; j++) {
    // 	sprintf(h_name,"h_transverse_prof_SS%02d_L%02d",i+1,j+1);
    // 	h_transverse_prof_layer[i][j] = new TProfile(h_name,h_name,24,0.0,24.0);
    //   }
    // }

  
    // d_trans_distance_all = d_transverse->mkdir("transverse_distance_all_SS");
    // d_trans_distance_all->cd();
  
    // for(int i = 0; i < 40; i++) {
    //   sprintf(SS_text,"SS_%02d",i+1);
    //   sprintf(h_name,"distance_%s",SS_text);
    //   d_transverse_SS_dis_all[i] = d_trans_distance_all->mkdir(h_name);
    //   d_transverse_SS_dis_all[i]->cd();
    
    //   for(int j = 0; j < 79; j++) {
      
    // 	sprintf(h_name,"h_transverse_%s_L%02d_distance_allSS",SS_text,j+1);
    // 	h_transverse_distance_allSS[i][j] = new TH1F(h_name,h_name,500,0.0,50.0);
    // 	h_transverse_distance_allSS[i][j]->Sumw2();

      
    //   }
    // }

    ////////////////////////////////
    ///// Transverse, inclusive  ///
    ////////////////////////////////

    d_transverse_inclusive = oFile->mkdir("transverse_inclusive");
    d_transverse_inclusive->cd();

    // d_trans_mips = d_transverse_inclusive->mkdir("transverse_mips");
    // d_trans_mips->cd();
  
    // for(int i = 0; i < 16; i++) {
    //   switch(i) {
    //   case 0: sprintf(SS_text,"SS_01_07"); break;
    //   case 1: sprintf(SS_text,"SS_08_14"); break;
    //   case 2: sprintf(SS_text,"SS_15_21"); break;
    //   case 3: sprintf(SS_text,"SS_22_28"); break;
    //   default: sprintf(SS_text,"SS_%02d",i+1+24); break;
    //   }
    //   sprintf(h_name,"transverse_%s",SS_text);
    //   d_transverse_SS_inc_mips[i] = d_trans_mips->mkdir(h_name);
    //   d_transverse_SS_inc_mips[i]->cd();
    
    //   for(int j = 0; j < 79; j++) {
      
    // 	sprintf(h_name,"h_transverse_dR0p56_%s_L%02d_mips",SS_text,j+1);
    // 	h_transverse_dR0p56_inc_mips[i][j] = new TH1F(h_name,h_name,1000,0,4000);
    // 	sprintf(h_name,"h_transverse_dR1_%s_L%02d_mips",SS_text,j+1);
    // 	h_transverse_dR1_inc_mips[i][j] = new TH1F(h_name,h_name,1000,0,4000);
    // 	sprintf(h_name,"h_transverse_dR2_%s_L%02d_mips",SS_text,j+1);
    // 	h_transverse_dR2_inc_mips[i][j] = new TH1F(h_name,h_name,1000,0,4000);
    // 	sprintf(h_name,"h_transverse_dR3_%s_L%02d_mips",SS_text,j+1);
    // 	h_transverse_dR3_inc_mips[i][j] = new TH1F(h_name,h_name,1000,0,4000);
    // 	sprintf(h_name,"h_transverse_dR5_%s_L%02d_mips",SS_text,j+1);
    // 	h_transverse_dR5_inc_mips[i][j] = new TH1F(h_name,h_name,1000,0,4000);
    // 	sprintf(h_name,"h_transverse_dR8_%s_L%02d_mips",SS_text,j+1);
    // 	h_transverse_dR8_inc_mips[i][j] = new TH1F(h_name,h_name,1000,0,4000);
    // 	sprintf(h_name,"h_transverse_dR10_%s_L%02d_mips",SS_text,j+1);
    // 	h_transverse_dR10_inc_mips[i][j] = new TH1F(h_name,h_name,1000,0,4000);
    // 	sprintf(h_name,"h_transverse_dR12_%s_L%02d_mips",SS_text,j+1);
    // 	h_transverse_dR12_inc_mips[i][j] = new TH1F(h_name,h_name,1000,0,4000);
    // 	sprintf(h_name,"h_transverse_dR18_%s_L%02d_mips",SS_text,j+1);
    // 	h_transverse_dR18_inc_mips[i][j] = new TH1F(h_name,h_name,1000,0,4000);
    // 	sprintf(h_name,"h_transverse_dR20_%s_L%02d_mips",SS_text,j+1);
    // 	h_transverse_dR20_inc_mips[i][j] = new TH1F(h_name,h_name,1000,0,4000);

      
    //   }
    // }


    // d_trans_mips_summed = d_transverse_inclusive->mkdir("transverse_mips_summed");
    // d_trans_mips_summed->cd();

    // for(int i = 0; i < 16; i++) {
    //   switch(i) {
    //   case 0: sprintf(SS_text,"SS_01_07"); break;
    //   case 1: sprintf(SS_text,"SS_08_14"); break;
    //   case 2: sprintf(SS_text,"SS_15_21"); break;
    //   case 3: sprintf(SS_text,"SS_22_28"); break;
    //   default: sprintf(SS_text,"SS_%02d",i+1+24); break;
    //   }
    //   sprintf(h_name,"transverse_%s_mips",SS_text);
    //   d_transverse_SS_inc_mips_summed[i] = d_trans_mips_summed->mkdir(h_name);
    //   d_transverse_SS_inc_mips_summed[i]->cd();
    
    //   for(int j = 0; j < 16; j++) {
    // 	char* layer_text = new char[100];
    // 	switch(j) {
    // 	case 0: sprintf(layer_text,"L01_07"); break;
    // 	case 1: sprintf(layer_text,"L08_14"); break;
    // 	case 2: sprintf(layer_text,"L15_21"); break;
    // 	case 3: sprintf(layer_text,"L22_28"); break;
    // 	default: sprintf(layer_text,"L%02d",j+1+24); break;
    // 	}

    // 	sprintf(h_name,"h_transverse_dR0p56_%s_%s_mips_summed",SS_text,layer_text);
    // 	h_transverse_dR0p56_inc_mips_summed[i][j] = new TH1F(h_name,h_name,2000,0.0,8000);
    // 	sprintf(h_name,"h_transverse_dR1_%s_%s_mips_summed",SS_text,layer_text);
    // 	h_transverse_dR1_inc_mips_summed[i][j] = new TH1F(h_name,h_name,2000,0.0,8000);
    // 	sprintf(h_name,"h_transverse_dR2_%s_%s_mips_summed",SS_text,layer_text);
    // 	h_transverse_dR2_inc_mips_summed[i][j] = new TH1F(h_name,h_name,2000,0.0,8000);
    // 	sprintf(h_name,"h_transverse_dR3_%s_%s_mips_summed",SS_text,layer_text);
    // 	h_transverse_dR3_inc_mips_summed[i][j] = new TH1F(h_name,h_name,2000,0.0,8000);
    // 	sprintf(h_name,"h_transverse_dR5_%s_%s_mips_summed",SS_text,layer_text);
    // 	h_transverse_dR5_inc_mips_summed[i][j] = new TH1F(h_name,h_name,2000,0.0,8000);
    // 	sprintf(h_name,"h_transverse_dR8_%s_%s_mips_summed",SS_text,layer_text);
    // 	h_transverse_dR8_inc_mips_summed[i][j] = new TH1F(h_name,h_name,2000,0.0,8000);
    // 	sprintf(h_name,"h_transverse_dR10_%s_%s_mips_summed",SS_text,layer_text);
    // 	h_transverse_dR10_inc_mips_summed[i][j] = new TH1F(h_name,h_name,2000,0.0,8000);
    // 	sprintf(h_name,"h_transverse_dR12_%s_%s_mips_summed",SS_text,layer_text);
    // 	h_transverse_dR12_inc_mips_summed[i][j] = new TH1F(h_name,h_name,2000,0.0,8000);
    // 	sprintf(h_name,"h_transverse_dR18_%s_%s_mips_summed",SS_text,layer_text);
    // 	h_transverse_dR18_inc_mips_summed[i][j] = new TH1F(h_name,h_name,2000,0.0,8000);
    // 	sprintf(h_name,"h_transverse_dR20_%s_%s_mips_summed",SS_text,layer_text);
    // 	h_transverse_dR20_inc_mips_summed[i][j] = new TH1F(h_name,h_name,2000,0.0,8000);
      
    //   }
    // }


  
    // d_trans_gev = d_transverse_inclusive->mkdir("transverse_gev");
    // d_trans_gev->cd();
  
    // for(int i = 0; i < 16; i++) {
    //   switch(i) {
    //   case 0: sprintf(SS_text,"SS_01_07"); break;
    //   case 1: sprintf(SS_text,"SS_08_14"); break;
    //   case 2: sprintf(SS_text,"SS_15_21"); break;
    //   case 3: sprintf(SS_text,"SS_22_28"); break;
    //   default: sprintf(SS_text,"SS_%02d",i+1+24); break;
    //   }
    //   sprintf(h_name,"transverse_%s",SS_text);
    //   d_transverse_SS_inc_gev[i] = d_trans_gev->mkdir(h_name);
    //   d_transverse_SS_inc_gev[i]->cd();
    
    //   for(int j = 0; j < 79; j++) {
      
    // 	sprintf(h_name,"h_transverse_dR0p56_%s_L%02d_gev",SS_text,j+1);
    // 	h_transverse_dR0p56_inc_gev[i][j] = new TH1F(h_name,h_name,500,0,50);
    // 	sprintf(h_name,"h_transverse_dR1_%s_L%02d_gev",SS_text,j+1);
    // 	h_transverse_dR1_inc_gev[i][j] = new TH1F(h_name,h_name,500,0,50);
    // 	sprintf(h_name,"h_transverse_dR2_%s_L%02d_gev",SS_text,j+1);
    // 	h_transverse_dR2_inc_gev[i][j] = new TH1F(h_name,h_name,500,0,50);
    // 	sprintf(h_name,"h_transverse_dR3_%s_L%02d_gev",SS_text,j+1);
    // 	h_transverse_dR3_inc_gev[i][j] = new TH1F(h_name,h_name,500,0,50);
    // 	sprintf(h_name,"h_transverse_dR5_%s_L%02d_gev",SS_text,j+1);
    // 	h_transverse_dR5_inc_gev[i][j] = new TH1F(h_name,h_name,500,0,50);
    // 	sprintf(h_name,"h_transverse_dR8_%s_L%02d_gev",SS_text,j+1);
    // 	h_transverse_dR8_inc_gev[i][j] = new TH1F(h_name,h_name,500,0,50);
    // 	sprintf(h_name,"h_transverse_dR10_%s_L%02d_gev",SS_text,j+1);
    // 	h_transverse_dR10_inc_gev[i][j] = new TH1F(h_name,h_name,500,0,50);
    // 	sprintf(h_name,"h_transverse_dR12_%s_L%02d_gev",SS_text,j+1);
    // 	h_transverse_dR12_inc_gev[i][j] = new TH1F(h_name,h_name,500,0,50);
    // 	sprintf(h_name,"h_transverse_dR18_%s_L%02d_gev",SS_text,j+1);
    // 	h_transverse_dR18_inc_gev[i][j] = new TH1F(h_name,h_name,500,0,50);
    // 	sprintf(h_name,"h_transverse_dR20_%s_L%02d_gev",SS_text,j+1);
    // 	h_transverse_dR20_inc_gev[i][j] = new TH1F(h_name,h_name,500,0,50);

      
    //   }
    // }


    // d_trans_gev_summed = d_transverse_inclusive->mkdir("transverse_gev_summed");
    // d_trans_gev_summed->cd();

    // for(int i = 0; i < 16; i++) {
    //   switch(i) {
    //   case 0: sprintf(SS_text,"SS_01_07"); break;
    //   case 1: sprintf(SS_text,"SS_08_14"); break;
    //   case 2: sprintf(SS_text,"SS_15_21"); break;
    //   case 3: sprintf(SS_text,"SS_22_28"); break;
    //   default: sprintf(SS_text,"SS_%02d",i+1+24); break;
    //   }
    //   sprintf(h_name,"transverse_%s_gev",SS_text);
    //   d_transverse_SS_inc_gev_summed[i] = d_trans_gev_summed->mkdir(h_name);
    //   d_transverse_SS_inc_gev_summed[i]->cd();
    
    //   for(int j = 0; j < 16; j++) {
    // 	char* layer_text = new char[100];
    // 	switch(j) {
    // 	case 0: sprintf(layer_text,"L01_07"); break;
    // 	case 1: sprintf(layer_text,"L08_14"); break;
    // 	case 2: sprintf(layer_text,"L15_21"); break;
    // 	case 3: sprintf(layer_text,"L22_28"); break;
    // 	default: sprintf(layer_text,"L%02d",j+1+24); break;
    // 	}

    // 	sprintf(h_name,"h_transverse_dR0p56_%s_%s_gev_summed",SS_text,layer_text);
    // 	h_transverse_dR0p56_inc_gev_summed[i][j] = new TH1F(h_name,h_name,2000,0,200);
    // 	sprintf(h_name,"h_transverse_dR1_%s_%s_gev_summed",SS_text,layer_text);
    // 	h_transverse_dR1_inc_gev_summed[i][j] = new TH1F(h_name,h_name,2000,0,200);
    // 	sprintf(h_name,"h_transverse_dR2_%s_%s_gev_summed",SS_text,layer_text);
    // 	h_transverse_dR2_inc_gev_summed[i][j] = new TH1F(h_name,h_name,2000,0,200);
    // 	sprintf(h_name,"h_transverse_dR3_%s_%s_gev_summed",SS_text,layer_text);
    // 	h_transverse_dR3_inc_gev_summed[i][j] = new TH1F(h_name,h_name,2000,0,200);
    // 	sprintf(h_name,"h_transverse_dR5_%s_%s_gev_summed",SS_text,layer_text);
    // 	h_transverse_dR5_inc_gev_summed[i][j] = new TH1F(h_name,h_name,2000,0,200);
    // 	sprintf(h_name,"h_transverse_dR8_%s_%s_gev_summed",SS_text,layer_text);
    // 	h_transverse_dR8_inc_gev_summed[i][j] = new TH1F(h_name,h_name,2000,0,200);
    // 	sprintf(h_name,"h_transverse_dR10_%s_%s_gev_summed",SS_text,layer_text);
    // 	h_transverse_dR10_inc_gev_summed[i][j] = new TH1F(h_name,h_name,2000,0,200);
    // 	sprintf(h_name,"h_transverse_dR12_%s_%s_gev_summed",SS_text,layer_text);
    // 	h_transverse_dR12_inc_gev_summed[i][j] = new TH1F(h_name,h_name,2000,0,200);
    // 	sprintf(h_name,"h_transverse_dR18_%s_%s_gev_summed",SS_text,layer_text);
    // 	h_transverse_dR18_inc_gev_summed[i][j] = new TH1F(h_name,h_name,2000,0,200);
    // 	sprintf(h_name,"h_transverse_dR20_%s_%s_gev_summed",SS_text,layer_text);
    // 	h_transverse_dR20_inc_gev_summed[i][j] = new TH1F(h_name,h_name,2000,0,200);
      
    //   }
    // }


    // d_trans_gev_summed_ratio = d_transverse_inclusive->mkdir("transverse_gev_summed_ratio");
    // d_trans_gev_summed_ratio->cd();

    // for(int i = 0; i < 16; i++) {
    //   switch(i) {
    //   case 0: sprintf(SS_text,"SS_01_07"); break;
    //   case 1: sprintf(SS_text,"SS_08_14"); break;
    //   case 2: sprintf(SS_text,"SS_15_21"); break;
    //   case 3: sprintf(SS_text,"SS_22_28"); break;
    //   default: sprintf(SS_text,"SS_%02d",i+1+24); break;
    //   }
    //   sprintf(h_name,"transverse_%s_gev",SS_text);
    //   d_transverse_SS_inc_gev_summed_ratio[i] = d_trans_gev_summed_ratio->mkdir(h_name);
    //   d_transverse_SS_inc_gev_summed_ratio[i]->cd();
    
    //   for(int j = 0; j < 16; j++) {
    // 	char* layer_text = new char[100];
    // 	switch(j) {
    // 	case 0: sprintf(layer_text,"L01_07"); break;
    // 	case 1: sprintf(layer_text,"L08_14"); break;
    // 	case 2: sprintf(layer_text,"L15_21"); break;
    // 	case 3: sprintf(layer_text,"L22_28"); break;
    // 	default: sprintf(layer_text,"L%02d",j+1+24); break;
    // 	}

    // 	sprintf(h_name,"h_transverse_dR0p56_%s_%s_gev_summed_ratio",SS_text,layer_text);
    // 	h_transverse_dR0p56_inc_gev_summed_ratio[i][j] = new TH1F(h_name,h_name,220,0.0,1.1);
    // 	sprintf(h_name,"h_transverse_dR1_%s_%s_gev_summed_ratio",SS_text,layer_text);
    // 	h_transverse_dR1_inc_gev_summed_ratio[i][j] = new TH1F(h_name,h_name,220,0.0,1.1);
    // 	sprintf(h_name,"h_transverse_dR2_%s_%s_gev_summed_ratio",SS_text,layer_text);
    // 	h_transverse_dR2_inc_gev_summed_ratio[i][j] = new TH1F(h_name,h_name,220,0.0,1.1);
    // 	sprintf(h_name,"h_transverse_dR3_%s_%s_gev_summed_ratio",SS_text,layer_text);
    // 	h_transverse_dR3_inc_gev_summed_ratio[i][j] = new TH1F(h_name,h_name,220,0.0,1.1);
    // 	sprintf(h_name,"h_transverse_dR5_%s_%s_gev_summed_ratio",SS_text,layer_text);
    // 	h_transverse_dR5_inc_gev_summed_ratio[i][j] = new TH1F(h_name,h_name,220,0.0,1.1);
    // 	sprintf(h_name,"h_transverse_dR8_%s_%s_gev_summed_ratio",SS_text,layer_text);
    // 	h_transverse_dR8_inc_gev_summed_ratio[i][j] = new TH1F(h_name,h_name,220,0.0,1.1);
    // 	sprintf(h_name,"h_transverse_dR10_%s_%s_gev_summed_ratio",SS_text,layer_text);
    // 	h_transverse_dR10_inc_gev_summed_ratio[i][j] = new TH1F(h_name,h_name,220,0.0,1.1);
    // 	sprintf(h_name,"h_transverse_dR12_%s_%s_gev_summed_ratio",SS_text,layer_text);
    // 	h_transverse_dR12_inc_gev_summed_ratio[i][j] = new TH1F(h_name,h_name,220,0.0,1.1);
    // 	sprintf(h_name,"h_transverse_dR18_%s_%s_gev_summed_ratio",SS_text,layer_text);
    // 	h_transverse_dR18_inc_gev_summed_ratio[i][j] = new TH1F(h_name,h_name,220,0.0,1.1);
    // 	sprintf(h_name,"h_transverse_dR20_%s_%s_gev_summed_ratio",SS_text,layer_text);
    // 	h_transverse_dR20_inc_gev_summed_ratio[i][j] = new TH1F(h_name,h_name,220,0.0,1.1);
      
    //   }
    // }


  
    // d_trans_ratio = d_transverse_inclusive->mkdir("transverse_ratio");
    // d_trans_ratio->cd();

    // for(int i = 0; i < 16; i++) {
    //   switch(i) {
    //   case 0: sprintf(SS_text,"SS_01_07"); break;
    //   case 1: sprintf(SS_text,"SS_08_14"); break;
    //   case 2: sprintf(SS_text,"SS_15_21"); break;
    //   case 3: sprintf(SS_text,"SS_22_28"); break;
    //   default: sprintf(SS_text,"SS_%02d",i+1+24); break;
    //   }
    //   sprintf(h_name,"transverse_%s",SS_text);
    //   d_transverse_SS_inc_ratio[i] = d_trans_ratio->mkdir(h_name);
    //   d_transverse_SS_inc_ratio[i]->cd();
    
    //   for(int j = 0; j < 79; j++) {
      
    // 	sprintf(h_name,"h_transverse_dR0p56_dR5_%s_L%02d_ratio",SS_text,j+1);
    // 	h_transverse_dR0p56_dR5_inc_ratio[i][j] = new TH1F(h_name,h_name,220,0.0,1.1);
    // 	sprintf(h_name,"h_transverse_dR1_dR5_%s_L%02d_ratio",SS_text,j+1);
    // 	h_transverse_dR1_dR5_inc_ratio[i][j] = new TH1F(h_name,h_name,220,0.0,1.1);
    // 	sprintf(h_name,"h_transverse_dR2_dR5_%s_L%02d_ratio",SS_text,j+1);
    // 	h_transverse_dR2_dR5_inc_ratio[i][j] = new TH1F(h_name,h_name,220,0.0,1.1);
    // 	sprintf(h_name,"h_transverse_dR3_dR5_%s_L%02d_ratio",SS_text,j+1);
    // 	h_transverse_dR3_dR5_inc_ratio[i][j] = new TH1F(h_name,h_name,220,0.0,1.1);

      
    //   }
    // }


    // d_trans_ratio_summed = d_transverse_inclusive->mkdir("transverse_ratio_summed");
    // d_trans_ratio_summed->cd();

    // for(int i = 0; i < 16; i++) {
    //   switch(i) {
    //   case 0: sprintf(SS_text,"SS_01_07"); break;
    //   case 1: sprintf(SS_text,"SS_08_14"); break;
    //   case 2: sprintf(SS_text,"SS_15_21"); break;
    //   case 3: sprintf(SS_text,"SS_22_28"); break;
    //   default: sprintf(SS_text,"SS_%02d",i+1+24); break;
    //   }
    //   sprintf(h_name,"transverse_%s",SS_text);
    //   d_transverse_SS_inc_ratio_summed[i] = d_trans_ratio_summed->mkdir(h_name);
    //   d_transverse_SS_inc_ratio_summed[i]->cd();
    
    //   for(int j = 0; j < 55; j++) {
    // 	char* layer_text = new char[100];
    // 	switch(j) {
    // 	case 0: sprintf(layer_text,"L01_07"); break;
    // 	case 1: sprintf(layer_text,"L08_14"); break;
    // 	case 2: sprintf(layer_text,"L15_21"); break;
    // 	case 3: sprintf(layer_text,"L22_28"); break;
    // 	default: sprintf(layer_text,"L%02d",j+1+24); break;
    // 	}

    // 	sprintf(h_name,"h_transverse_dR0p56_dR5_%s_%s_ratio_summed",SS_text,layer_text);
    // 	h_transverse_dR0p56_dR5_inc_ratio_summed[i][j] = new TH1F(h_name,h_name,220,0.0,1.1);
    // 	sprintf(h_name,"h_transverse_dR1_dR5_%s_%s_ratio_summed",SS_text,layer_text);
    // 	h_transverse_dR1_dR5_inc_ratio_summed[i][j] = new TH1F(h_name,h_name,220,0.0,1.1);
    // 	sprintf(h_name,"h_transverse_dR2_dR5_%s_%s_ratio_summed",SS_text,layer_text);
    // 	h_transverse_dR2_dR5_inc_ratio_summed[i][j] = new TH1F(h_name,h_name,220,0.0,1.1);
    // 	sprintf(h_name,"h_transverse_dR3_dR5_%s_%s_ratio_summed",SS_text,layer_text);
    // 	h_transverse_dR3_dR5_inc_ratio_summed[i][j] = new TH1F(h_name,h_name,220,0.0,1.1);

      
    //   }
    // }


    // d_trans_distance = d_transverse_inclusive->mkdir("transverse_distance");
    // d_trans_distance->cd();
  
    // for(int i = 0; i < 16; i++) {
    //   switch(i) {
    //   case 0: sprintf(SS_text,"SS_01_07"); break;
    //   case 1: sprintf(SS_text,"SS_08_14"); break;
    //   case 2: sprintf(SS_text,"SS_15_21"); break;
    //   case 3: sprintf(SS_text,"SS_22_28"); break;
    //   default: sprintf(SS_text,"SS_%02d",i+1+24); break;
    //   }
    //   sprintf(h_name,"distance_%s",SS_text);
    //   d_transverse_SS_dis[i] = d_trans_distance->mkdir(h_name);
    //   d_transverse_SS_dis[i]->cd();
    
    //   sprintf(h_name,"h_transverse_%s_distance_inc",SS_text);
    //   h_transverse_distance_inclusive[i] = new TH1F(h_name,h_name,500,0.0,50.0);
    //   h_transverse_distance_inclusive[i]->Sumw2();
    
    //   for(int j = 0; j < 79; j++) {
    // 	if(j < 3) {
    // 	  sprintf(h_name,"h_transverse_%s_distance_inc_compart_%d",SS_text,j+1);
    // 	  h_transverse_distance_inclusive_compart[i][j] = new TH1F(h_name,h_name,500,0.0,100.0);
    // 	  h_transverse_distance_inclusive_compart[i][j]->Sumw2();

    // 	}
    // 	sprintf(h_name,"h_transverse_%s_L%02d_distance",SS_text,j+1);
    // 	h_transverse_distance[i][j] = new TH1F(h_name,h_name,500,0.0,50.0);
    // 	h_transverse_distance[i][j]->Sumw2();

      
    //   }
    // }

    // d_trans_distance_zoom = d_transverse_inclusive->mkdir("transverse_distance_zoom");
    // d_trans_distance_zoom->cd();
  
    // for(int i = 0; i < 16; i++) {
    //   switch(i) {
    //   case 0: sprintf(SS_text,"SS_01_07"); break;
    //   case 1: sprintf(SS_text,"SS_08_14"); break;
    //   case 2: sprintf(SS_text,"SS_15_21"); break;
    //   case 3: sprintf(SS_text,"SS_22_28"); break;
    //   default: sprintf(SS_text,"SS_%02d",i+1+24); break;
    //   }
    //   sprintf(h_name,"distance_%s",SS_text);
    //   d_transverse_SS_dis_zoom[i] = d_trans_distance_zoom->mkdir(h_name);
    //   d_transverse_SS_dis_zoom[i]->cd();
    
    //   sprintf(h_name,"h_transverse_%s_distance_inc_zoom",SS_text);
    //   h_transverse_distance_inclusive_zoom[i] = new TH1F(h_name,h_name,50,0.0,5.0);
    //   h_transverse_distance_inclusive_zoom[i]->Sumw2();
    
    //   for(int j = 0; j < 79; j++) {
    // 	if(j < 3) {
    // 	  sprintf(h_name,"h_transverse_%s_distance_inc_compart_%d_zoom",SS_text,j+1);
    // 	  h_transverse_distance_inclusive_compart_zoom[i][j] = new TH1F(h_name,h_name,50,0.0,5.0);
    // 	  h_transverse_distance_inclusive_compart_zoom[i][j]->Sumw2();

    // 	}
    // 	sprintf(h_name,"h_transverse_%s_L%02d_distance_zoom",SS_text,j+1);
    // 	h_transverse_distance_zoom[i][j] = new TH1F(h_name,h_name,50,0.0,5.0);
    // 	h_transverse_distance_zoom[i][j]->Sumw2();

      
    //   }
    // }



    // d_trans_distance_prof = d_transverse_inclusive->mkdir("transverse_distance_prof");
    // d_trans_distance_prof->cd();
  
    // for(int i = 0; i < 16; i++) {
    //   switch(i) {
    //   case 0: sprintf(SS_text,"SS_01_07"); break;
    //   case 1: sprintf(SS_text,"SS_08_14"); break;
    //   case 2: sprintf(SS_text,"SS_15_21"); break;
    //   case 3: sprintf(SS_text,"SS_22_28"); break;
    //   default: sprintf(SS_text,"SS_%02d",i+1+24); break;
    //   }
    //   sprintf(h_name,"distance_%s",SS_text);
    //   d_transverse_SS_dis_prof[i] = d_trans_distance_prof->mkdir(h_name);
    //   d_transverse_SS_dis_prof[i]->cd();
    
    //   sprintf(h_name,"h_transverse_%s_distance_inc_prof",SS_text);
    //   h_transverse_distance_inclusive_prof[i] = new TProfile(h_name,h_name,500,0.0,50.0);
    //   //h_transverse_distance_inclusive[i]->Sumw2();
    
    //   for(int j = 0; j < 79; j++) {
    // 	if(j < 3) {
    // 	  sprintf(h_name,"h_transverse_%s_distance_inc_compart_%d_prof",SS_text,j+1);
    // 	  h_transverse_distance_inclusive_compart_prof[i][j] = new TProfile(h_name,h_name,500,0.0,50.0);
	

    // 	}
    // 	sprintf(h_name,"h_transverse_%s_L%02d_distance_prof",SS_text,j+1);
    // 	h_transverse_distance_prof[i][j] = new TProfile(h_name,h_name,500,0.0,50.0);


      
    //   }
    // }


    d_transverse_inclusive->cd();
    //p_htemp_check14 = new TProfile("p_htemp_check14","SS < 7 @L14 frac",500,0,50);
    // p_htemp_check38 = new TProfile("p_htemp_check38","SS = 30 @L38 frac",500,0,50);
    

    d_trans_prof_frac = d_transverse_inclusive->mkdir("transverse_prof_frac");
    d_trans_prof_frac->cd();

    
    for(int i = 0; i < 3; i++) {
      if(i == 0) sprintf(SS_text,"SS_03_07");
      else if(i == 1) sprintf(SS_text,"SS_29");
      else sprintf(SS_text,"SS_37");
      d_trans_prof_frac_SS[i] = d_trans_prof_frac->mkdir(SS_text);
      d_trans_prof_frac_SS[i]->cd();

      sprintf(h_name,"h_transverse_%s_prof_frac_EE",SS_text);
      h_transverse_prof_frac_Inc[i][0] = new TProfile(h_name,h_name,500,0.0,50.0);
      
      sprintf(h_name,"h_transverse_%s_prof_frac_FH",SS_text);
      h_transverse_prof_frac_Inc[i][1] = new TProfile(h_name,h_name,500,0.0,50.0);
      
      sprintf(h_name,"h_transverse_%s_prof_frac_AH",SS_text);
      h_transverse_prof_frac_Inc[i][2] = new TProfile(h_name,h_name,500,0.0,100.0);
      
      for(int j = 0; j < 79; j++) {
	sprintf(h_name,"h_transverse_%s_prof_frac_L%d",SS_text,j+1);
	h_transverse_prof_frac[i][j] = new TProfile(h_name,h_name,500,0.0,50.0);
      }
    }

    ///////////////////////////////////////

  }

  ///////////////////////////////////////////
  // E N D    O F    T R A N S V E R S E ///
  ///////////////////////////////////////////







    //h_trackX_trackY_EE1 = new TH2F("h_trackX_trackY_EE1","h_trackX_trackY_EE1",16,-8.0,8.0,16,-8.0,8.0);

    d_trackInfo = oFile->mkdir("track_information");
    d_trackInfo->cd();
    for(int i = 0; i < 40; i++) {

      /* sprintf(hname,"HGCAL_layer_%02d",i+1); */
      /* d_trackLayer[i] = d_trackInfo->mkdir(hname); */
    
      sprintf(hname,"h_trackX_layer_%02d",i+1);
      h_trackX[i] = new TH1F(hname,hname,400,-10.0,10.0);
      sprintf(hname,"h_trackY_layer_%02d",i+1);
      h_trackY[i] = new TH1F(hname,hname,400,-10.0,10.0);
      sprintf(hname,"h_trackX_trackY_layer_%02d",i+1);
      //h_trackX_trackY[i] = new TH2F(hname,hname,100,-10.0,10.0,100,-10.0,10.0);
      h_trackX_trackY[i] = new TH2F(hname,hname,400,-10.0,10.0,400,-10.0,10.0);
      h_trackX_trackY[i]->GetXaxis()->SetTitle("TrackX\'");
      h_trackX_trackY[i]->GetYaxis()->SetTitle("TrackY\'");


    }


    d_COG = oFile->mkdir("cog");
    d_COG->cd();
    for(int i = 0; i < 40; i++) {
      char* temp = new char[200];
      sprintf(temp,"h_cogX_layer%d",i+1);
      h_cogX[i] = new TH1F(temp,temp,400,-10.0,10.0);
      sprintf(temp,"h_cogY_layer%d",i+1);
      h_cogY[i] = new TH1F(temp,temp,400,-10.0,10.0);
      sprintf(temp,"h_cogX_cogY_layer%d",i+1);
      h_cogX_cogY[i] = new TH2F(temp,temp,400,-10.0,10.0,400,-10.0,10.0);
      h_cogX_cogY[i]->GetXaxis()->SetTitle("cogX");
      h_cogX_cogY[i]->GetYaxis()->SetTitle("cogY");
    }

  
    d_lateral = oFile->mkdir("lateral");
    d_lateral->cd();
    h_E1_E7_EE = new TH1F("h_E1_E7_EE","h_E1_E7_EE",220.0,0.0,1.1);
    h_E7_E19_EE = new TH1F("h_E7_E19_EE","h_E7_E19_EE",220.0,0.0,1.1);
    h_E1_E19_EE = new TH1F("h_E1_E19_EE","h_E1_E19_EE",220.0,0.0,1.1);

    h_E1_E7_FH = new TH1F("h_E1_E7_FH","h_E1_E7_FH",220.0,0.0,1.1);
    h_E7_E19_FH = new TH1F("h_E7_E19_FH","h_E7_E19_FH",220.0,0.0,1.1);
    h_E1_E19_FH = new TH1F("h_E1_E19_FH","h_E1_E19_FH",220.0,0.0,1.1);

    h_E1_E7_HGCAL = new TH1F("h_E1_E7_HGCAL","h_E1_E7_HGCAL",220.0,0.0,1.1);
    h_E7_E19_HGCAL = new TH1F("h_E7_E19_HGCAL","h_E7_E19_HGCAL",220.0,0.0,1.1);
    h_E1_E19_HGCAL = new TH1F("h_E1_E19_HGCAL","h_E1_E19_HGCAL",220.0,0.0,1.1);

    h_S1_S9_AH = new TH1F("h_S1_S9_AH","h_S1_S9_AH",220.0,0.0,1.1);
    h_S9_S25_AH = new TH1F("h_S9_S25_AH","h_S9_S25_AH",220.0,0.0,1.1);
    h_S1_S25_AH = new TH1F("h_S1_S25_AH","h_S1_S25_AH",220.0,0.0,1.1);


    d_debug = d_lateral->mkdir("debug");
    d_debug->cd();
    h_E1_E7_SS_EE = new TH1F("h_E1_E7_SS_EE","h_E1_E7_SS_EE",220.0,0.0,1.1);
    h_E7_E19_SS_EE = new TH1F("h_E7_E19_SS_EE","h_E7_E19_SS_EE",220.0,0.0,1.1);
    h_E1_E19_SS_EE = new TH1F("h_E1_E19_SS_EE","h_E1_E19_SS_EE",220.0,0.0,1.1);

    h_E1_E7_MIPs_in_EE = new TH1F("h_E1_E7_MIPs_in_EE","h_E1_E7_MIPs_in_EE",220.0,0.0,1.1);
    h_E7_E19_MIPs_in_EE = new TH1F("h_E7_E19_MIPs_in_EE","h_E7_E19_MIPs_in_EE",220.0,0.0,1.1);
    h_E1_E19_MIPs_in_EE = new TH1F("h_E1_E19_MIPs_in_EE","h_E1_E19_MIPs_in_EE",220.0,0.0,1.1);


    h_E1_E7_SS_FH = new TH1F("h_E1_E7_SS_FH","h_E1_E7_SS_FH",220.0,0.0,1.1);
    h_E7_E19_SS_FH = new TH1F("h_E7_E19_SS_FH","h_E7_E19_SS_FH",220.0,0.0,1.1);
    h_E1_E19_SS_FH = new TH1F("h_E1_E19_SS_FH","h_E1_E19_SS_FH",220.0,0.0,1.1);

    h_E1_E7_MIPs_in_FH = new TH1F("h_E1_E7_MIPs_in_FH","h_E1_E7_MIPs_in_FH",220.0,0.0,1.1);
    h_E7_E19_MIPs_in_FH = new TH1F("h_E7_E19_MIPs_in_FH","h_E7_E19_MIPs_in_FH",220.0,0.0,1.1);
    h_E1_E19_MIPs_in_FH = new TH1F("h_E1_E19_MIPs_in_FH","h_E1_E19_MIPs_in_FH",220.0,0.0,1.1);


    h_S1_S9_AH_debug = new TH1F("h_S1_S9_AH_debug","h_S1_S9_AH_debug",220.0,0.0,1.1);
    h_S9_S25_AH_debug = new TH1F("h_S9_S25_AH_debug","h_S9_S25_AH_debug",220.0,0.0,1.1);
    h_S1_S25_AH_debug = new TH1F("h_S1_S25_AH_debug","h_S1_S25_AH_debug",220.0,0.0,1.1);

    h_S1_S9_AH_debug_inv = new TH1F("h_S1_S9_AH_debug_inv","h_S1_S9_AH_debug_inv",220.0,0.0,1.1);
    h_S9_S25_AH_debug_inv = new TH1F("h_S9_S25_AH_debug_inv","h_S9_S25_AH_debug_inv",220.0,0.0,1.1);
    h_S1_S25_AH_debug_inv = new TH1F("h_S1_S25_AH_debug_inv","h_S1_S25_AH_debug_inv",220.0,0.0,1.1);


    /* h_S1_S9_SS_AH = new TH1F("h_S1_S9_SS_AH","h_S1_S9_SS_AH",220.0,0.0,1.1); */
    /* h_S9_S16_SS_AH = new TH1F("h_S9_S16_SS_AH","h_S9_S16_SS_AH",220.0,0.0,1.1); */
    /* h_S1_S16_SS_AH = new TH1F("h_S1_S16_SS_AH","h_S1_S16_SS_AH",220.0,0.0,1.1); */

    /* h_S1_S9_MIPs_in_AH = new TH1F("h_S1_S9_MIPs_in_AH","h_S1_S9_MIPs_in_AH",220.0,0.0,1.1); */
    /* h_S9_S16_MIPs_in_AH = new TH1F("h_S9_S16_MIPs_in_AH","h_S9_S16_MIPs_in_AH",220.0,0.0,1.1); */
    /* h_S1_S16_MIPs_in_AH = new TH1F("h_S1_S16_MIPs_in_AH","h_S1_S16_MIPs_in_AH",220.0,0.0,1.1); */



    d_layerwise_distribution = oFile->mkdir("layerwise_distribution");
    d_layerwise_distribution->cd();
    d_ld_nrechit = d_layerwise_distribution->mkdir("LD_nrechit");
    d_ld_energy = d_layerwise_distribution->mkdir("LD_energy");


    for(int i = 0; i < 79; i++) {
      d_ld_nrechit->cd();
      char* htemp = new char[2000];
      sprintf(htemp,"h_ld_nrechit_L%d",i+1);
      h_ld_nrechit[i] = new TH1F(htemp,htemp,300,0,300);	
      d_ld_energy->cd();
      sprintf(htemp,"h_ld_energy_L%d",i+1);
      h_ld_energy[i] = new TH1F(htemp,htemp,1700,0,3400);



    }

    
    if(SHOWER_RECO_HIST) {
  

        ////////////////////////////
    //energy distribution
    ////////////////////////////

    d_rechitEn_layer = oFile->mkdir("energy_distribution");
    d_rechitEn_layer->cd();
    sprintf(h_name,"h_ShowerInEE_EE_fraction");
    h_ShowerInEE_EE_fraction = new TH1F(h_name,h_name,220,0.0,1.1);
    sprintf(h_name,"h_ShowerInEE_FH_fraction");
    h_ShowerInEE_FH_fraction = new TH1F(h_name,h_name,220,0.0,1.1);
    sprintf(h_name,"h_ShowerInEE_AH_fraction");
    h_ShowerInEE_AH_fraction = new TH1F(h_name,h_name,220,0.0,1.1);

    for(int i = 0; i < 40; i++) {
      sprintf(h_name,"SS_%02d",i+1);
      d_prof_layers[i] = d_rechitEn_layer->mkdir(h_name);
      d_prof_layers[i]->cd();
      sprintf(h_name,"h_EE_fraction_SS%02d",i+1);
      h_EE_fraction[i] = new TH1F(h_name,h_name,220,0.0,1.1);
      sprintf(h_name,"h_FH_fraction_SS%02d",i+1);
      h_FH_fraction[i] = new TH1F(h_name,h_name,220,0.0,1.1);
      sprintf(h_name,"h_AH_fraction_SS%02d",i+1);
      h_AH_fraction[i] = new TH1F(h_name,h_name,220,0.0,1.1);
    
      for(int j = 0; j < 79; j++) {
    sprintf(h_name,"h_rechit_energy_SS_%02d_layer_%02d",i+1,j+1);
    //h_Rechits_En_SS[i][j] = new TH1F(h_name,h_name,per_layer_xbin,0.0,per_layer_xmax);
    h_Rechits_En_SS[i][j] = new TH1F(h_name,h_name,200,0.0,2000.0);
      
      }
    }
    // Energy distribution for SS = -1
    d_prof_layers[40] = d_rechitEn_layer->mkdir("SS_-1");
    d_prof_layers[40]->cd();
    for(int j = 0; j < 79; j++) {
      //sprintf(h_name,"h_rechit_energy_SS_%02d_layer_%02d",i+1,j+1);
      sprintf(h_name,"h_rechit_energy_SS_-1_layer_%02d",j+1);
      //h_Rechits_En_SS[i][j] = new TH1F(h_name,h_name,per_layer_xbin,0.0,per_layer_xmax);
      h_Rechits_En_SS[40][j] = new TH1F(h_name,h_name,250.0,0.0,500.0);
      
    }

    ////////////////////////////
    //Nrechit distribution
    ////////////////////////////

    d_rechitNrechit_layer = oFile->mkdir("nrechit_distribution");
    d_rechitNrechit_layer->cd();
    for(int i = 0; i < 40; i++) {
      sprintf(h_name,"SS_%02d",i+1);
      d_prof_layers_Nrec[i] = d_rechitNrechit_layer->mkdir(h_name);
      d_prof_layers_Nrec[i]->cd();
      for(int j = 0; j < 79; j++) {
    sprintf(h_name,"h_nrechit_SS_%02d_layer_%02d",i+1,j+1);
    //h_Rechits_nrec_SS[i][j] = new TH1F(h_name,h_name,per_layer_xbin,0.0,per_layer_xmax);
    h_Rechits_nrec_SS[i][j] = new TH1F(h_name,h_name,200.0,0.0,200.0);
      
      }
    }
    // Energy distribution for SS = -1
    d_prof_layers_Nrec[40] = d_rechitNrechit_layer->mkdir("SS_-1");
    d_prof_layers_Nrec[40]->cd();
    for(int j = 0; j < 79; j++) {
      //sprintf(h_name,"h_rechit_energy_SS_%02d_layer_%02d",i+1,j+1);
      sprintf(h_name,"h_nrechit_energy_SS_-1_layer_%02d",j+1);
      //h_Rechits_nrec_SS[i][j] = new TH1F(h_name,h_name,per_layer_xbin,0.0,per_layer_xmax);
      h_Rechits_nrec_SS[40][j] = new TH1F(h_name,h_name,200.0,0.0,200.0);
      
    }

    /////////////////////////////////


    d_shower_energy_reco = oFile->mkdir("shower_energy_reco");
    d_shower_energy_reco->cd();

    h_EE_inclusive = new TH1F("h_EE_inclusive","EE mips",xbin,0.0,x_2D_max);
    h_EE_inclusive->GetXaxis()->SetTitle("Energy in units of MIPs");
    h_FH_inclusive = new TH1F("h_FH_inclusive","FH mips",xbin,0.0,x_2D_max);
    h_FH_inclusive->GetXaxis()->SetTitle("Energy in units of MIPs");
    h_AH_inclusive = new TH1F("h_AH_inclusive","AH mips",xbin,0.0,x_2D_max);
    h_AH_inclusive->GetXaxis()->SetTitle("Energy in units of MIPs");
    h_all_inclusive = new TH1F("h_all_inclusive","all mips",xbin,0.0,x_2D_max);
    h_all_inclusive->GetXaxis()->SetTitle("Energy in units of MIPs");
    h_inclusive_elecpi_scale = new TH1F("h_inclusive_elecpi_scale","h_inclusive_elecpi_scale",y_max,0.0,y_max);
    h_inclusive_elecpi_scale->GetXaxis()->SetTitle("Energy [GeV]");

  
    d_showering_in_EE = d_shower_energy_reco->mkdir("showering_in_EE");
    d_showering_in_EE->cd();
  
    h_EE_showerinEE = new TH1F("h_EE_showerinEE","EE mips",xbin,0.0,x_2D_max);
    h_EE_showerinEE->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_FH_showerinEE = new TH1F("h_FH_showerinEE","FH mips",xbin,0.0,x_2D_max);
    h_FH_showerinEE->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_AH_showerinEE = new TH1F("h_AH_showerinEE","AH mips",xbin,0.0,x_2D_max);
    h_AH_showerinEE->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_all_showerinEE = new TH1F("h_all_showerinEE","all mips",xbin,0.0,x_2D_max);
    h_all_showerinEE->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_showerinEE_elecpi_scale = new TH1F("h_showerinEE_elecpi_scale","h_showerinEE_elecpi_scale",y_max,0.0,y_max);
    h_showerinEE_elecpi_scale->GetXaxis()->SetTitle("Energy [GeV]");

    h_showerinEE_elecpi_scale_withoutAH = new TH1F("h_showerinEE_elecpi_scale_withoutAH","h_showerinEE_elecpi_scale_withoutAH",y_max,0.0,y_max);
    h_showerinEE_elecpi_scale_withoutAH->GetXaxis()->SetTitle("Energy [GeV]");



    d_MIPs_in_EE = d_shower_energy_reco->mkdir("MIPs_in_EE");
    d_MIPs_in_EE->cd();
  
    h_EE_MIPsinEE = new TH1F("h_EE_MIPsinEE","EE mips",xbin,0.0,x_2D_max);
    h_EE_MIPsinEE->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_FH_MIPsinEE = new TH1F("h_FH_MIPsinEE","FH mips",xbin,0.0,x_2D_max);
    h_FH_MIPsinEE->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_AH_MIPsinEE = new TH1F("h_AH_MIPsinEE","AH mips",xbin,0.0,x_2D_max);
    h_AH_MIPsinEE->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_all_MIPsinEE = new TH1F("h_all_MIPsinEE","all mips",xbin,0.0,x_2D_max);
    h_all_MIPsinEE->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_MIPsinEE_elecpi_scale = new TH1F("h_MIPsinEE_elecpi_scale","h_MIPsinEE_elecpi_scale",y_max,0.0,y_max);
    h_MIPsinEE_elecpi_scale->GetXaxis()->SetTitle("Energy(GeV)");

    h_MIPsinEE_elecpi_scale_withoutAH = new TH1F("h_MIPsinEE_elecpi_scale_withoutAH","h_MIPsinEE_elecpi_scale_withoutAH",y_max,0.0,y_max);
    h_MIPsinEE_elecpi_scale_withoutAH->GetXaxis()->SetTitle("Energy(GeV)");


    d_MIPs_in_FH = d_shower_energy_reco->mkdir("MIPs_in_FH");
    d_MIPs_in_FH->cd();
  
    h_EE_MIPsinFH = new TH1F("h_EE_MIPsinFH","EE mips",xbin,0.0,x_2D_max);
    h_EE_MIPsinFH->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_FH_MIPsinFH = new TH1F("h_FH_MIPsinFH","FH mips",xbin,0.0,x_2D_max);
    h_FH_MIPsinFH->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_AH_MIPsinFH = new TH1F("h_AH_MIPsinFH","AH mips",xbin,0.0,x_2D_max);
    h_AH_MIPsinFH->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_all_MIPsinFH = new TH1F("h_all_MIPsinFH","all mips",xbin,0.0,x_2D_max);
    h_all_MIPsinFH->GetXaxis()->SetTitle("Energy in units of MIPs");

    h_MIPsinFH_elecpi_scale = new TH1F("h_MIPsinFH_elecpi_scale","h_MIPsinFH_elecpi_scale",y_max,0.0,y_max);
    h_MIPsinFH_elecpi_scale->GetXaxis()->SetTitle("Energy(GeV)");



    // d_FH_debug = oFile->mkdir("FH_debug");
    // d_FH_debug->cd();

    // for(int i = 0; i < 12; i++) {
    //   for(int j = 0; j < 7; j++) {
    //     char* htemp = new char[2000];
    //     sprintf(htemp,"h_nrechit_FH_L%d_P%d",i+1,j+1);
    //     h_nrechit_FH[i][j] = new TH1F(htemp,htemp,130.0,0.0,130.0);
    //     sprintf(htemp,"h_energy_FH_L%d_P%d",i+1,j+1);
    //     h_energy_FH[i][j] = new TH1F(htemp,htemp,1000.0,0.0,2000.0);
      
    //   }
    // }
  



    d_selection_cut_check_showering_in_EE = oFile->mkdir("selection_cut_check_showering_in_EE");
    d_selection_cut_check_showering_in_EE->cd();
    h_baseline = new TH1F("h_baseline","h_baseline",y_max,0.0,y_max);
    h_baseline->GetXaxis()->SetTitle("Energy(GeV)");
    h_SS1_reject = new TH1F("h_SS1_reject","h_SS1_reject",y_max,0.0,y_max);
    h_SS1_reject->GetXaxis()->SetTitle("Energy(GeV)");

    h_SS2_reject = new TH1F("h_SS2_reject","h_SS2_reject",y_max,0.0,y_max);
    h_SS2_reject->GetXaxis()->SetTitle("Energy(GeV)");

    h_trackwindow = new TH1F("h_trackwindow","h_trackwindow",y_max,0.0,y_max);
    h_trackwindow->GetXaxis()->SetTitle("Energy(GeV)");

    h_trackwindow_SS1_reject = new TH1F("h_trackwindow_SS1_reject","h_trackwindow_SS1_reject",y_max,0.0,y_max);
    h_trackwindow_SS1_reject->GetXaxis()->SetTitle("Energy(GeV)");

    h_trackwindow_SS2_reject = new TH1F("h_trackwindow_SS2_reject","h_trackwindow_SS2_reject",y_max,0.0,y_max);
    h_trackwindow_SS2_reject->GetXaxis()->SetTitle("Energy(GeV)");


    d_selection_cut_check_MIPs_in_EE = oFile->mkdir("selection_cut_check_MIPs_in_EE");
    d_selection_cut_check_MIPs_in_EE->cd();
    h_baseline_FH = new TH1F("h_baseline_FH","h_baseline_FH",y_max,0.0,y_max);
    h_baseline_FH->GetXaxis()->SetTitle("Energy(GeV)");
    h_SS1_reject_FH = new TH1F("h_SS1_reject_FH","h_SS1_reject_FH",y_max,0.0,y_max);
    h_SS1_reject_FH->GetXaxis()->SetTitle("Energy(GeV)");

    h_SS2_reject_FH = new TH1F("h_SS2_reject_FH","h_SS2_reject_FH",y_max,0.0,y_max);
    h_SS2_reject_FH->GetXaxis()->SetTitle("Energy(GeV)");

    h_trackwindow_FH = new TH1F("h_trackwindow_FH","h_trackwindow_FH",y_max,0.0,y_max);
    h_trackwindow_FH->GetXaxis()->SetTitle("Energy(GeV)");

    h_trackwindow_SS1_reject_FH = new TH1F("h_trackwindow_SS1_reject_FH","h_trackwindow_SS1_reject_FH",y_max,0.0,y_max);
    h_trackwindow_SS1_reject_FH->GetXaxis()->SetTitle("Energy(GeV)");

    h_trackwindow_SS2_reject_FH = new TH1F("h_trackwindow_SS2_reject_FH","h_trackwindow_SS2_reject_FH",y_max,0.0,y_max);
    h_trackwindow_SS2_reject_FH->GetXaxis()->SetTitle("Energy(GeV)");


    d_selection_cut_chi2 = oFile->mkdir("selection_cut_chi2");
    d_selection_cut_chi2->cd();

    h_trackwindow_SS2_reject_chi2_all = new TH1F("h_trackwindow_SS2_reject_chi2_all","h_trackwindow_SS2_reject_chi2_all",y_max,0.0,y_max);
    h_trackwindow_SS2_reject_chi2_all->GetXaxis()->SetTitle("Energy(GeV)");

    h_trackwindow_SS2_reject_chi2_all_recoE = new TH1F("h_trackwindow_SS2_reject_chi2_all_recoE","h_trackwindow_SS2_reject_chi2_all_recoE",y_max,0.0,y_max);
    h_trackwindow_SS2_reject_chi2_all_recoE->GetXaxis()->SetTitle("Energy(GeV)");

    d_selcection_cut_EE = d_selection_cut_chi2->mkdir("chi2_EE");
    d_selcection_cut_EE->cd();

    h_baseline_chi2_EE = new TH1F("h_baseline_chi2_EE","h_baseline_chi2_EE",y_max,0.0,y_max);
    h_baseline_chi2_EE->GetXaxis()->SetTitle("Energy(GeV)");
    h_SS1_reject_chi2_EE = new TH1F("h_SS1_reject_chi2_EE","h_SS1_reject_chi2_EE",y_max,0.0,y_max);
    h_SS1_reject_chi2_EE->GetXaxis()->SetTitle("Energy(GeV)");

    h_SS2_reject_chi2_EE = new TH1F("h_SS2_reject_chi2_EE","h_SS2_reject_chi2_EE",y_max,0.0,y_max);
    h_SS2_reject_chi2_EE->GetXaxis()->SetTitle("Energy(GeV)");

    h_trackwindow_chi2_EE = new TH1F("h_trackwindow_chi2_EE","h_trackwindow_chi2_EE",y_max,0.0,y_max);
    h_trackwindow->GetXaxis()->SetTitle("Energy(GeV)");

    h_trackwindow_SS1_reject_chi2_EE = new TH1F("h_trackwindow_SS1_reject_chi2_EE","h_trackwindow_SS1_reject_chi2_EE",y_max,0.0,y_max);
    h_trackwindow_SS1_reject_chi2_EE->GetXaxis()->SetTitle("Energy(GeV)");

    h_trackwindow_SS2_reject_chi2_EE = new TH1F("h_trackwindow_SS2_reject_chi2_EE","h_trackwindow_SS2_reject_chi2_EE",y_max,0.0,y_max);
    h_trackwindow_SS2_reject_chi2_EE->GetXaxis()->SetTitle("Energy(GeV)");

    h_trackwindow_SS2_reject_chi2_EE_recoE = new TH1F("h_trackwindow_SS2_reject_chi2_EE_recoE","h_trackwindow_SS2_reject_chi2_EE_recoE",y_max,0.0,y_max);
    h_trackwindow_SS2_reject_chi2_EE_recoE->GetXaxis()->SetTitle("Energy(GeV)");
  
    d_selcection_cut_FH = d_selection_cut_chi2->mkdir("chi2_FH");
    d_selcection_cut_FH->cd();

    h_baseline_chi2_FH = new TH1F("h_baseline_chi2_FH","h_baseline_chi2_FH",y_max,0.0,y_max);
    h_baseline_chi2_FH->GetXaxis()->SetTitle("Energy(GeV)");

  
    h_trackwindow_chi2_FH = new TH1F("h_trackwindow_chi2_FH","h_trackwindow_chi2_FH",y_max,0.0,y_max);
    h_trackwindow->GetXaxis()->SetTitle("Energy(GeV)");

    h_trackwindow_chi2_FH_recoE = new TH1F("h_trackwindow_chi2_FH_recoE","h_trackwindow_chi2_FH_recoE",y_max,0.0,y_max);
    h_trackwindow->GetXaxis()->SetTitle("Energy(GeV)");

    h_trackwindow_SS2_reject_chi2_FH = new TH1F("h_trackwindow_SS2_reject_chi2_FH","h_trackwindow_SS2_reject_chi2_FH",y_max,0.0,y_max);
    h_trackwindow_SS2_reject_chi2_FH->GetXaxis()->SetTitle("Energy(GeV)");

    h_trackwindow_SS2_reject_chi2_FH_recoE = new TH1F("h_trackwindow_SS2_reject_chi2_FH_recoE","h_trackwindow_SS2_reject_chi2_FH_recoE",y_max,0.0,y_max);
    h_trackwindow_SS2_reject_chi2_FH_recoE->GetXaxis()->SetTitle("Energy(GeV)");


  }
  

  
}

void AnalyzeHGCOctTB::Alignment_Map_Init() {
  char* f_name = new char[200];
  sprintf(f_name,"../txt_maps/Alignment_Map.txt");
  std::ifstream in(f_name);
  //int layer;
  if(!in) {
    cout<<"ERROR => "<<f_name<<" Not found"<<endl;
    //return;
    exit(0);
  }

  std::pair<float,float> dx_dy;
  std::pair<int, std::pair<float,float> > temp;
  int layer;
  float dx,dy;
  while(in>>layer>>dx>>dy) {
    //run_layer = std::make_pair(run,layer);
    dx_dy = std::make_pair(dx,dy);
    temp = std::make_pair(layer,dx_dy);
    align_map.insert(temp);
  }

  std::cout<<"INFO: Alignment MAP initialized successfully!!!"<<endl;
}

void AnalyzeHGCOctTB::Noise_Map_Init() {
  char* f_name = new char[200];
  sprintf(f_name,"../txt_maps/Noise_Map.txt");
  std::ifstream in(f_name);
  //int layer;
  if(!in) {
    cout<<"ERROR => "<<f_name<<" Not found"<<endl;
    //return;
    exit(0);
  }
  std::pair<int,int> mod_chip;
  std::pair<std::pair<int,int>, float> temp;
  int layer,mod_id,mod_pos,chip;
  float noise;
  while(in>>layer>>mod_id>>mod_pos>>chip>>noise) {
    //run_layer = std::make_pair(run,layer);
    mod_chip = std::make_pair(mod_id,chip);
    temp = std::make_pair(mod_chip,noise);
    noise_map.insert(temp);
  }

  std::cout<<"INFO: Noise MAP initialized successfully!!!"<<endl;
}
void AnalyzeHGCOctTB::Module_Map_Init(const char* config) {
  char *f_name = new char[200];

  if(strcmp(config,"alpha")==0 || strcmp(config,"config1")==0) {
    sprintf(f_name,"../txt_maps/config_maps/moduleMAP_config1.txt");
    cout<<"\n\nINFO: Mapping module configuration ALPHA (oct10-oct17) "<<endl;
    cout<<"INFO: Mapping EE[0]/FH[1]::Layer #[1-40]::Position on Layer[0 for EE]&[1-7 for FH] consult figure for Daisy structure configuration!!!"<<endl;

  }
  else if(strcmp(config,"bravo")==0 || strcmp(config,"config2")==0) {
    sprintf(f_name,"../txt_maps/config_maps/moduleMAP_config2.txt");
    cout<<"\n\nINFO: Mapping module configuration BRAVO (17oct-22oct) "<<endl;
    cout<<"INFO: Mapping EE[0]/FH[1]::Layer #[1-40]::Position on Layer[0 for EE]&[1-7 for FH] consult figure for Daisy structure configuration!!!"<<endl;

  }
  else if(strcmp(config,"charlie")==0  || strcmp(config,"config3")==0) {
    sprintf(f_name,"../txt_maps/config_maps/moduleMAP_config3.txt");
    cout<<"\n\nINFO: Mapping module configuration CHARLIE (23Oct-4Nov) "<<endl;
    cout<<"INFO: Mapping EE[0]/FH[1]::Layer #[1-40]::Position on Layer[0 for EE]&[1-7 for FH] consult figure for Daisy structure configuration!!!"<<endl;

  }
  else {
    cout<<"\n\nERROR: Incorrect configuration entered "<<endl;
    cout<<" Allowed configuration :\n alpha = Configuration 1 (10Oct-17Nov) \n bravo = Configuration 2 (17Oct-22Oct) \n charlie = Configuration 3 (23Oct-17Nov)"<<endl;
    return;
    
  }

  std::ifstream in(f_name);
  if(!in) {
    cout<<"ERROR => "<<f_name<<" Not found"<<endl;
    //return;
    exit(0);
  }

  /* if(!in){ */
  /*   cout<<"Could not find "<<f_name<<endl; */
  /*   return; */
  /* } */
  int modID_, part_, layer_, pos_;
  cout<<"File name = "<<f_name<<endl;
  while(in>>modID_>>part_>>layer_>>pos_){
    std::pair<int, std::vector<int>> temp_pair;
    std::vector<int> temp_vector;
    temp_vector.push_back(part_);
    temp_vector.push_back(layer_);
    temp_vector.push_back(pos_);
    temp_pair = std::make_pair(modID_,temp_vector);
    module_map.insert(temp_pair);
  }

  cout<<"INFO: Module Mapping Done!!! "<<endl<<endl;


}


void AnalyzeHGCOctTB::Layer_Position_Init() {
  char *f_name = new char[200];
  sprintf(f_name,"../txt_maps/config1_lengths.txt");

  std::ifstream in(f_name);
  if(!in) {
    cout<<"ERROR => "<<f_name<<" Not found"<<endl;
    //return;   
    exit(0);
  }
  int layer_;
  float  cm_, x0, nuc_int_, pi_int_;

  cout<<"File name = "<<f_name<<endl;
  while(in>>layer_>>cm_>>x0>>nuc_int_>>pi_int_){
    std::pair<int, std::vector<float>> temp_pair;
    std::vector<float> temp_vector;
    temp_vector.push_back(cm_/10.0);
    temp_vector.push_back(x0);
    temp_vector.push_back(nuc_int_);
    temp_vector.push_back(pi_int_);

    temp_pair = std::make_pair(layer_,temp_vector);
    layer_positions.insert(temp_pair);
  }

  cout<<"INFO: Layer Position Mapping Done!!! "<<endl<<endl;


}

void AnalyzeHGCOctTB::Weight_Map_Init() {
  char *f_name = new char[200];
  bool UseSimWeights = false;
  if(UseSimWeights) 
    sprintf(f_name,"../txt_maps/relative_weights_forSim.txt");
  else
    sprintf(f_name,"../txt_maps/relative_weights_forData.txt");

  std::ifstream in(f_name);
  if(!in) { 
    cout<<"ERROR => "<<f_name<<" Not found"<<endl;
    //return;
    exit(0);
  }
  int beamEnergy;
  float  alpha, beta, gamma_ee, gamma_fh;
  
  cout<<"File name = "<<f_name<<endl;
  while(in>>beamEnergy>>alpha>>beta>>gamma_ee>>gamma_fh){
    std::pair<int, std::vector<float>> temp_pair;
    std::vector<float> temp_vector;
    temp_vector.push_back(alpha);
    temp_vector.push_back(beta);
    temp_vector.push_back(gamma_ee);
    temp_vector.push_back(gamma_fh);
    
    temp_pair = std::make_pair(beamEnergy,temp_vector);
    rel_weights.insert(temp_pair);
  }

  if(UseSimWeights)
    cout<<BOLDRED<<"IMPORTANT: Weight have initialized based on SIM SCAN "<<RESET<<endl<<endl;
  else
    cout<<BOLDRED<<"IMPORTANT: Weight have initialized based on DATA SCAN "<<RESET<<endl<<endl;


}

void AnalyzeHGCOctTB::Chi2_Weight_Map_Init() {
  char *f_name_EH = new char[200];
  char *f_name_H = new char[200];
  bool UseSimWeights = false;
  sprintf(f_name_EH,"../txt_maps/chi2_calib_factors_EH_hadrons_DATA.txt");
  sprintf(f_name_H,"../txt_maps/chi2_calib_factors_H_hadrons_DATA.txt");

  std::ifstream in_EH(f_name_EH);
  std::ifstream in_H(f_name_H);
  if(!in_EH) { 
    cout<<"ERROR => "<<f_name_EH<<" Not found"<<endl;
    //return;
    exit(0);
  }
  if(!in_H) { 
    cout<<"ERROR => "<<f_name_H<<" Not found"<<endl;
    //return;
    exit(0);
  }

  int beamEnergy;
  float  w1, w2, w3;
  
  cout<<"File name = "<<f_name_EH<<endl;
  while(in_EH>>beamEnergy>>w1>>w2>>w3){
    std::pair<int, std::vector<float>> temp_pair;
    std::vector<float> temp_vector;
    temp_vector.push_back(w1);
    temp_vector.push_back(w2);
    temp_vector.push_back(w3);

    
    temp_pair = std::make_pair(beamEnergy,temp_vector);
    chi2_weights_EH.insert(temp_pair);
  }
  cout<<BOLDGREEN<<"INFO: Chi2 calibration Map initialized for EH hadrons!! "<<RESET<<endl;

  beamEnergy = -1.0; w1 = -1.0; w2 = -1.0; w3 = -1.0; 
  

  cout<<"File name = "<<f_name_H<<endl;
  while(in_H>>beamEnergy>>w1>>w2>>w3){
    std::pair<int, std::vector<float>> temp_pair;
    std::vector<float> temp_vector;
    temp_vector.push_back(w1);
    temp_vector.push_back(w2);
    temp_vector.push_back(w3);

    
    temp_pair = std::make_pair(beamEnergy,temp_vector);
    chi2_weights_H.insert(temp_pair);
  }
  cout<<BOLDGREEN<<"INFO: Chi2 calibration Map initialized for H hadrons!! "<<RESET<<endl<<endl;

}


void AnalyzeHGCOctTB::Scaling_Factor_Init() {
  scaling_factor_EH[20] = 1.087;
  scaling_factor_EH[50] = 1.051;
  scaling_factor_EH[80] = 1.032;
  scaling_factor_EH[100] = 1.035;
  scaling_factor_EH[120] = 1.029;
  scaling_factor_EH[200] = 1.026;
  scaling_factor_EH[250] = 1.021;
  scaling_factor_EH[300] = 1.019;

  scaling_factor_H[20] = 1.093;
  scaling_factor_H[50] = 1.047;
  scaling_factor_H[80] = 1.028;
  scaling_factor_H[100] = 1.027;
  scaling_factor_H[120] = 1.016;
  scaling_factor_H[200] = 1.014;
  scaling_factor_H[250] = 1.012;
  scaling_factor_H[300] = 1.009;
  
}


// void AnalyzeHGCOctTB::offical_calib_init() {
//   char *f_name = new char[200];
//   sprintf(f_name,"/home/shubham/work/HGCAL/CERNTB/CERN_5_oct_2018/txt_files/official_calib.txt");
//   std::ifstream in(f_name);
//   if(!in){
//     cout<<"Could not find "<<f_name<<endl;
//     return;
//   }
//   int layer_, module_, chip_, channel_;
//   long en_chan;
//   float adc_;
//   while(in>>layer_>>module_>>chip_>>channel_>>adc_){
//     en_chan = chip_*1000+channel_;
//     std::pair<int, int> temp;
//     temp = std::make_pair(layer_,en_chan);
//     std::pair<std::pair<int,int>, float> temp1;
//     temp1 = std::make_pair(temp, adc_);
//     offical_calib_map.insert(temp1);
//   }
// }

// void AnalyzeHGCOctTB::my_calib_init() {
//   char *f_name = new char[200];
//   sprintf(f_name,"/home/shubham/work/HGCAL/CERNTB/CERN_5_oct_2018/txt_files/ADC_MIP_v11.txt");
//   std::ifstream in(f_name);
//   if(!in){
//     cout<<"Could not find "<<f_name<<endl;
//     return;
//   }
//   int layer_, chip_, channel_,entry;
//   long en_chan;
//   float adc_,chi2,mip_err;
//   while(in>>layer_>>chip_>>channel_>>adc_>>chi2>>mip_err>>entry){
//     en_chan = chip_*1000+channel_;
//     std::pair<int,int> temp;
//     temp = std::make_pair(layer_+1,en_chan);
//     std::pair<std::pair<int,int>, float> temp1;
//     temp1 = std::make_pair(temp, adc_);
//     my_calib_map.insert(temp1);
//   }
// }

AnalyzeHGCOctTB::AnalyzeHGCOctTB(const TString &inputFileList, const char *outFileName, const char* dataset, const char* config, const char* energy, bool shower, bool longi, bool trans) {

  TChain *tree = new TChain("pion_variables");


  // if( ! FillChain(tree, tree2, tree3, inputFileList) ) {
  if( ! FillChain(tree, inputFileList) ) {
    
    std::cerr << "Cannot get the tree " << std::endl;
  } else {
    std::cout << "Initiating analysis of dataset " << dataset << std::endl;
  }

  /* if( ! FillChain(tree, inputFileList) ) { */
  /*   std::cerr << "Cannot get the tree " << std::endl; */
  /* } else { */
  /*   std::cout << "Initiating analysis of dataset " << dataset << std::endl; */
  /* } */

  HGCNtupleVariables::Init(tree);
  // HGCNtupleVariables::Init(tree, tree2, tree3);

  SHOWER_RECO_HIST = shower;
  LONGI_PROFILE_HIST = longi;
  TRANSVERSE_PROFILE_HIST = trans;
  
  BookHistogram(outFileName, config, energy);
  Module_Map_Init(config);
  Alignment_Map_Init();
  Noise_Map_Init();
  Layer_Position_Init();
  Weight_Map_Init();
  Chi2_Weight_Map_Init();
  Scaling_Factor_Init();
  // offical_calib_init();
  // my_calib_init();
  
}
// Bool_t AnalyzeHGCOctTB::FillChain(TChain *chain, TChain *chain2, TChain *chain3, const TString &inputFileList) {
Bool_t AnalyzeHGCOctTB::FillChain(TChain *chain, const TString &inputFileList) {

  ifstream infile(inputFileList, ifstream::in);
  std::string buffer;

  if(!infile.is_open()) {
    std::cerr << "** ERROR: Can't open '" << inputFileList << "' for input" << std::endl;
    return kFALSE;
  }

  std::cout << "TreeUtilities : FillChain " << std::endl;
  while(1) {
    infile >> buffer;
    if(!infile.good()) break;
    //std::cout << "Adding tree from " << buffer.c_str() << std::endl;                                                              
    chain->Add(buffer.c_str());
    // chain2->Add(buffer.c_str());
    // chain3->Add(buffer.c_str());

  }
  std::cout << "No. of Entries in chain  : " << chain->GetEntries() << std::endl;
  // std::cout << "No. of Entries in chain2 : " << chain2->GetEntries() << std::endl;
  // std::cout << "No. of Entries in chain3 : " << chain3->GetEntries() << std::endl;

  return kTRUE;
}

Long64_t AnalyzeHGCOctTB::LoadTree(Long64_t entry) {

  // Set the environment to read one entry                                                                                          
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (!fChain->InheritsFrom(TChain::Class()))  return centry;
  TChain *chain = (TChain*)fChain;
  if (chain->GetTreeNumber() != fCurrent) {
    fCurrent = chain->GetTreeNumber();
    //    Notify();
  }

  // if (!fChain2) return -5;
  // Long64_t centry2 = fChain2->LoadTree(entry);
  // if (centry2 < 0) return centry2;
  // if (!fChain2->InheritsFrom(TChain::Class()))  return centry2;
  // TChain *chain2 = (TChain*)fChain2;
  // if (chain2->GetTreeNumber() != fCurrent) {
  //   fCurrent = chain->GetTreeNumber();
  //   //    Notify();
  // }

  // if (!fChain3) return -5;
  // Long64_t centry3 = fChain3->LoadTree(entry);
  // if (centry3 < 0) return centry3;
  // if (!fChain3->InheritsFrom(TChain::Class()))  return centry3;
  // TChain *chain3 = (TChain*)fChain3;
  // if (chain3->GetTreeNumber() != fCurrent) {
  //   fCurrent = chain->GetTreeNumber();
  //   //    Notify();
  // }
  
  //if (centry==centry2)
  
  return centry;
  
  // cout<<"centry = "<<centry<<endl;
  // if(centry>0)
  //   return centry;
  // else return -1;
}

AnalyzeHGCOctTB::~AnalyzeHGCOctTB() { 

  // if (!fChain || !fChain2) return;
  // delete fChain->GetCurrentFile();
  // delete fChain2->GetCurrentFile();
  // oFile->cd();
  // oFile->Write();
  // oFile->Close();

  TStopwatch sw_dest;
  sw_dest.Start();
  double time_dest = 0.0;
  
  if (!fChain) return;
  delete fChain->GetCurrentFile();
  oFile->cd();
  oFile->Write();
  oFile->Close();

  time_dest += sw_dest.RealTime();
  sw_dest.Stop();
  cout<<endl<<endl<<"** time spent in file ops at destructor : "<<time_dest<<endl<<endl;

}

#endif

/*  LocalWords:  Nrechit EE R1 FH GetXaxis SetTitle Sumw2 TH2F reg3 NRechits
 */
/*  LocalWords:  GetYaxis SetTitleOffset
 */
