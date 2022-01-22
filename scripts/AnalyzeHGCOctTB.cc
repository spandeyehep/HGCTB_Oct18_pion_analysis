#define AnalyzeHGCOctTB_cxx

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
#include <vector>
#include <cstring>
#include "AnalyzeHGCOctTB.h"
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include <math.h>
#include <TF1.h>

using namespace std;



// chip 3022,44,3028




int main(int argc, char* argv[])
{

  if (argc < 4) {
    cerr << "Please give 5 arguments " << "runList " << " " << "outputFileName" << " " << "dataset" <<" "  <<" " << "energy" <<" "<< endl;
    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];
  
  const char *energy = argv[4];
  // const char *threshold = argv[6];
  std::map<int,double> E_thres;
  const char *config          = "alpha";
  
  E_thres[20] = 12.0; E_thres[50] = 20.0; E_thres[80] = 25.0; E_thres[100] = 30.0; E_thres[120] = 30.0; E_thres[200] = 40.0; E_thres[250] = 40.0; E_thres[300] = 40.0;


  float ene = (float)(std::atoi(energy));
  float thres = (float)(E_thres[ene]);

  
  
  //AnalyzeHGCOctTB hgcOctTB(inputFileList, outFileName, data, config, energy, threshold);
  AnalyzeHGCOctTB hgcOctTB(inputFileList, outFileName, data, config, energy);
  cout << "dataset " << data << " " << endl;
  cout << "configuration " << config << " " << endl;
  cout << "energy " << energy << " " << endl;
  // cout << "threshold " << threshold << " " << endl;

  hgcOctTB.EventLoop(data,thres);
  return 0;
}

void AnalyzeHGCOctTB::EventLoop(const char *data, float ratio_thres) {


  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();
  // Long64_t nentries3 = fChain3->GetEntriesFast();
  Long64_t hgc_jentry = 0;

  cout << "nentries " << nentries << endl;
  cout << "Analyzing dataset " << data << " " << endl;

  Long64_t nbytes = 0, nb = 0;
  Long64_t nbytes2 = 0, nb2 = 0;
  Long64_t nbytes3 = 0, nb3 = 0;
  // Long64_t nbytes3 = 0, nb3 = 0;

  bool ScaleSim = true;
  // float ee_rescaling = 1.035;
  // float fh_rescaling = 1.095;
  // float ah_rescaling = 1.0;

  float ee_rescaling = 1.035;
  float fh_rescaling = 1.095;
  float ah_rescaling = 1.095;

  if(ScaleSim) {
    std::cout<<BOLDGREEN<<"Sim rescaling ON, EE with factor : "<<ee_rescaling<<" ; FH with factor : "<<fh_rescaling<<" ; AH with factor : "<<ah_rescaling<<RESET<<std::endl;
  }
  else {
    std::cout<<BOLDRED<<"Sim rescaling OFF"<<RESET<<std::endl;
    ee_rescaling = 1.0;
    fh_rescaling = 1.0;
    ah_rescaling = 1.0;
  }

  Long64_t region_1_classified_events = 0;
  Long64_t region_2_classified_events = 0;
  Long64_t region_3_classified_events = 0;
  Long64_t non_classified_events = 0;

  Long64_t hadronic_int = 0;
  Long64_t hard_hadronic_int = 0;
  Long64_t notFound_hadronic_int = 0;

  long total_leading_pi0 = 0;

  double TQ_SS_notFound = 0.0;
  double TQ_SS_EE = 0.0;
  double TQ_SS_FH = 0.0;

  double SP_SS_notFound = 0.0;
  double SP_SS_EE = 0.0;
  double SP_SS_FH = 0.0;

  double true_SS_notFound = 0.0;
  double true_SS_EE = 0.0;
  double true_SS_FH = 0.0;

  double total_SS = 0.0;

  double within1Layer_cog = 0.0;
  double within2Layer_cog = 0.0;
  double total_hadInt_cog = 0.0;

  double within1Layer_cog_EE = 0.0;
  double within2Layer_cog_EE = 0.0;
  double within3Layer_cog_EE = 0.0;
  double total_hadInt_cog_EE = 0.0;

  double within1Layer_cog_FH = 0.0;
  double within2Layer_cog_FH = 0.0;
  double total_hadInt_cog_FH = 0.0;


  double within1Layer = 0.0;
  double within2Layer = 0.0;
  double total_hadInt = 0.0;

  double within1Layer_EE = 0.0;
  double within2Layer_EE = 0.0;
  double within3Layer_EE = 0.0;
  double total_hadInt_EE = 0.0;

  double within1Layer_interface = 0.0;
  double within2Layer_interface = 0.0;
  double total_hadInt_interface = 0.0;

  double within1Layer_FH = 0.0;
  double within2Layer_FH = 0.0;
  double total_hadInt_FH = 0.0;


  

  int decade = 0;
  int ahc_zeroHit = 0;
  double ee_zeroHit = 0.0;

  TF1* f_EH_w1 = new TF1("f_EH_w1","[0]*[0]+[1]*[1]/sqrt(x)", 20, 300);
  f_EH_w1->FixParameter(0,1.025153);
  f_EH_w1->FixParameter(1,1.624093);
  TF1* f_EH_w2 = new TF1("f_EH_w2","[0]*[0]+[1]*[1]/sqrt(x)", 20, 300);
  f_EH_w2->FixParameter(0,0.958412);
  f_EH_w2->FixParameter(1,1.226055);
  TF1* f_EH_w3 = new TF1("f_EH_w3","sqrt([0]*[0]+[1]*[1]/x)", 20, 300);
  f_EH_w3->FixParameter(0,1.019201);
  f_EH_w3->FixParameter(1,3.796437);

  TF1* f_H_w1 = new TF1("f_H_w1","[0]*x", 20, 300);
  f_H_w1->FixParameter(0,0.0);
  TF1* f_H_w2 = new TF1("f_H_w2","[0]*[0]+[1]*[1]/sqrt(x)", 20, 300);
  f_H_w2->FixParameter(0,0.908137);
  f_H_w2->FixParameter(1,0.995139);
  TF1* f_H_w3 = new TF1("f_H_w3","sqrt([0]*[0]+[1]*[1]/x)", 20, 300);
  f_H_w3->FixParameter(0,0.977450);
  f_H_w3->FixParameter(1,2.996701);

  
  
  float lambda[79];
  for(int i = 0; i < 79; i++) {
    lambda[i] = layer_positions[i+1].at(2); // pick up lambda_int
  }

  

  //in mm
  float ahc_pos[39] = {27.45, 53.65, 79.85, 106.05, 132.25, 
		       158.45, 184.65,210.85,237.05, 263.25,
		       289.45, 315.65, 341.85, 368.05, 394.25,
		       420.45, 446.65, 472.85, 499.05, 525.25,
		       551.45, 577.65, 603.85, 630.05, 656.25, 
		       682.45, 708.65, 734.85, 761.05, 787.25,
		       813.45, 839.65, 865.85, 892.05, 918.25, 
		       944.45, 970.65, 996.85, 1075.45};
  //in cm
  float ahc_front = 169.9;


  bool DEBUG = false;
  int debug_count = 0;
  Long64_t cut_count[28];
  Long64_t nEvents = 0;
  Long64_t MIP_pions = 0;
  int TOTAL_ACTIVE_LAYER = -1;
  int EE_LAYER = -1;
  int FH_LAYER = -1;
  if(!strcmp(conf_,"alpha") || !strcmp(conf_,"config1")) {
    TOTAL_ACTIVE_LAYER = 40;
    EE_LAYER = 28;
    FH_LAYER = 12;
  }
  else if(!strcmp(conf_,"bravo") || !strcmp(conf_,"config2")){
    TOTAL_ACTIVE_LAYER = 39;
    EE_LAYER = 28;
    FH_LAYER = 11;
  }
  else if(!strcmp(conf_,"charlie") || !strcmp(conf_,"config3")) {
    TOTAL_ACTIVE_LAYER = 20;
    EE_LAYER = 8;
    FH_LAYER = 12;
  }
  else {
    cout<<"ERROR: Unknown configuration!!!!"<<endl;
    return;
  }



  double E_beam = -1.0;




  if(DEBUG) cout<<"DEBUG: Configuration = "<<conf_<<endl;
  if(DEBUG) cout<<"DEBUG: TOTAL_ACTIVE_LAYER = "<<TOTAL_ACTIVE_LAYER<<endl;
  if(DEBUG) cout<<"DEBUG: EE_LAYER = "<<EE_LAYER<<endl;
  if(DEBUG) cout<<"DEBUG: FH_LAYER = "<<FH_LAYER<<endl;

  if(DEBUG) cout<<"DEBUG: Entering event Loop"<<endl;
  for(int i = 0; i< 28; i++){ cut_count[i] = 0;}
  Long64_t jentry = 0;;
  for (jentry=0; jentry<nentries;jentry++) {
    // ==============print number of events done == == == == == == == =
    double progress = 10.0 * jentry / (1.0 * nentries);
    int k = int (progress);
    if (k > decade)
      cout << 10 * k << " %" << endl;
    decade = k;
    
    // ===============read this entry == == == == == == == == == == ==

    Long64_t ientry = LoadTree(jentry);
    if (ientry < 0) { break; cout<<"Breaking"<<endl;}
        // cout<<"****"<<endl;
    nb = fChain->GetEntry(jentry);   nbytes += nb;
    nb2 = fChain2->GetEntry(jentry); nbytes2 += nb2;
    nb3 = fChain3->GetEntry(jentry); nbytes3 += nb3;

    
    // ===============Synch HGCAL & AHCAL with same entry == == == == == == == == == == ==

    // while(run != runNumber) {
    //   hgc_jentry++;
    //   nb = fChain->GetEntry(hgc_jentry);   nbytes += nb;
    //   nb2 = fChain2->GetEntry(hgc_jentry); nbytes2 += nb2;
    // }



    if(NRechits == 0) continue;
    h_nTracks->Fill(ntracks);
    if(ntracks != 1) continue;
    E_beam = beamEnergy;
    
    h_particleID->Fill(pdgID);
    h_runNumber->Fill(run);
    h_beamEnergy->Fill(trueBeamEnergy);


    const int TOTAL_LAYERS = 79;
    Float_t track_x[TOTAL_LAYERS];
    Float_t track_y[TOTAL_LAYERS];

    track_x[0] = impactX_HGCal_layer_1;
    track_y[0] = impactY_HGCal_layer_1;
    track_x[1] = impactX_HGCal_layer_2;
    track_y[1] = impactY_HGCal_layer_2;
    track_x[2] = impactX_HGCal_layer_3;
    track_y[2] = impactY_HGCal_layer_3;
    track_x[3] = impactX_HGCal_layer_4;
    track_y[3] = impactY_HGCal_layer_4;
    track_x[4] = impactX_HGCal_layer_5;
    track_y[4] = impactY_HGCal_layer_5;
    track_x[5] = impactX_HGCal_layer_6;
    track_y[5] = impactY_HGCal_layer_6;
    track_x[6] = impactX_HGCal_layer_7;
    track_y[6] = impactY_HGCal_layer_7;
    track_x[7] = impactX_HGCal_layer_8;
    track_y[7] = impactY_HGCal_layer_8;
    track_x[8] = impactX_HGCal_layer_9;
    track_y[8] = impactY_HGCal_layer_9;
    track_x[9] = impactX_HGCal_layer_10;
    track_y[9] = impactY_HGCal_layer_10;
    track_x[10] = impactX_HGCal_layer_11;
    track_y[10] = impactY_HGCal_layer_11;
    track_x[11] = impactX_HGCal_layer_12;
    track_y[11] = impactY_HGCal_layer_12;
    track_x[12] = impactX_HGCal_layer_13;
    track_y[12] = impactY_HGCal_layer_13;
    track_x[13] = impactX_HGCal_layer_14;
    track_y[13] = impactY_HGCal_layer_14;
    track_x[14] = impactX_HGCal_layer_15;
    track_y[14] = impactY_HGCal_layer_15;
    track_x[15] = impactX_HGCal_layer_16;
    track_y[15] = impactY_HGCal_layer_16;
    track_x[16] = impactX_HGCal_layer_17;
    track_y[16] = impactY_HGCal_layer_17;
    track_x[17] = impactX_HGCal_layer_18;
    track_y[17] = impactY_HGCal_layer_18;
    track_x[18] = impactX_HGCal_layer_19;
    track_y[18] = impactY_HGCal_layer_19;
    track_x[19] = impactX_HGCal_layer_20;
    track_y[19] = impactY_HGCal_layer_20;
    track_x[20] = impactX_HGCal_layer_21;
    track_y[20] = impactY_HGCal_layer_21;
    track_x[21] = impactX_HGCal_layer_22;
    track_y[21] = impactY_HGCal_layer_22;
    track_x[22] = impactX_HGCal_layer_23;
    track_y[22] = impactY_HGCal_layer_23;
    track_x[23] = impactX_HGCal_layer_24;
    track_y[23] = impactY_HGCal_layer_24;
    track_x[24] = impactX_HGCal_layer_25;
    track_y[24] = impactY_HGCal_layer_25;
    track_x[25] = impactX_HGCal_layer_26;
    track_y[25] = impactY_HGCal_layer_26;
    track_x[26] = impactX_HGCal_layer_27;
    track_y[26] = impactY_HGCal_layer_27;
    track_x[27] = impactX_HGCal_layer_28;
    track_y[27] = impactY_HGCal_layer_28;
    track_x[28] = impactX_HGCal_layer_29;
    track_y[28] = impactY_HGCal_layer_29;
    track_x[29] = impactX_HGCal_layer_30;
    track_y[29] = impactY_HGCal_layer_30;
    track_x[30] = impactX_HGCal_layer_31;
    track_y[30] = impactY_HGCal_layer_31;
    track_x[31] = impactX_HGCal_layer_32;
    track_y[31] = impactY_HGCal_layer_32;
    track_x[32] = impactX_HGCal_layer_33;
    track_y[32] = impactY_HGCal_layer_33;
    track_x[33] = impactX_HGCal_layer_34;
    track_y[33] = impactY_HGCal_layer_34;
    track_x[34] = impactX_HGCal_layer_35;
    track_y[34] = impactY_HGCal_layer_35;
    track_x[35] = impactX_HGCal_layer_36;
    track_y[35] = impactY_HGCal_layer_36;
    track_x[36] = impactX_HGCal_layer_37;
    track_y[36] = impactY_HGCal_layer_37;
    track_x[37] = impactX_HGCal_layer_38;
    track_y[37] = impactY_HGCal_layer_38;
    track_x[38] = impactX_HGCal_layer_39;
    track_y[38] = impactY_HGCal_layer_39;
    track_x[39] = impactX_HGCal_layer_40;
    track_y[39] = impactY_HGCal_layer_40;

    for(int j = 0; j < 39; j++) {
      // double x = (rechit_z->at(j) - b_x)/m_x;
      // double y = (rechit_z->at(j) - b_y)/m_y;
      double z = ahc_front + (ahc_pos[j]/10.0) ;
      double x = (z * m_x) + b_x;
      double y = (z * m_y) + b_y;

      track_x[j+40] = x;
      track_y[j+40] = y;
      // cout<<"layer="<<rechit_layer->at(j)<<", z="<<rechit_z->at(j)
      // 	  <<" ;my (x,y)=("<<x<<","<<y<<") ;track(x,y)=("
      // 	  <<track_x[rechit_layer->at(j)-1]<<","
      // 	  <<track_y[rechit_layer->at(j)-1]<<")"<<endl;
    }


    Double_t rechitEnergySum = 0.0;
    Double_t un_cali = 0.0;
    Long_t Nrechit_layer[40];

    Long_t NRechits_EE[28];
    Long_t NRechits_FH[12];
    Long_t NRechits_AH[39];
    Long_t Nrechit_EE_debug = 0;
    Double_t RechitsEn_EE[28];
    Double_t RechitsEn_FH[12];
    Double_t RechitsEn_AH[39];
    Double_t RechitEn_layer[40];

    int module_part_ = -1;
    int module_layer_ = -1;
    int module_position_ = -1;

    double energy_sum_dR2_[40];
    double cogX_[40];
    double cogY_[40];
    
    for(int ii=0;ii<40;ii++){
      if(ii<28) {
      	NRechits_EE[ii]=0;
      	RechitsEn_EE[ii] = 0.0;
      }
      else{
      	NRechits_FH[ii-28]=0;
      	RechitsEn_FH[ii-28] = 0.0;
      }
      Nrechit_layer[ii]=0;
      energy_sum_dR2_[ii] = 0.0;
      RechitEn_layer[ii] = 0.0;
      cogX_[ii] = 0.0;
      cogY_[ii] = 0.0;
      if(ii < 39) {
      	RechitsEn_AH[ii] = 0.0;
      	NRechits_AH[ii] = 0;
      }

    }
    /// FIll Tree
    if(DEBUG) cout<<"DEBUG: Start Analylizing RecHits!!"<<endl;
    if(DEBUG) cout<<"DEBUG: NRechits = "<<NRechits<<endl;


    bool DoMIPrescaling = false;

    std::vector<int> temp_moduleID;
    int Nrech_L1[4] = {0,0,0,0};
    for(int i = 0 ; i < NRechits; i++){
      temp_moduleID.clear();
      int temp_layer = rechit_layer->at(i);
      int temp_chip = rechit_chip->at(i);
      int temp_channel = rechit_channel->at(i);
      int en_chan = temp_chip*1000+temp_channel;
      
      double recx = rechit_x->at(i);
      double recy = rechit_y->at(i);
      int rechit_modID = (int)rechit_module->at(i);
      std::pair<float, float> temp_mod_chip(rechit_modID,temp_chip);

      //rescaling
      if(temp_layer<=28) {
	rechit_energy->at(i) = rechit_energy->at(i)/ee_rescaling;
      }
      else {
	rechit_energy->at(i) = rechit_energy->at(i)/fh_rescaling;
      }

      
      //channel masking
      if(en_chan == 3022 || en_chan == 3028 || en_chan == 44) continue;
      if(temp_layer==1 && temp_chip==0) continue;
      if(temp_layer==1) Nrech_L1[temp_chip]++;
      
      // noise cut
      
      float noise_chip = getNoise(temp_mod_chip);
      if(DEBUG && false) cout<<"Module, layer, chip , noise = "<<rechit_module->at(i)<<", "<<rechit_layer->at(i)<<", "<<rechit_chip->at(i)<<", "<<noise_chip<<endl;
      if(temp_layer <= 28 && rechit_amplitudeHigh->at(i) < 3*noise_chip) continue;
      else if(temp_layer > 28 && rechit_amplitudeHigh->at(i) < 4*noise_chip) continue;


      
      float mip_ratio_chip = getMIPRatio(temp_mod_chip);
      bool Is200umCell = false;
      if(rechit_modID == 144 || rechit_modID == 145 || rechit_modID == 146 || rechit_modID == 147) Is200umCell = true;

      // // energy scaling using MIPs
      if(DoMIPrescaling) {
      	if(mip_ratio_chip > -1.0 && !Is200umCell) {
      	  rechit_energy->at(i) = rechit_energy->at(i)/mip_ratio_chip;
      	  if(rechit_modID == 144 || rechit_modID == 145 || rechit_modID == 146 || rechit_modID == 147) 
      	    cout<<"BUG found"<<endl;
      	}
      }

      //5% rescaling 
      // rechit_energy->at(i) = rechit_energy->at(i)*0.95;

      if(temp_layer <= 28) Nrechit_EE_debug++;
      
      double trackx = track_x[temp_layer-1];
      double tracky = track_y[temp_layer-1];
      if(temp_layer == 2)
	       dX_dY_layer_2->Fill(recx-trackx,recy-tracky);
      
      h_dX_dY_track_vs_rechit[temp_layer-1]->Fill(recx-trackx,recy-tracky);

      if(!strcmp(data,"data")) {
      	trackx = -1*trackx;
      	tracky = -1*tracky;
      	std::pair<float, float> dxy_al = dxy_alignment(temp_layer);
      	float dx_corr = dxy_al.first;
      	float dy_corr = dxy_al.second; 
      	recx -= dx_corr;
      	recy -= dy_corr;
      	
      	
      }
      Nrechit_layer[temp_layer-1]++;
      RechitEn_layer[temp_layer-1]+=rechit_energy->at(i);
      
      temp_moduleID = getModuleLocation(rechit_module->at(i));
      
      cogX_[temp_layer-1] += (rechit_x->at(i)*rechit_energy->at(i));
      cogY_[temp_layer-1] += (rechit_y->at(i)*rechit_energy->at(i));



      if((temp_layer==38 && rechit_module->at(i)!=62) || (temp_layer==39 && rechit_module->at(i)!=54) || (temp_layer==40 && rechit_module->at(i)!=43)){
      	cout<< jentry+1 << " = " << rechit_module->at(i)<<endl;
      }


      if(!temp_moduleID.size() || temp_moduleID.size()<3) {
      	cout<<"ERROR: Could NOT locate MODULE location for module "<<rechit_module->at(i)<<endl;
      	cout<<"\t more info temp_moduleID.size() = "<<temp_moduleID.size()<<endl;
      	return;
      }
      module_part_ = temp_moduleID.at(0);
      module_layer_ = temp_moduleID.at(1);
      module_position_ = temp_moduleID.at(2);
      
      
      h_moduleID->Fill(module_part_);
      


      
      
      /////////////////// For older shower finding algo ////////
      
      double dR = deltaR(recx,recy,trackx,tracky);
      // cout<<rechit_layer->at(i)<<" => rec(X,Y) = ("<<recx<<","<<recy
      // 	  <<"), track(X,Y) = ("<<trackx<<","<<tracky
      // 	  <<")\t dR = "<<dR<<endl;
      if(rechit_layer->at(i) > 28)
      {
        if(dR < 2.0) 	energy_sum_dR2_[temp_layer-1]+=rechit_energy->at(i);
      }
      else 
      {
        if(dR < 2.0)  energy_sum_dR2_[temp_layer-1]+=rechit_energy->at(i);
      }
      
            
      ///////////////////////////////////////////////
      
      
      if(module_part_ == 0) {
      	RechitsEn_EE[module_layer_-1]+=rechit_energy->at(i);
      	NRechits_EE[module_layer_-1]++;
      }
      else if(module_part_ == 1) {
      	RechitsEn_FH[module_layer_-1]+=rechit_energy->at(i);
      	NRechits_FH[module_layer_-1]++;
	
      }
      else {
      	cout<<"ERROR: Unknown Module Part detected!!!!"<<endl;
      	return;
      	
      }
      rechitEnergySum+=rechit_energy->at(i);
      
    }
    if(Nrechit_EE_debug == 0)  {
      continue;
    }

    ////////////////////////////////////////////
    //            AHCAL Part                  //
    ////////////////////////////////////////////
    
    Double_t rechitEnergySumAHCAL = sim_energyAH;
    // for(int i = 0 ; i < ahc_nHits; i++) {
    //   int temp_layer = ahc_hitK->at(i);
    //   RechitsEn_AH[temp_layer -1] += ahc_hitEnergy->at(i);
    //   rechitEnergySumAHCAL += ahc_hitEnergy->at(i);
    //   NRechits_AH[temp_layer -1]++;
    // }

    /////////////////////////////////////////////

    if(DEBUG) cout<<"DEBUG: For shower start"<<endl;
    if(DEBUG && false) {
      for (int ll = 0; ll < 40; ll++) {
      	cout<<ll+1<<"\t"<<energy_sum_dR2_[ll]<<endl;
      }
    }

    
    /////////////////// shower finding algo ////////
    bool MIP = true;
    int shower_start_index = -1;
    float shower_lambda_ = -1.0;
    float shower_weight_ = 1.0;
    float energy_thres = 20.0;
    if(energy_sum_dR2_[0] > energy_thres) {
      // cout<<"Shower started in layer =1"<<endl;
      shower_start_index = 0;
      shower_lambda_ = lambda[0];
      shower_weight_ = lambda[0];
      MIP = false;
    }
    
    else if(energy_sum_dR2_[1] > energy_thres && energy_sum_dR2_[1] > 2*energy_sum_dR2_[0]) {
      // cout<<"Shower started in layer =2"<<endl;
      shower_start_index = 1;
      shower_lambda_ = lambda[1];
      shower_weight_ = lambda[1]-lambda[0];
      MIP = false;
    }
    
    else {
      
      for(int i = 2; i < 40; i++) {
      	if((energy_sum_dR2_[i] > energy_thres) && (energy_sum_dR2_[i] > 2*energy_sum_dR2_[i-1]) && (energy_sum_dR2_[i] > 2*energy_sum_dR2_[i-2])) {

      	  // if(i+1==27 || i+1==28) {
      	  //   cout<<"Shower started in layer ="<<i+1<<" ,lambda_ = "<<lambda[i]<<endl;
      	  // }
      	  shower_start_index = i;
      	  shower_lambda_ = lambda[i];
      	  shower_weight_ = lambda[i]-lambda[i-1];
      	  MIP = false;
      	  break;
      	}
      }
    }


    int closest_layer = getNextLayer(HadInt_z);
    if(closest_layer > 40) closest_layer = -1;
    // cout<<"jentry : HadInt_z : closest_layer :: "<<jentry<<" : "<<" : "<<HadInt_z<<" : "<<closest_layer<<endl;

    if(pdgID == -13)  { closest_layer = 13; shower_start_index = 12;}

    
    // if(foundHadInt && closest_layer > 10)
    //   cout<<jentry+1<<" Shower Start : HadInt_z : NextClosestLayer is \t "<<shower_start_index<<" : "<<HadInt_z<<" : "<<closest_layer<<endl;

    
    // cout<<"\t First hadronic interaction :  "<<HadInt_z<<endl;

    double nsc_total_KE = 0.0;
    vector<double> KE_list;
    double leadingHadKE = 0.0;
    double E_gen_kin = 0.0;
    KE_list.clear();
    double E_gen_kin_thres = 0.0;
   
    vector<double> KE_pi0;
    int npi0 = 0;
    double maxKE = 0.0;
    double lead_pi0 = 0.0;
    double pi0_Esum = 0.0;
    if(!foundHadInt) notFound_hadronic_int++;
    if(foundHadInt) {
      hadronic_int++;
      for(int i = 0 ; i < nsec; i++) {
        nsc_total_KE += sec_kin->at(i);
        if(abs(sec_pdgID->at(i)) == 211 && sec_charge->at(i) == -1) KE_list.push_back(sec_kin->at(i));
	if(sec_pdgID->at(i) == 111) {npi0++; KE_pi0.push_back(sec_kin->at(i)); pi0_Esum += sec_kin->at(i);}

	if(sec_kin->at(i) > maxKE) maxKE = sec_kin->at(i);
      }
      leadingHadKE = getLeadingKE(KE_list);
      E_gen_kin = nsc_total_KE - leadingHadKE;


      //////////////////////////////
      // for neutral pions (pi0) ///
      //////////////////////////////
      //if(E_gen_kin/E_beam > 0.4 && pdgID != -13) {
      if(pdgID != -13) {
	//double lead_pi0 = getLeadingKE(KE_pi0);
	lead_pi0 = getLeadingKE(KE_pi0);
	if(lead_pi0 >= maxKE) {
	  total_leading_pi0++;
	  
	  h_pi0_energy->Fill(lead_pi0/trueBeamEnergy);
	  
	  if(lead_pi0/trueBeamEnergy > 0.98 && false)
	    cout<<jentry<<" : "<<npi0<< " : "<< " : "<<lead_pi0<<endl;
	  
	}
	h_npi0->Fill(npi0);
	h_pi0_energysum->Fill(pi0_Esum/trueBeamEnergy);
      }
      //////////////////////////////
      if(E_gen_kin/E_beam > 0.4 && pdgID != -13) {
            if(E_gen_kin > 10 && false) {
              cout<<jentry<<" NextClosestLayer is \t "<<closest_layer<<endl;
            }
            
           
            hard_hadronic_int++;
            
            if(jentry < 20 && false) {
              cout<<jentry<<" : "<<nsec<<" : "<<nsc_total_KE<< " : "<<leadingHadKE<<" : "<<E_gen_kin<<endl;
            }
            h_2D_nsc_correlation->Fill(E_gen_kin,nsec);
            h_nsc->Fill(nsec);
            h_seckin->Fill(nsc_total_KE);
            h_had_leading_KE->Fill(leadingHadKE);
            h_gen_kin->Fill(E_gen_kin);
            h_2D_nsc_correlation_frac->Fill(E_gen_kin/E_beam,nsec);
            h_2D_hlead_sec_frac->Fill(E_gen_kin/E_beam,leadingHadKE/E_beam);

            if(closest_layer >= 3 && closest_layer <= 28) {
              h_EE_layer_i->Fill(energy_sum_dR2_[closest_layer-1]);
              h_EE_layer_i_1->Fill(energy_sum_dR2_[closest_layer-2]);
              h_EE_layer_i_2->Fill(energy_sum_dR2_[closest_layer-3]);

              h_EE_layer_i_1_ratio->Fill(energy_sum_dR2_[closest_layer-1]/energy_sum_dR2_[closest_layer-2]);
              h_EE_layer_i_2_ratio->Fill(energy_sum_dR2_[closest_layer-1]/energy_sum_dR2_[closest_layer-3]);
            }
            else if(closest_layer > 28 && closest_layer <= 31) {
              h_FH_interface_layer_i->Fill(energy_sum_dR2_[closest_layer-1]);
              h_FH_interface_layer_i_1->Fill(energy_sum_dR2_[closest_layer-2]);
              h_FH_interface_layer_i_2->Fill(energy_sum_dR2_[closest_layer-3]);
              h_FH_interface_layer_i_1_ratio->Fill(energy_sum_dR2_[closest_layer-1]/energy_sum_dR2_[closest_layer-2]);
              h_FH_interface_layer_i_2_ratio->Fill(energy_sum_dR2_[closest_layer-1]/energy_sum_dR2_[closest_layer-3]);
            }
            else if(closest_layer >= 32 && closest_layer <= 40) {
             h_FH_layer_i->Fill(energy_sum_dR2_[closest_layer-1]);
             h_FH_layer_i_1->Fill(energy_sum_dR2_[closest_layer-2]);
             h_FH_layer_i_2->Fill(energy_sum_dR2_[closest_layer-3]); 
             h_FH_layer_i_1_ratio->Fill(energy_sum_dR2_[closest_layer-1]/energy_sum_dR2_[closest_layer-2]);
             h_FH_layer_i_2_ratio->Fill(energy_sum_dR2_[closest_layer-1]/energy_sum_dR2_[closest_layer-3]);
            }

            if(closest_layer >= 1 && closest_layer <= 40) {
              if(abs(closest_layer - (shower_start_index +1)) <= 1 && (shower_start_index +1) != 0) within1Layer++;
              if(abs(closest_layer - (shower_start_index +1)) <= 2 && (shower_start_index +1) != 0) within2Layer++;

              total_hadInt++;

              if(closest_layer >= 1 && closest_layer <= 28) {
                if(abs(closest_layer - (shower_start_index +1)) <= 1 && (shower_start_index +1) != 0) within1Layer_EE++;
                if(abs(closest_layer - (shower_start_index +1)) <= 2 && (shower_start_index +1) != 0) within2Layer_EE++;
                if(abs(closest_layer - (shower_start_index +1)) <= 3 && (shower_start_index +1) != 0) within3Layer_EE++;
                total_hadInt_EE++;
              }
              else if(closest_layer >= 29 && closest_layer <= 40) {
                if(abs(closest_layer - (shower_start_index +1)) <= 1) within1Layer_FH++;
                if(abs(closest_layer - (shower_start_index +1)) <= 2) within2Layer_FH++;
                else {
                  if(shower_start_index+1 <  closest_layer && false) {
                  bool cond1 = (energy_sum_dR2_[closest_layer-1] > 20);
                  bool cond2 = (energy_sum_dR2_[closest_layer-1] > 2*energy_sum_dR2_[closest_layer-2]);
                  bool cond3 = (energy_sum_dR2_[closest_layer-1] > 2*energy_sum_dR2_[closest_layer-2]);
                 cout<<jentry<<" : "<<closest_layer<<" : "<<(shower_start_index +1)<<" :: "<<cond1<<" : "<<cond2
                    <<" : "<<cond3<<endl;
                  cout<<"\tcond1 : cond2 : cond3 :: "<<energy_sum_dR2_[closest_layer-1]
                    <<" : "<<(energy_sum_dR2_[shower_start_index]/energy_sum_dR2_[closest_layer-2])
                    <<" : " <<(energy_sum_dR2_[shower_start_index]/energy_sum_dR2_[closest_layer-2])<<endl;
                  }
                }
                total_hadInt_FH++;
              }

            }
            
      }
      if(closest_layer >= 1 && closest_layer <= 40 )
	h_delta_layer->Fill(closest_layer - (shower_start_index +1));
      h_truth_reco_SS_corr->Fill(closest_layer,(shower_start_index +1));
      if(shower_start_index < 0) {
      	h_truth_reco_SS_corr_neg1->Fill(closest_layer,(shower_start_index));
      	h_delta_layer_neg1->Fill(closest_layer - (shower_start_index));
      }
      else {
      	h_truth_reco_SS_corr_neg1->Fill(closest_layer,(shower_start_index+1));
      	h_delta_layer_neg1->Fill(closest_layer - (shower_start_index +1));
      }
    }



    h_shower_start_check->Fill(shower_lambda_);
    if(closest_layer < 0.0) h_hadronic_Int->Fill(closest_layer);
    else h_hadronic_Int->Fill(lambda[closest_layer-1]);
	

    
    h_shower_start->Fill(shower_lambda_);
    h_shower_start_dN_dLambda->Fill(shower_lambda_,1/shower_weight_);

    ///////////////////////////////////////////////
    //                Collapsed EE               //
    ///////////////////////////////////////////////
    if(shower_start_index+1 <= 28 && (shower_start_index+1) > 0) {
      h_shower_start_full_collapsed_EE->Fill(lambda[0]);
      if((shower_start_index+1)%2 == 0) {
    	h_shower_start_part_collapsed_EE->Fill(lambda[shower_start_index-1]);
      }
      else {
    	h_shower_start_part_collapsed_EE->Fill(lambda[shower_start_index]);
      }
    }
    else {
      h_shower_start_full_collapsed_EE->Fill(shower_lambda_);
      h_shower_start_part_collapsed_EE->Fill(shower_lambda_);
    }

    if(DEBUG) cout<<"DEBUG: Shower start index = "<<shower_start_index<<endl;

    ///////////////////////////////////////////////

    

    Long_t Nrechit_EE = 0;
    Long_t Nrechit_FH = 0;
    Double_t rechitEnergySum_EE = 0.0;
    Double_t rechitEnergySum_FH = 0.0;
    Double_t test = 0.0;
    bool zero_rh[40];
    
    bool DoWeight = false;
    double beta_ = -1.0;

    switch((int)beamEnergy) {
    case 20: beta_ = 4.207; break;
    case 50: beta_ = 4.752; break;
    case 80: beta_ = 5.009; break;
    case 100: beta_ = 5.174; break;
    case 120: beta_ = 5.074; break;
    case 200: beta_ = 5.094; break;
    case 250: beta_ = 5.305; break;
    case 30: beta_ = 5.253; break;
    }
    if(!DoWeight) beta_ = 1.0;

    
    for(int iL = 0; iL < TOTAL_ACTIVE_LAYER; iL++){
      if(foundHadInt && lead_pi0 < maxKE/2) continue;
      if(iL < EE_LAYER) {
      	Nrechit_EE+=NRechits_EE[iL];
      	rechitEnergySum_EE+=RechitsEn_EE[iL];

      	cogX_[iL] = cogX_[iL]/RechitsEn_EE[iL];
      	cogY_[iL] = cogY_[iL]/RechitsEn_EE[iL];

      	h_Rechits_En_inclusive[iL]->Fill(RechitsEn_EE[iL]);
      	h_longi_profile_inclusive->Fill(lambda[iL],RechitsEn_EE[iL]);
      	if(!MIP) {
      	  h_longi_profile_raw[shower_start_index]->Fill(lambda[iL],RechitsEn_EE[iL]);
      	  h_longi_profile_region3->Fill(lambda[iL],RechitsEn_EE[iL]);
      	  h_Rechits_En_SS[shower_start_index][iL]->Fill(RechitsEn_EE[iL]);

      	}
      }
      else {
      	Nrechit_FH+=NRechits_FH[iL-EE_LAYER];
      	rechitEnergySum_FH+=RechitsEn_FH[iL-EE_LAYER];
      	cogX_[iL] = cogX_[iL]/RechitsEn_FH[iL-EE_LAYER];
      	cogY_[iL] = cogY_[iL]/RechitsEn_FH[iL-EE_LAYER];
      	h_Rechits_En_inclusive[iL]->Fill(RechitsEn_FH[iL-EE_LAYER]);
      	h_longi_profile_inclusive->Fill(lambda[iL],beta_*RechitsEn_FH[iL-EE_LAYER]);
      	if(!MIP) {
      	  h_longi_profile_raw[shower_start_index]->Fill(lambda[iL],beta_*RechitsEn_FH[iL-EE_LAYER]);
      	  h_Rechits_En_SS[shower_start_index][iL]->Fill(RechitsEn_FH[iL-EE_LAYER]);

      	}
      }
      
      h_trackX_trackY[iL]->Fill(track_x[iL],track_y[iL]);
      h_cogX_cogY[iL]->Fill(cogX_[iL],cogY_[iL]);
      h_cogX_trackX[iL]->Fill(cogX_[iL],track_x[iL]);
      h_cogY_trackY[iL]->Fill(cogY_[iL],track_y[iL]);
      h_dX_dY[iL]->Fill(cogX_[iL]-track_x[iL],cogY_[iL]-track_y[iL]);
      h_dX[iL]->Fill(cogX_[iL]-track_x[iL]);
      h_dY[iL]->Fill(cogY_[iL]-track_y[iL]);

      //cout<<iL+1<<":("<<track_x[iL]<<","<<track_y[iL]<<"):("<<cogX_[iL]<<","<<cogY_[iL]<<")"<<endl;
    }
    
   
    ////////////////////////////////////////////////////////////////////////////
    //////////////////// For shower start finder optimization    ///////////////
    ////////////////////////////////////////////////////////////////////////////

    double SS_indR2[40];
    double SS_indR5[40];
    double SS_indR10[40];
    double SS_outdR2[40];


    double SS_indR10_fraction[40];
    double SS_indR10_MW_fraction[40];
    double R210_avg[40];
    
    double SS_indR10_MW_total[40];

    int SS_inNrec_dR5[40];
    int SS_inNrec_dR10[40];
    //int SS_inNrec_dR5[40];
    for(int ii=0;ii<40;ii++){ SS_indR2[ii] = 0.0; SS_outdR2[ii] = 0.0; SS_indR5[ii] = 0.0; SS_indR10[ii] = 0.0; SS_inNrec_dR5[ii] = 0; SS_inNrec_dR10[ii] = 0; SS_indR10_fraction[ii] = 0.0; SS_indR10_MW_fraction[ii] = 0.0; SS_indR10_MW_total[ii] = 0.0; R210_avg[ii] = 0.0;}

    for(int i = 0 ; i < NRechits; i++){
      temp_moduleID.clear();
      int temp_layer = rechit_layer->at(i);
      int temp_chip = rechit_chip->at(i);
      int temp_channel = rechit_channel->at(i);
      int en_chan = temp_chip*1000+temp_channel;
      
      double recx = rechit_x->at(i);
      double recy = rechit_y->at(i);
      int rechit_modID = (int)rechit_module->at(i);
      std::pair<float, float> temp_mod_chip(rechit_modID,temp_chip);
      //channel masking
      if(en_chan == 3022 || en_chan == 3028 || en_chan == 44) continue;
      if(temp_layer==1 && temp_chip==0) continue;
     
      // noise cut
      float noise_chip = getNoise(temp_mod_chip);
      if(temp_layer <= 28 && rechit_amplitudeHigh->at(i) < 3*noise_chip) continue;
      else if(temp_layer > 28 && rechit_amplitudeHigh->at(i) < 4*noise_chip) continue;

      if(!strcmp(data,"data")) {
        std::pair<float, float> dxy_al = dxy_alignment(temp_layer);
        float dx_corr = dxy_al.first;
        float dy_corr = dxy_al.second; 
        recx -= dx_corr;
        recy -= dy_corr;
      }
      /////////////////// shower finding algo ////////
      
      double dR = deltaR(recx,recy,cogX_[temp_layer-1],cogY_[temp_layer-1]);


      if(dR < 2.0) SS_indR2[temp_layer-1] += rechit_energy->at(i);
      else SS_outdR2[temp_layer-1] += rechit_energy->at(i);
      if(dR < 5.0) { SS_indR5[temp_layer-1] += rechit_energy->at(i); SS_inNrec_dR5[temp_layer-1]++;}
      if(dR < 10.0) { SS_indR10[temp_layer-1] += rechit_energy->at(i); SS_inNrec_dR10[temp_layer-1]++;}

    }


    for(int ii=0; ii < 40; ii++) {
      if(ii < 28) {
        SS_indR10_fraction[ii] = SS_indR10[ii]/(rechitEnergySum_EE+rechitEnergySum_FH);
        if(ii ==  0) { 
          SS_indR10_MW_total[0] = (SS_indR10[0] + SS_indR10[1])/2;
          SS_indR10_MW_fraction[0] = SS_indR10_MW_total[0]/(rechitEnergySum_EE+rechitEnergySum_FH);
          
        }
        else if(ii == 27) { 
          SS_indR10_MW_total[ii] = (SS_indR10[ii-1] + SS_indR10[ii])/2;
          SS_indR10_MW_fraction[ii] = SS_indR10_MW_total[ii]/(rechitEnergySum_EE+rechitEnergySum_FH);
        }
        else { 
          SS_indR10_MW_total[ii] = (SS_indR10[ii-1] + SS_indR10[ii] + SS_indR10[ii+1])/3;
          SS_indR10_MW_fraction[ii] = SS_indR10_MW_total[ii]/(rechitEnergySum_EE+rechitEnergySum_FH);
        }

      }
      else {
        SS_indR10_fraction[ii] = SS_indR10[ii]/(rechitEnergySum_EE+rechitEnergySum_FH);
        if(ii-28 ==  0) { 
          SS_indR10_MW_total[ii] = (SS_indR10[ii] + SS_indR10[ii+1])/2;
          SS_indR10_MW_fraction[ii] = SS_indR10_MW_total[ii]/(rechitEnergySum_EE+rechitEnergySum_FH);
        }
        else if(ii-28 == 11) {
          SS_indR10_MW_total[ii] = (SS_indR10[ii-1] + SS_indR10[ii])/2;
          SS_indR10_MW_fraction[ii] = SS_indR10_MW_total[ii]/(rechitEnergySum_EE+rechitEnergySum_FH);
       }
        else { 
          SS_indR10_MW_total[ii] = (SS_indR10[ii-1] + SS_indR10[ii] + SS_indR10[ii+1])/3;
          SS_indR10_MW_fraction[ii] = SS_indR10_MW_total[ii]/(rechitEnergySum_EE+rechitEnergySum_FH);
        }

      }
      
      if(ii == 38) {
        double R2 = SS_indR2[ii] + SS_indR2[ii+1];
        double R10 = SS_indR10[ii] + SS_indR10[ii+1];
        if(R10 > 0.0) {
          R210_avg[ii] = (R2/R10);
        }
      }
      else if(ii == 39) {
        double R2 = SS_indR2[ii];
        double R10 = SS_indR10[ii];
        if(R10 > 0.0) {
          R210_avg[ii] = (R2/R10);
        }
      }
      else {
        double R2 = SS_indR2[ii] + SS_indR2[ii+1] + SS_indR2[ii+2];
        double R10 = SS_indR10[ii] + SS_indR10[ii+1] + SS_indR10[ii+2];
        if(R10 > 0.0) {
          R210_avg[ii] = (R2/R10);
        }
      }

    }
   

    bool foundReco = false;
    bool show = false;
    int new_SS = -1.0;
    
    if(E_gen_kin/E_beam > 0.4 || pdgID == -13) {

      for(int ii=0; ii < 40; ii++){ 
        double R210_i = ((SS_indR2[ii])/(SS_indR10[ii]));
        double R210_i_p1 = ((SS_indR2[ii+1])/(SS_indR10[ii+1]));
        double R210_i_p2 = ((SS_indR2[ii+2])/(SS_indR10[ii+2]));

        bool cond1 = (SS_inNrec_dR10[ii] > 3);

        bool cond1_1 = SS_inNrec_dR10[ii+1] > SS_inNrec_dR10[ii]; // Nrcht_10[ii+1] > Nrcht_10[ii]
        bool cond1_2 = SS_inNrec_dR10[ii] > SS_inNrec_dR10[ii-1]; // Nrcht_10[ii] > Nrcht_10[ii-1]
        bool cond2 = R210_i < 0.95; // R(2/10)[ii] < ratio_thres
        bool cond2_2 = R210_i_p1 <  R210_i;   // R(2/10)[ii+1] < R(2/10)[ii]
        bool cond2_3 = R210_i_p2 <  0.95;   // R(2/10)[ii+1] < R(2/10)[ii]
        bool cond2_4 = R210_i_p2 <  R210_i_p1;   // R(2/10)[ii+1] < R(2/10)[ii]
        bool cond2_5 = R210_i_p2 <=  R210_i;   // R(2/10)[ii+1] < R(2/10)[ii]

        bool cond3 = SS_indR10[ii+1] > SS_indR10[ii] ; // En(10)[ii+1] > En(10)[ii]
        bool cond3_3 = SS_indR10[ii] > SS_indR10[ii-1] ; // En(10)[ii] > En(10)[ii+1]

	bool cond4 = SS_indR10[ii] > ratio_thres ;   // En(2)[ii] > 15.0
        bool cond4_1 = SS_indR10[ii] > SS_indR10[ii-1] ;  // En(2)[ii] > En(2)[ii-1]
        bool cond4_2 = SS_indR10[ii] > SS_indR10[ii-2] ; // En(2)[ii] > En(2)[ii-2]
        
        bool cond4_3 = SS_indR10[ii+1] > 15.0;
        bool cond4_4 = SS_indR10[ii+2] > 15.0;
        
        bool cond5 = SS_indR10[ii] > 30.0;


        bool cond6 = SS_indR10_MW_fraction[ii] > 0.08;

        bool cond6_1 = SS_indR10_MW_fraction[ii] > 1.5*SS_indR10_MW_fraction[ii-1];
        bool cond6_2 = SS_indR10_MW_fraction[ii+1] > 1.2*SS_indR10_MW_fraction[ii];

        bool cond7 = (SS_indR10_MW_total[ii] + SS_indR10_MW_total[ii+1]) > 30.0;
	bool cond7_1 = (SS_indR10_MW_total[ii] + SS_indR10_MW_total[ii+1]) > ratio_thres;

        bool cond8 = SS_indR10_fraction[ii] > 1.5*SS_indR10_MW_fraction[ii-1];
        bool cond8_1 = SS_indR10_fraction[ii+1] > 1.5*SS_indR10_MW_fraction[ii];

	bool cond9 = R210_avg[ii] < 0.96;
        if((cond1 && cond9 && cond4) && !foundReco)  {
        
          new_SS = ii+1;
          foundReco = true;

        }
              
      }


      bool IsOut = false;
      if(closest_layer >= 1 && closest_layer <= 28) {
        if(abs(closest_layer - new_SS) > 3)  IsOut = true;
      }
      else {
        if(abs(closest_layer - new_SS) > 2)  IsOut = true; 
      }
      

      //if(show  && debug_count < 200 && IsOut) {
      if(show  && debug_count < 200) {
        char* temp = new char[100];
        sprintf(temp,"%0.2f",E_gen_kin);
        cout<<endl<<" **********   EVENT = "<<jentry<<" :: closest_layer = "<<closest_layer
            <<" :: Nsec = "<<nsec<<" :: Kin_sec = "<<temp<<"*****************"<<endl;
       
        cout<<endl<<"Layer NRhit(10cm)  R(2/5) En(2) En(5)\t\t R(2/10) En(2) En(10)\t\tFrac_1Layer(10cm)  Frac_3Layer(10cm) Total_3Layer(10cm)"<<endl<<endl;;

        for(int ii = 0; ii < 40; ii++) {

          char* toShow = new char[2000];
          sprintf(toShow,"%02d\t %d\t   %0.2f\t %0.2f\t %0.2f\t\t %0.2f\t %0.2f\t %0.2f\t\t %0.4f\t %0.4f\t %0.2f",ii+1,SS_inNrec_dR10[ii],((SS_indR2[ii])/(SS_indR5[ii])),SS_indR2[ii],SS_indR5[ii],((SS_indR2[ii])/(SS_indR10[ii])),SS_indR2[ii],SS_indR10[ii],SS_indR10_fraction[ii],SS_indR10_MW_fraction[ii],SS_indR10_MW_total[ii]);
          if(SS_indR10[ii] > 0.0 && SS_indR5[ii] > 0.0 && show) { 
            cout<<toShow<<endl;
          }
          else  {
            cout<<ii+1<<" : -1"<<endl;
          }
        }
        
       

        cout<<"Reconstructed layer = "<<new_SS<<endl;
        debug_count++;


      }
      
      if(IsOut) {
        double R210_i = ((SS_indR2[closest_layer-1])/(SS_indR10[closest_layer-1]));
        bool cond1 = (SS_inNrec_dR10[closest_layer-1] > 3);
        bool cond2 = R210_i < ratio_thres; // R(2/10)[ii] < ratio_thres
        bool cond4 = SS_indR10[closest_layer-1] > 15.0;
        int dec = ((int)cond1*1 + (int)cond2*2 + (int)cond4*4);
        if(closest_layer >= 1 && closest_layer <= 28) { 
          h_nsc_inv_EE->Fill(nsec); h_kin_sec_inv_EE->Fill(E_gen_kin); 
          h_cond_EE->Fill(dec);

        }
        else if(closest_layer >= 29) { 
          h_nsc_inv_FH->Fill(nsec); h_kin_sec_inv_FH->Fill(E_gen_kin); 
          h_cond_FH->Fill(dec);
        }
      }
      

      double R25_i = (SS_indR2[closest_layer-1]/SS_indR5[closest_layer-1]);
      double R25_i_m1 = (SS_indR2[closest_layer-2]/SS_indR5[closest_layer-2]);
      double R25_i_m2 = (SS_indR2[closest_layer-3]/SS_indR5[closest_layer-3]);
      double R25_i_p1 = (SS_indR2[closest_layer]/SS_indR5[closest_layer]);
      double R25_i_p2 = (SS_indR2[closest_layer+1]/SS_indR5[closest_layer+1]);

      double R210_i = (SS_indR2[closest_layer-1]/SS_indR10[closest_layer-1]);
      double R210_i_m1 = (SS_indR2[closest_layer-2]/SS_indR10[closest_layer-2]);
      double R210_i_m2 = (SS_indR2[closest_layer-3]/SS_indR10[closest_layer-3]);
      double R210_i_p1 = (SS_indR2[closest_layer]/SS_indR10[closest_layer]);
      double R210_i_p2 = (SS_indR2[closest_layer+1]/SS_indR10[closest_layer+1]);

      if(closest_layer >= 1 && closest_layer <= 40 ) {
         if(abs(closest_layer - new_SS) <= 1 && new_SS != -1) within1Layer_cog++;
         if(abs(closest_layer - new_SS) <= 2 && new_SS != -1) within2Layer_cog++;
         h_delta_layer_inc_cog->Fill((closest_layer - new_SS));
         total_hadInt_cog++;
         
      }
      if(closest_layer >= 1 && closest_layer <= 28) {
          h_EE_En2_layer_i->Fill(SS_indR2[closest_layer-1]); 
          h_EE_En2_layer_i_m1->Fill(SS_indR2[closest_layer-2]); 
          h_EE_En2_layer_i_m2->Fill(SS_indR2[closest_layer-3]);
          h_EE_En2_layer_i_p1->Fill(SS_indR2[closest_layer]); 
          h_EE_En2_layer_i_p2->Fill(SS_indR2[closest_layer+1]);

          h_EE_En5_layer_i->Fill(SS_indR5[closest_layer-1]); 
          h_EE_En5_layer_i_m1->Fill(SS_indR5[closest_layer-2]); 
          h_EE_En5_layer_i_m2->Fill(SS_indR5[closest_layer-3]);
          h_EE_En5_layer_i_p1->Fill(SS_indR5[closest_layer]); 
          h_EE_En5_layer_i_p2->Fill(SS_indR5[closest_layer+1]);

          h_EE_En10_layer_i->Fill(SS_indR10[closest_layer-1]); 
          h_EE_En10_layer_i_m1->Fill(SS_indR10[closest_layer-2]); 
          h_EE_En10_layer_i_m2->Fill(SS_indR10[closest_layer-3]);
          h_EE_En10_layer_i_p1->Fill(SS_indR10[closest_layer]); 
          h_EE_En10_layer_i_p2->Fill(SS_indR10[closest_layer+1]);


          h_EE_En10R_layer_i->Fill(SS_indR10[closest_layer-1]/SS_indR10[closest_layer-1]); 
          h_EE_En10R_layer_i_m1->Fill(SS_indR10[closest_layer-1]/SS_indR10[closest_layer-2]); 
          h_EE_En10R_layer_i_m2->Fill(SS_indR10[closest_layer-1]/SS_indR10[closest_layer-3]);
          h_EE_En10R_layer_i_p1->Fill(SS_indR10[closest_layer]/SS_indR10[closest_layer]); 
          h_EE_En10R_layer_i_p2->Fill(SS_indR10[closest_layer-1]/SS_indR10[closest_layer+1]);

          h_EE_Fr10R_layer_i->Fill(SS_indR10_MW_fraction[closest_layer-1]/SS_indR10_MW_fraction[closest_layer-1]); 
          h_EE_Fr10R_layer_i_m1->Fill(SS_indR10_MW_fraction[closest_layer-1]/SS_indR10_MW_fraction[closest_layer-2]); 
          h_EE_Fr10R_layer_i_m2->Fill(SS_indR10_MW_fraction[closest_layer-1]/SS_indR10_MW_fraction[closest_layer-3]);
          h_EE_Fr10R_layer_i_p1->Fill(SS_indR10_MW_fraction[closest_layer]/SS_indR10_MW_fraction[closest_layer-1]); 
          h_EE_Fr10R_layer_i_p2->Fill(SS_indR10_MW_fraction[closest_layer+1]/SS_indR10_MW_fraction[closest_layer-1]);

          h_EE_R25_layer_i->Fill(R25_i);
          h_EE_R25_layer_i_m1->Fill(R25_i_m1);
          h_EE_R25_layer_i_m2->Fill(R25_i_m2);
          h_EE_R25_layer_i_p1->Fill(R25_i_p1);
          h_EE_R25_layer_i_p2->Fill(R25_i_p2);

          h_EE_R210_layer_i->Fill(R210_i);
          h_EE_R210_layer_i_m1->Fill(R210_i_m1);
          h_EE_R210_layer_i_m2->Fill(R210_i_m2);
          h_EE_R210_layer_i_p1->Fill(R210_i_p1);
          h_EE_R210_layer_i_p2->Fill(R210_i_p2);

          h_EE_Nrechit10_layer_i->Fill(SS_inNrec_dR10[closest_layer-1]); 
          h_EE_Nrechit10_layer_i_m1->Fill(SS_inNrec_dR10[closest_layer-2]); 
          h_EE_Nrechit10_layer_i_m2->Fill(SS_inNrec_dR10[closest_layer-3]);
          h_EE_Nrechit10_layer_i_p1->Fill(SS_inNrec_dR10[closest_layer]); 
          h_EE_Nrechit10_layer_i_p2->Fill(SS_inNrec_dR10[closest_layer+1]);

          if(abs(closest_layer - new_SS) <= 1 && new_SS != -1) within1Layer_cog_EE++;
          if(abs(closest_layer - new_SS) <= 2 && new_SS != -1) within2Layer_cog_EE++;
          if(abs(closest_layer - new_SS) <= 3) within3Layer_cog_EE++;
          total_hadInt_cog_EE++;
          h_delta_layer_EE_cog->Fill((closest_layer - new_SS));
      }
      if(closest_layer >= 29 && closest_layer <= 40) {
        h_FH_En2_layer_i->Fill(SS_indR2[closest_layer-1]); 
        h_FH_En2_layer_i_m1->Fill(SS_indR2[closest_layer-2]); 
        h_FH_En2_layer_i_m2->Fill(SS_indR2[closest_layer-3]);
        h_FH_En2_layer_i_p1->Fill(SS_indR2[closest_layer]); 
        h_FH_En2_layer_i_p2->Fill(SS_indR2[closest_layer+1]);

        h_FH_En5_layer_i->Fill(SS_indR5[closest_layer-1]); 
        h_FH_En5_layer_i_m1->Fill(SS_indR5[closest_layer-2]); 
        h_FH_En5_layer_i_m2->Fill(SS_indR5[closest_layer-3]);
        h_FH_En5_layer_i_p1->Fill(SS_indR5[closest_layer]); 
        h_FH_En5_layer_i_p2->Fill(SS_indR5[closest_layer+1]);

        h_FH_En10_layer_i->Fill(SS_indR10[closest_layer-1]); 
        h_FH_En10_layer_i_m1->Fill(SS_indR10[closest_layer-2]); 
        h_FH_En10_layer_i_m2->Fill(SS_indR10[closest_layer-3]);
        h_FH_En10_layer_i_p1->Fill(SS_indR10[closest_layer]); 
        h_FH_En10_layer_i_p2->Fill(SS_indR10[closest_layer+1]);

        h_FH_En10R_layer_i_m1->Fill(SS_indR10[closest_layer-2]/SS_indR10[closest_layer-1]); 
        h_FH_En10R_layer_i_m2->Fill(SS_indR10[closest_layer-3]/SS_indR10[closest_layer-1]);
        h_FH_En10R_layer_i_p1->Fill(SS_indR10[closest_layer]/SS_indR10[closest_layer-1]); 
        h_FH_En10R_layer_i_p2->Fill(SS_indR10[closest_layer+1]/SS_indR10[closest_layer-1]);


        h_FH_Fr10R_layer_i->Fill(SS_indR10_MW_fraction[closest_layer-1]/SS_indR10_MW_fraction[closest_layer-1]); 
        h_FH_Fr10R_layer_i_m1->Fill(SS_indR10_MW_fraction[closest_layer-1]/SS_indR10_MW_fraction[closest_layer-2]); 
        h_FH_Fr10R_layer_i_m2->Fill(SS_indR10_MW_fraction[closest_layer-1]/SS_indR10_MW_fraction[closest_layer-3]);
        h_FH_Fr10R_layer_i_p1->Fill(SS_indR10_MW_fraction[closest_layer-1]/SS_indR10_MW_fraction[closest_layer]); 
        h_FH_Fr10R_layer_i_p2->Fill(SS_indR10_MW_fraction[closest_layer-1]/SS_indR10_MW_fraction[closest_layer+1]);

        h_FH_R25_layer_i->Fill(R25_i);
        h_FH_R25_layer_i_m1->Fill(R25_i_m1);
        h_FH_R25_layer_i_m2->Fill(R25_i_m2);
        h_FH_R25_layer_i_p1->Fill(R25_i_p1);
        h_FH_R25_layer_i_p2->Fill(R25_i_p2);

        h_FH_R210_layer_i->Fill(R210_i);
        h_FH_R210_layer_i_m1->Fill(R210_i_m1);
        h_FH_R210_layer_i_m2->Fill(R210_i_m2);
        h_FH_R210_layer_i_p1->Fill(R210_i_p1);
        h_FH_R210_layer_i_p2->Fill(R210_i_p2);

        h_FH_Nrechit10_layer_i->Fill(SS_inNrec_dR10[closest_layer-1]); 
        h_FH_Nrechit10_layer_i_m1->Fill(SS_inNrec_dR10[closest_layer-2]); 
        h_FH_Nrechit10_layer_i_m2->Fill(SS_inNrec_dR10[closest_layer-3]);
        h_FH_Nrechit10_layer_i_p1->Fill(SS_inNrec_dR10[closest_layer]); 
        h_FH_Nrechit10_layer_i_p2->Fill(SS_inNrec_dR10[closest_layer+1]);

        if(abs(closest_layer - new_SS) <= 1) within1Layer_cog_FH++;
        if(abs(closest_layer - new_SS) <= 2) within2Layer_cog_FH++;
        total_hadInt_cog_FH++;
        h_delta_layer_FH_cog->Fill((closest_layer - new_SS));
      }


      if((shower_start_index +1) >=1 && (shower_start_index +1) <= 28) TQ_SS_EE++;
      else if((shower_start_index +1) >=29 && (shower_start_index +1) <= 40) TQ_SS_FH++;
      else TQ_SS_notFound++;

      if(new_SS >= 1 && new_SS <= 28) SP_SS_EE++;
      else if(new_SS >= 29 && new_SS <= 40) SP_SS_FH++;
      else SP_SS_notFound++;

      if(closest_layer >= 1 && closest_layer <= 28) true_SS_EE++;
      else if(closest_layer >= 29 && closest_layer <= 40) true_SS_FH++;
      else true_SS_notFound++;
      
      total_SS++;
    }

    ////////////////////////////////////////////////////////////////////////////

    // if(E_gen_kin/E_beam > 0.4)  cout<<"jentry : closest_layer : new_SS :: "<<jentry<<" : "<<closest_layer<<" : "<<new_SS<<endl;

    // if(npi0 == 0) cout<<"jentry : new_SS :: "<<jentry<<" : "<<new_SS<<endl;
    ////////////////////////////////////////////////////////////////////
    //////////////// For shower shape debug  ///////////////////////////
    ////////////////////////////////////////////////////////////////////
    
    
    if(new_SS >= 3 && new_SS <= 7) {
      //cout<<"sim_Energy in AH : MIP_energy in AH "<<rechitEnergySumAHCAL<<endl;
      if(npi0 == 0) cout<<"jentry : npi0 :: "<<jentry<<" : "<<npi0<<endl;

      if(new_SS == 3) {
	h_npi0_SS01->Fill(npi0);
	h_all_pi0_energysum_SS01->Fill(pi0_Esum/trueBeamEnergy);
	h_lead_pi0_energy_SS01->Fill(getLeadingKE(KE_pi0)/trueBeamEnergy);
      }

      h_npi0_SS03_07->Fill(npi0);
      h_all_pi0_energysum_SS03_07->Fill(pi0_Esum/trueBeamEnergy);
      h_lead_pi0_energy_SS03_07->Fill(getLeadingKE(KE_pi0)/trueBeamEnergy);

      double pi0_frac = pi0_Esum/trueBeamEnergy;
      double w1 = 0.0106; // GeV/MIP
      double w2 = 0.0789; // GeV/MIP
      
      for(int iL = 0; iL < TOTAL_ACTIVE_LAYER; iL++){
	if(iL < EE_LAYER) {
	  if(pi0_frac < 0.2) {
	    if(new_SS == 1) {
	      h_longi_profile_SS01_NoPi0_raw->Fill(lambda[iL],RechitsEn_EE[iL]);
	      h_longi_profile_SS01_NoPi0_gev->Fill(lambda[iL],w1*RechitsEn_EE[iL]);
	    }
	    h_longi_profile_SS03_07_low_raw->Fill(lambda[iL],RechitsEn_EE[iL]);
	    h_longi_profile_SS03_07_low_gev->Fill(lambda[iL],w1*RechitsEn_EE[iL]);
	  }
	  else if(pi0_frac < 0.4) {
	    if(new_SS == 1) {
	      h_longi_profile_SS01_WithPi0_raw->Fill(lambda[iL],RechitsEn_EE[iL]);
	      h_longi_profile_SS01_WithPi0_gev->Fill(lambda[iL],w1*RechitsEn_EE[iL]);
	    }
	    h_longi_profile_SS03_07_mid_raw->Fill(lambda[iL],RechitsEn_EE[iL]);
	    h_longi_profile_SS03_07_mid_gev->Fill(lambda[iL],w1*RechitsEn_EE[iL]);
	    
	  }
	  else {
	    h_longi_profile_SS03_07_high_raw->Fill(lambda[iL],RechitsEn_EE[iL]);
	    h_longi_profile_SS03_07_high_gev->Fill(lambda[iL],w1*RechitsEn_EE[iL]);
	    
	  }
	  h_longi_profile_SS01_inclusive_raw->Fill(lambda[iL],RechitsEn_EE[iL]);
	  h_longi_profile_SS01_inclusive_gev->Fill(lambda[iL],w1*RechitsEn_EE[iL]);

	  h_longi_profile_SS03_07_inclusive_raw->Fill(lambda[iL],RechitsEn_EE[iL]);
	  h_longi_profile_SS03_07_inclusive_gev->Fill(lambda[iL],w1*RechitsEn_EE[iL]);

	  
	}
	else {
	  if(pi0_frac < 0.2) {
	    if(new_SS == 1) {
	      h_longi_profile_SS01_NoPi0_raw->Fill(lambda[iL],RechitsEn_FH[iL-EE_LAYER]);
	      h_longi_profile_SS01_NoPi0_gev->Fill(lambda[iL],w2*RechitsEn_FH[iL-EE_LAYER]);
	    }
	    h_longi_profile_SS03_07_low_raw->Fill(lambda[iL],RechitsEn_FH[iL-EE_LAYER]);
	    h_longi_profile_SS03_07_low_gev->Fill(lambda[iL],w2*RechitsEn_FH[iL-EE_LAYER]);

	  }
	  else if(pi0_frac < 0.4){
	    if(new_SS == 1) {
	      h_longi_profile_SS01_WithPi0_raw->Fill(lambda[iL],RechitsEn_FH[iL-EE_LAYER]);
	      h_longi_profile_SS01_WithPi0_gev->Fill(lambda[iL],w2*RechitsEn_FH[iL-EE_LAYER]);
	    }
	    h_longi_profile_SS03_07_mid_raw->Fill(lambda[iL],RechitsEn_FH[iL-EE_LAYER]);
	    h_longi_profile_SS03_07_mid_gev->Fill(lambda[iL],w2*RechitsEn_FH[iL-EE_LAYER]);

	  }
	  else {
	    h_longi_profile_SS03_07_high_raw->Fill(lambda[iL],RechitsEn_FH[iL-EE_LAYER]);
	    h_longi_profile_SS03_07_high_gev->Fill(lambda[iL],w2*RechitsEn_FH[iL-EE_LAYER]);

	  }
	  h_longi_profile_SS01_inclusive_raw->Fill(lambda[iL],RechitsEn_FH[iL-EE_LAYER]);
	  h_longi_profile_SS01_inclusive_gev->Fill(lambda[iL],w2*RechitsEn_FH[iL-EE_LAYER]);

	  h_longi_profile_SS03_07_inclusive_raw->Fill(lambda[iL],RechitsEn_FH[iL-EE_LAYER]);
	  h_longi_profile_SS03_07_inclusive_gev->Fill(lambda[iL],w2*RechitsEn_FH[iL-EE_LAYER]);


	}

      }
    }




    ////////////////////////////////////////////////////////////////////


    ////////////////////////////////////////////////////////////////////
    ////////////////        Summed up EE     ///////////////////////////
    ////////////////////////////////////////////////////////////////////
    
    
    if(new_SS >= 3 && new_SS <= 7) {
      double EE_01_07_low = 0.0;
      double EE_08_14_low = 0.0;
      double EE_15_21_low = 0.0;
      double EE_22_28_low = 0.0;

      double EE_01_07_mid = 0.0;
      double EE_08_14_mid = 0.0;
      double EE_15_21_mid = 0.0;
      double EE_22_28_mid = 0.0;

      double EE_01_07_high = 0.0;
      double EE_08_14_high = 0.0;
      double EE_15_21_high = 0.0;
      double EE_22_28_high = 0.0;

      double EE_01_07_inc = 0.0;
      double EE_08_14_inc = 0.0;
      double EE_15_21_inc = 0.0;
      double EE_22_28_inc = 0.0;

      double pi0_frac = pi0_Esum/trueBeamEnergy;
      double w1 = 0.0106; // GeV/MIP
      double w2 = 0.0789; // GeV/MIP
      
      // if(pi0_frac > 0.2) {
      // 	// cout<<"jentry : pi0_frac :: "<<jentry<<" : "<<pi0_frac<<endl;
      // 	// return;
      // 	continue;
      // }
      //cout<<"******* pi0_frac = "<<pi0_frac<<endl;
      for(int iL = 0; iL < TOTAL_ACTIVE_LAYER; iL++){
	if(iL < EE_LAYER) {
	  if(pi0_frac < 0.2) {
	    if(iL+1 <= 7) EE_01_07_low += RechitsEn_EE[iL];
	    else if(iL+1 <= 14) EE_08_14_low += RechitsEn_EE[iL];
	    else if(iL+1 <= 21) EE_15_21_low += RechitsEn_EE[iL];
	    else EE_22_28_low += RechitsEn_EE[iL];

	    switch(iL+1) {
	    case 7:
	      h_longi_profile_SS03_07_summed_EE_low_raw->Fill(lambda[iL],EE_01_07_low);
	      h_longi_profile_SS03_07_summed_EE_low_gev->Fill(lambda[iL],w1*EE_01_07_low);
	      break;
	    case 14:
	      h_longi_profile_SS03_07_summed_EE_low_raw->Fill(lambda[iL],EE_08_14_low);
	      h_longi_profile_SS03_07_summed_EE_low_gev->Fill(lambda[iL],w1*EE_08_14_low);
	      break;
	    case 21:
	      h_longi_profile_SS03_07_summed_EE_low_raw->Fill(lambda[iL],EE_15_21_low);
	      h_longi_profile_SS03_07_summed_EE_low_gev->Fill(lambda[iL],w1*EE_15_21_low);
	      break;
	    case 28:
	      h_longi_profile_SS03_07_summed_EE_low_raw->Fill(lambda[iL],EE_22_28_low);
	      h_longi_profile_SS03_07_summed_EE_low_gev->Fill(lambda[iL],w1*EE_22_28_low);
	      break;
	    default:
	      break;
	    }

	    
	  }
	  else if(pi0_frac < 0.4) {
	    if(iL+1 <= 7) EE_01_07_mid += RechitsEn_EE[iL];
	    else if(iL+1 <= 14) EE_08_14_mid += RechitsEn_EE[iL];
	    else if(iL+1 <= 21) EE_15_21_mid += RechitsEn_EE[iL];
	    else EE_22_28_mid += RechitsEn_EE[iL];

	    switch(iL+1) {
	    case 7:
	      h_longi_profile_SS03_07_summed_EE_mid_raw->Fill(lambda[iL],EE_01_07_mid);
	      h_longi_profile_SS03_07_summed_EE_mid_gev->Fill(lambda[iL],w1*EE_01_07_mid);
	      break;
	    case 14:
	      h_longi_profile_SS03_07_summed_EE_mid_raw->Fill(lambda[iL],EE_08_14_mid);
	      h_longi_profile_SS03_07_summed_EE_mid_gev->Fill(lambda[iL],w1*EE_08_14_mid);
	      break;
	    case 21:
	      h_longi_profile_SS03_07_summed_EE_mid_raw->Fill(lambda[iL],EE_15_21_mid);
	      h_longi_profile_SS03_07_summed_EE_mid_gev->Fill(lambda[iL],w1*EE_15_21_mid);
	      break;
	    case 28:
	      h_longi_profile_SS03_07_summed_EE_mid_raw->Fill(lambda[iL],EE_22_28_mid);
	      h_longi_profile_SS03_07_summed_EE_mid_gev->Fill(lambda[iL],w1*EE_22_28_mid);
	      break;
	    default:
	      break;
	    }

	  }
	  else {
	    if(iL+1 <= 7) EE_01_07_high += RechitsEn_EE[iL];
	    else if(iL+1 <= 14) EE_08_14_high += RechitsEn_EE[iL];
	    else if(iL+1 <= 21) EE_15_21_high += RechitsEn_EE[iL];
	    else EE_22_28_high += RechitsEn_EE[iL];

	    switch(iL+1) {
	    case 7:
	      h_longi_profile_SS03_07_summed_EE_high_raw->Fill(lambda[iL],EE_01_07_high);
	      h_longi_profile_SS03_07_summed_EE_high_gev->Fill(lambda[iL],w1*EE_01_07_high);
	      break;
	    case 14:
	      h_longi_profile_SS03_07_summed_EE_high_raw->Fill(lambda[iL],EE_08_14_high);
	      h_longi_profile_SS03_07_summed_EE_high_gev->Fill(lambda[iL],w1*EE_08_14_high);
	      break;
	    case 21:
	      h_longi_profile_SS03_07_summed_EE_high_raw->Fill(lambda[iL],EE_15_21_high);
	      h_longi_profile_SS03_07_summed_EE_high_gev->Fill(lambda[iL],w1*EE_15_21_high);
	      break;
	    case 28:
	      h_longi_profile_SS03_07_summed_EE_high_raw->Fill(lambda[iL],EE_22_28_high);
	      h_longi_profile_SS03_07_summed_EE_high_gev->Fill(lambda[iL],w1*EE_22_28_high);
	      break;
	    default:
	      break;
	    }

	  }

	  if(iL+1 <= 7) EE_01_07_inc += RechitsEn_EE[iL];
	  else if(iL+1 <= 14) EE_08_14_inc += RechitsEn_EE[iL];
	  else if(iL+1 <= 21) EE_15_21_inc += RechitsEn_EE[iL];
	  else EE_22_28_inc += RechitsEn_EE[iL];


	  switch(iL+1) {
	  case 7:
	    h_longi_profile_SS03_07_summed_EE_inclusive_raw->Fill(lambda[iL],EE_01_07_inc);
	    h_longi_profile_SS03_07_summed_EE_inclusive_gev->Fill(lambda[iL],w1*EE_01_07_inc);
	    break;
	  case 14:
	    h_longi_profile_SS03_07_summed_EE_inclusive_raw->Fill(lambda[iL],EE_08_14_inc);
	    h_longi_profile_SS03_07_summed_EE_inclusive_gev->Fill(lambda[iL],w1*EE_08_14_inc);
	    break;
	  case 21:
	    h_longi_profile_SS03_07_summed_EE_inclusive_raw->Fill(lambda[iL],EE_15_21_inc);
	    h_longi_profile_SS03_07_summed_EE_inclusive_gev->Fill(lambda[iL],w1*EE_15_21_inc);
	    break;
	  case 28:
	    h_longi_profile_SS03_07_summed_EE_inclusive_raw->Fill(lambda[iL],EE_22_28_inc);
	    h_longi_profile_SS03_07_summed_EE_inclusive_gev->Fill(lambda[iL],w1*EE_22_28_inc);
	    break;
	  default:
	    break;
	  }
	  // if(iL+1 <= 7)
	  //   cout<<"iL+1 : Energy Layer : EE_01_07_low : EE_01_07_mid : EE_01_07_high :: "<<iL+1<<" : "<<RechitsEn_EE[iL]<<" : "<<EE_01_07_low<<" : "<<EE_01_07_mid<<" : "<<EE_01_07_high<<endl;
	  // else return;
	}
	else {
	  if(pi0_frac < 0.2) {
	    h_longi_profile_SS03_07_summed_EE_low_raw->Fill(lambda[iL],RechitsEn_FH[iL-EE_LAYER]);
	    h_longi_profile_SS03_07_summed_EE_low_gev->Fill(lambda[iL],w2*RechitsEn_FH[iL-EE_LAYER]);

	  }
	  else if(pi0_frac < 0.4){
	    h_longi_profile_SS03_07_summed_EE_mid_raw->Fill(lambda[iL],RechitsEn_FH[iL-EE_LAYER]);
	    h_longi_profile_SS03_07_summed_EE_mid_gev->Fill(lambda[iL],w2*RechitsEn_FH[iL-EE_LAYER]);

	  }
	  else {
	    h_longi_profile_SS03_07_summed_EE_high_raw->Fill(lambda[iL],RechitsEn_FH[iL-EE_LAYER]);
	    h_longi_profile_SS03_07_summed_EE_high_gev->Fill(lambda[iL],w2*RechitsEn_FH[iL-EE_LAYER]);

	  }
	  h_longi_profile_SS01_inclusive_raw->Fill(lambda[iL],RechitsEn_FH[iL-EE_LAYER]);
	  h_longi_profile_SS01_inclusive_gev->Fill(lambda[iL],w2*RechitsEn_FH[iL-EE_LAYER]);

	  h_longi_profile_SS03_07_summed_EE_inclusive_raw->Fill(lambda[iL],RechitsEn_FH[iL-EE_LAYER]);
	  h_longi_profile_SS03_07_summed_EE_inclusive_gev->Fill(lambda[iL],w2*RechitsEn_FH[iL-EE_LAYER]);


	}

      }
    }

    /////////////////////////////////////////////////////
    ////////////////////////////////////////////////////
    
    if(Nrechit_EE == 0) {
      h_hadInteraction_withzeroEEhits->Fill(HadInt_z/100.0);
      ee_zeroHit++;
    }
    else {
      h_hadInteraction_withNonzeroEEhits->Fill(HadInt_z/100.0);
    }

    //if(jentry > 5) return;
    if(DEBUG) cout<<"DEBUG: NRechits in EE = "<<Nrechit_EE<<"\t Rechit_energy in EE = "<<rechitEnergySum_EE<<endl;
    if(DEBUG) cout<<"DEBUG: NRechits in FH = "<<Nrechit_FH<<"\t Rechit_energy in FH = "<<rechitEnergySum_FH<<endl;

    double isMuonLike = false;
    bool isRegion1 = false;
    bool isRegion2 = false;
    bool isRegion3 = false;




    ///////////////////////////////////////////////////////////
    //// Original Event Classification                       //
    ///////////////////////////////////////////////////////////

    // if(Nrechit_EE < 80 && Nrechit_FH < 60) {  //MIP LIKE
    //   isRegion1 = true;
    //   isRegion2 = false;
    //   isRegion3 = false;

    // }
    // else if(Nrechit_EE < 80) { // H-Hadron
    //   //classify as H-hadron
    //   isRegion2 = true;
    //   isRegion1 = false;
    //   isRegion3 = false;
    // }
    // else {  // EH-Hadron
    //   isRegion3 = true;
    //   isRegion1 = false;
    //   isRegion2 = false;
	
    // }



    ///////////////////////////////////////////////////////////
    //// WITH ENERGY RATIOS  (IMPROVED Event Classification) //
    ///////////////////////////////////////////////////////////

    // if(Nrechit_EE < 80 && Nrechit_FH < 60 && shower_start_index == -1) {  //MIP LIKE
    //   isRegion1 = true;
    //   isRegion2 = false;
    //   isRegion3 = false;
    //   h_shower_start_reg1->Fill(shower_lambda_);
    //   region_1_classified_events++;
    // }
    // else if(Nrechit_EE < 80 && shower_start_index > 27) { //H hadrons
    //   isRegion1 = false;
    //   isRegion2 = true;
    //   isRegion3 = false;
    //   h_shower_start_reg2->Fill(shower_lambda_);
    //   region_2_classified_events++;
    // }
    // else if(Nrechit_EE > 80 && shower_start_index < 28  && shower_start_index > -1){  //EH hadron
    //   isRegion1 = false;
    //   isRegion2 = false;
    //   isRegion3 = true;
    //   h_shower_start_reg3->Fill(shower_lambda_);
    //   region_3_classified_events++;
    // }
    // else {
    //   h_shower_start_CHECK_FOR_UNCLASSIFIED->Fill(shower_lambda_);
    //   h_Nrechit_EE_vs_FH_CHECK_FOR_UNCLASSIFIED->Fill(Nrechit_EE,Nrechit_FH);
    //   if(shower_start_index == -1)
    // 	h_Nrechit_EE_vs_FH_CHECK_FOR_UNCLASSIFIED_R1->Fill(Nrechit_EE,Nrechit_FH);
    //   else if (shower_start_index < 28  && shower_start_index > -1)
    // 	h_Nrechit_EE_vs_FH_CHECK_FOR_UNCLASSIFIED_R3->Fill(Nrechit_EE,Nrechit_FH);
    //   else
    // 	h_Nrechit_EE_vs_FH_CHECK_FOR_UNCLASSIFIED_R2->Fill(Nrechit_EE,Nrechit_FH);
	     
    //   non_classified_events++;
    // }

    /////////////////////////////////////////////////////////////


    /////////////////////////////////////////////////////////////
    //// IMPROVED Event Classification only with shower start //
    ////////////////////////////////////////////////////////////

    if(shower_start_index == -1) {  //MIP LIKE
      isRegion1 = true;
      isRegion2 = false;
      isRegion3 = false;
      h_shower_start_reg1->Fill(shower_lambda_);
      region_1_classified_events++;
    }
    else if(shower_start_index > 27) { //H hadrons
      isRegion1 = false;
      isRegion2 = true;
      isRegion3 = false;
      h_shower_start_reg2->Fill(shower_lambda_);
      region_2_classified_events++;
    }
    else {  //EH hadron
      isRegion1 = false;
      isRegion2 = false;
      isRegion3 = true;
      h_shower_start_reg3->Fill(shower_lambda_);
      region_3_classified_events++;
    }


    // else {  //EH hadron
    //   isRegion1 = false;
    //   isRegion2 = false;
    //   isRegion3 = true;
    // }
    




    h_Nrechit_EE->Fill(Nrechit_EE);
    h_Nrechit_FH->Fill(Nrechit_FH);
    // h_Nrechit_AH->Fill(ahc_nHits);
    //rechitEnergySum_EE
    h_rechit_energy_raw_EE->Fill(rechitEnergySum_EE);    
    h_rechit_energy_raw_FH->Fill(rechitEnergySum_FH);
    h_rechit_energy_raw_all->Fill(rechitEnergySum_EE+rechitEnergySum_FH);    

    h_rechit_energy_raw_low_EE->Fill(rechitEnergySum_EE);    
    h_rechit_energy_raw_low_FH->Fill(rechitEnergySum_FH);
    h_rechit_energy_raw_low_all->Fill(rechitEnergySum_EE+rechitEnergySum_FH);    


    h_rechit_energy_raw_mid_EE->Fill(rechitEnergySum_EE);    
    h_rechit_energy_raw_mid_FH->Fill(rechitEnergySum_FH);
    h_rechit_energy_raw_mid_all->Fill(rechitEnergySum_EE+rechitEnergySum_FH);    

    h_rechit_energy_raw_high_EE->Fill(rechitEnergySum_EE);    
    h_rechit_energy_raw_high_FH->Fill(rechitEnergySum_FH);
    h_rechit_energy_raw_high_all->Fill(rechitEnergySum_EE+rechitEnergySum_FH);    




    h_Nrechit_EE_vs_FH->Fill(Nrechit_EE,Nrechit_FH);
    ///////////////////////////////////////////////////////////
    //////////////////////////////////////////////////////////
    // h_rechit_energy_FH_AHCAL_FB->Fill(rechitEnergySum_FH,rechitEnergySumAHCAL);
    // h_rechit_energy_EE_FH->Fill(rechitEnergySum_EE,rechitEnergySum_FH);

    if(isRegion2) {
    // if(!isRegion3) {
      // if(ahc_nHits > 35 && ahc_nHits < 45 && rechitEnergySumAHCAL < 250) {
      // 	cout<<"Do NOT caluculate"<<endl;
      // 	continue;
      // }
      // h_Nrechit_EE->Fill(Nrechit_EE);
      // h_Nrechit_FH->Fill(Nrechit_FH);
      // h_Nrechit_AH->Fill(ahc_nHits);

      // h_Nrechit_EEFH_vs_AH->Fill(Nrechit_EE+Nrechit_FH,ahc_nHits);
      // double a = 0.0351296;
      // double b = 0.031534;
      // h_Nrechits_EE_FH_FB->Fill(Nrechit_EE,Nrechit_FH);
      h_rechit_energy_inEE_region2_mips->Fill(rechitEnergySum_EE);
      double a = -1.0;
      double b = -1.0;
      // double a = 0.0611934;
      // double b = 0.0183451;
      double alpha = -1.0;  // FH & AHCAL relative weight
      double gamma = -1.0;
      //double MIP_to_GeV_conversion_factor = 12.83;
      
      if(beamEnergy == 20) {
      	// a = 0.0611934; b = 0.0183451;
      	//a = 0.0630539 ;b = 0.0120932;
      	a = 0.0604318 ;b = 0.0307894;
      	//alpha = 0.45; gamma = 13.205;
      	alpha = 0.402; gamma = 12.826;
      }
      if(beamEnergy == 50) {
      	// a = 0.0608286 ;b = 0.0136046;
      	//a = 0.076075 ;b = 0.0116685;
      	a = 0.0679195 ;b = 0.0324864;
      	//alpha = 0.40; gamma = 12.856;
      	alpha = 0.404; gamma = 12.811;
      }
      if(beamEnergy == 80) {
      	// a = 0.0622612 ;b = 0.0152219;
      	//a = 0.0788491 ;b = 0.013591;
      	a = 0.0683456 ;b = 0.0320886;
      	//alpha = 0.40; gamma = 13.1375;
      	alpha = 0.405; gamma = 13.225;
      }
      if(beamEnergy == 100) {
      	//a = 0.0786438 ;b = 0.0151615;
      	a = 0.0677498 ;b = 0.031738;
      	//alpha = 0.40; gamma = 13.38;
      	alpha = 0.389; gamma = 13.323;
      }
      if(beamEnergy == 120) {
      	// a = 0.0600868 ;b = 0.018315;
      	//a = 0.0794755 ;b = 0.0174122;
      	a = 0.068558 ;b = 0.0314515;
      	//alpha = 0.40; gamma = 13.38;
      	alpha = 0.408; gamma = 13.556;
      }
      if(beamEnergy == 200) {
      	//a = 0.0802285 ;b = 0.0178579;
      	//a = 0.0681734 ;b = 0.031085;
      	a = 0.0678221 ;b = 0.0308716;
      	//alpha = 0.40; gamma = 13.635;
      	alpha = 0.422; gamma = 13.941;
      }
      if(beamEnergy == 250) {
      	// a = 0.0802285 ;b = 0.0178579;
      	// a = 0.0804709 ;b = 0.0182122;
      	a = 0.0678221 ;b = 0.0308716;
      	//alpha = 0.40; gamma = 13.756;
      	alpha = 0.432; gamma = 14.196;
      }
      if(beamEnergy == 300) {
      	//a = 0.0844364 ;b = 0.0148193;
      	a = 0.0703497 ;b = 0.0293021;
      	//alpha = 0.40; gamma = 13.7266;
      	alpha = 0.434; gamma = 14.193;
      }
      
      // double tot_ = (rechitEnergySum+rechitEnergySumAHCAL);

      //For 50 GeV
      // alpha = 0.40; gamma = 12.856;
      double tot_ = beamEnergy;
      // int AH_frac = ((a*rechitEnergySumAHCAL)/tot_)*100;
      // int FH_frac = ((b*rechitEnergySum_FH)/tot_)*100;
      // int AH_frac = ((a*rechitEnergySumAHCAL));
      // int FH_frac = ((b*rechitEnergySum_FH));

      //double sigma2_FB = (rechitEnergySumAHCAL+rechitEnergySum_FH);
      // double sigma2_FB = 1.0;

      h_shower_start_reg2->Fill(shower_lambda_);

      // h_rechit_energy_FB_FH_alone_raw->Fill(rechitEnergySum_FH);

      //double MIP_to_GeV_conversion_factor = 11.214;
      // double MIP_to_GeV_conversion_factor = 9.51;
      // double MIP_to_GeV_conversion_factor = 13.7;

      double MIP_to_GeV_conversion_factor = gamma;
      double FH_Gev_Scale = (rechitEnergySum_FH/MIP_to_GeV_conversion_factor);
      //double AH_Gev_Scale = (rechitEnergySumAHCAL/MIP_to_GeV_conversion_factor);
      double AH_Gev_Scale = rechitEnergySumAHCAL;
      //double full_energy = FH_Gev_Scale+(alpha*AH_Gev_Scale);
      double full_energy = FH_Gev_Scale+AH_Gev_Scale;

      

      // h_rechit_energy_FH_FB_weighted->Fill(FH_Gev_Scale);
      // h_rechit_energy_FH_AHCAL_FB_raw->Fill(rechitEnergySum_FH, rechitEnergySumAHCAL);
      // h_rechit_energy_FH_AHCAL_FB_weighted->Fill(FH_Gev_Scale,(alpha*AH_Gev_Scale));

      h_energy_region2_gev->Fill(full_energy);
      h_energy_region2_mips->Fill(rechitEnergySum_FH+alpha*rechitEnergySumAHCAL);
      h_energy_region2_gev_withoutAHCAL->Fill(FH_Gev_Scale);
      h_energy_all_gev->Fill(full_energy);
      h_energy_FH_vs_AHCAL_region2_gev->Fill(FH_Gev_Scale,AH_Gev_Scale);

      h_sim_energy_inGev_region2->Fill(energyLostFH+energyLostBH);
      h_sim_energy_all_gev->Fill(energyLostFH+energyLostBH);



      for(int iL = 0; iL < TOTAL_ACTIVE_LAYER; iL++){
  	    if(iL < EE_LAYER) {
  	      // Nrechit_EE+=NRechits_EE[iL];
  	      rechitEnergySum_EE+=RechitsEn_EE[iL];
  	      // cogX_[iL] = cogX_[iL]/RechitsEn_EE[iL];
  	      // cogY_[iL] = cogY_[iL]/RechitsEn_EE[iL];
  	    }
  	    else {
  	      Nrechit_FH+=NRechits_FH[iL-EE_LAYER];
  	      rechitEnergySum_FH+=RechitsEn_FH[iL-EE_LAYER];
  	      cogX_[iL] = cogX_[iL]/RechitsEn_FH[iL-EE_LAYER];
  	      cogY_[iL] = cogY_[iL]/RechitsEn_FH[iL-EE_LAYER];
  	      
  	    }

  	  }
	  }
    
    

    //if(isRegion3 && shower_start_index>0) {
    if(isRegion3) {
      // double MIP_to_GeV = 57.6;
      // double MIP_to_GeV = 60.56;
      // double alpha = 0.4;
      // // double beta = 4.4;
      // double beta = 4.9;

      double alpha = -1.0; // EE & FH + AHCAL relative weight
      double beta = -1.0; // FH & AHCAL relative weight
      double gamma = -1.0; // FH & AHCAL relative weight
      
      if(beamEnergy == 20) {
      	//beta = 4.8; alpha = 0.45; gamma = 54.15;
      	beta = 4.207; alpha = 0.402; gamma = 51.223;
      }
      if(beamEnergy == 50) {
      	//beta = 5.3; alpha = 0.40; gamma = 62.92;
      	beta = 4.752; alpha = 0.404; gamma = 59.03;
      }
      if(beamEnergy == 80) {
      	//beta = 5.0; alpha = 0.40; gamma = 63.737;
      	beta = 5.009; alpha = 0.405; gamma = 63.354;
      }
      if(beamEnergy == 100) {
      	//beta = 4.9; alpha = 0.40; gamma = 64.41;
      	beta = 5.174; alpha = 0.398; gamma = 65.967;
      }
      if(beamEnergy == 120) {
      	//beta = 5.0; alpha = 0.40; gamma = 65.812;
      	beta = 5.074; alpha = 0.408; gamma = 65.947;
      }
      if(beamEnergy == 200) {
      	//beta = 4.5; alpha = 0.40; gamma = 63.3;
      	beta = 5.094; alpha = 0.422; gamma = 68.321;
      }
      if(beamEnergy == 250) {
      	//beta = 5.6; alpha = 0.40; gamma = 73.76;
      	beta = 5.305; alpha = 0.432; gamma = 71.396;
      }
      if(beamEnergy == 300) {
      	//beta = 5.5; alpha = 0.40; gamma = 73.567;
      	beta = 5.253; alpha = 0.434; gamma = 71.443;
      }

      //For 50 GeV
      // beta = 5.3; alpha = 0.40; gamma = 62.92;
      double MIP_to_GeV = gamma;
      // double tot_E_MIPs = rechitEnergySum_EE + beta*(rechitEnergySum_FH+(alpha*rechitEnergySumAHCAL));
      double EE_FH_MIPs = rechitEnergySum_EE + beta*(rechitEnergySum_FH);
      //double tot_E_gev = tot_E_MIPs/MIP_to_GeV;
      double EE_FH_gev = EE_FH_MIPs/MIP_to_GeV;
      double EE_gev = rechitEnergySum_EE/MIP_to_GeV;
      double FH_gev = (beta*rechitEnergySum_FH)/MIP_to_GeV;
      double AH_gev = rechitEnergySumAHCAL ;
      double tot_E_gev = EE_gev + FH_gev + AH_gev;
      
      // h_rechit_energy_raw_region3->Fill(rechitEnergySum_EE+rechitEnergySum_FH+rechitEnergySumAHCAL);
      h_rechit_energy_raw_EE_alone_region3->Fill(rechitEnergySum_EE);
      h_rechit_energy_raw_FH_AH_alone_region3->Fill(beta*rechitEnergySum_FH);
      
      // h_rechit_energy_part_weighted_region3->Fill(rechitEnergySum_EE+rechitEnergySum_FH+(alpha*rechitEnergySumAHCAL));
      //For sim-data comparison
      h_rechit_energy_part_weighted_region3->Fill(rechitEnergySum_EE+rechitEnergySum_FH);

      // h_rechit_energy_full_weighted_region3->Fill(rechitEnergySum_EE+beta*(rechitEnergySum_FH+(alpha*rechitEnergySumAHCAL)));
      h_Nrechit_EE_vs_FH_region3->Fill(Nrechit_EE,Nrechit_FH);

      h_rechit_energy_full_weighted_region3->Fill(rechitEnergySum_EE+ beta*rechitEnergySum_FH);

      h_rechit_energy_full_weighted_inGev_region3->Fill(tot_E_gev);
      h_rechit_energy_full_weighted_inGev_region3_withoutAHCAL->Fill(EE_FH_gev);
      // h_rechit_energy_EE_FHAHCAL_EH_mips->Fill(rechitEnergySum_EE,(beta*(rechitEnergySum_FH+(alpha*rechitEnergySumAHCAL))));
      //h_rechit_energy_EE_FHAHCAL_EH_gev->Fill(rechitEnergySum_EE/MIP_to_GeV,(beta*(rechitEnergySum_FH+(alpha*rechitEnergySumAHCAL)))/MIP_to_GeV);
      h_rechit_energy_EE_FHAHCAL_EH_gev->Fill(EE_gev,(FH_gev+AH_gev));

      if(tot_E_gev < 100) {
      	h_rechit_energy_EE_FHAHCAL_EH_gev_DEBUG->Fill(EE_gev,(FH_gev+AH_gev));
      	h_shower_start_reg3_DEBUG->Fill(shower_lambda_);
      }

      h_energy_all_gev->Fill(tot_E_gev);
      h_sim_energy_inGev_region3->Fill(energyLostEE+energyLostFH+energyLostBH);
      h_sim_energy_all_gev->Fill(energyLostEE+energyLostFH+energyLostBH);


      h_shower_start_reg3->Fill(shower_lambda_);
      // if(!strcmp(data,"data")) {
      // 	float rel_weight = 0.0;
      // 	//float start = 4.0;
      // 	float start = 3.0;
      // 	float step = 0.1;
      // 	for(int ii = 0; ii < 50; ii++) {
      // 	  rel_weight = start+(step*ii);
      // 	  h_rechit_energy_EH_rel_weightScan[ii]->Fill(rechitEnergySum_EE+rel_weight*(rechitEnergySum_FH+(alpha*rechitEnergySumAHCAL)));
      // 	}
      // }
      

    }

    if(shower_lambda_ > -1 || true) {
      //h_total_sim_energy->Fill(energyLostEE+energyLostFH+energyLostBH+energyLostBeam+energyLostOutside);
      h_total_sim_energy->Fill(energyLostEE+energyLostFH+energyLostBH+energyLostBeam);
    }
    if(DEBUG) cout<<"DEBUG: End of Event = "<<jentry+1<<endl;
    if(DEBUG) cout<<"DEBUG: ****************** "<<endl;
    //if(DEBUG && jentry > 10) return;
    if(DEBUG) return;

    // if(jentry > 2) break;
    // if(shower_start_index < 0)
    //   cout<<jentry+1<<endl;
    // if(region_1_classified_events > 20) break;
  } // loop over entries

  cout<<"Got Out "<<jentry<<endl;
  Long64_t total_events = (region_1_classified_events+region_2_classified_events+region_3_classified_events+non_classified_events);
  
  // cout<<"Events with zero hits in AHCAL = "<<ahc_zeroHit<<endl;
  // cout<<"Events with zero hits in CE-E = "<<ee_zeroHit<<endl;
  // cout<<"MIP like events = "<<((float)region_1_classified_events*100.0)/total_events<<"%"<<endl;
  // cout<<"shower start in EE = "<<((float)region_3_classified_events*100.0)/total_events<<"%"<<endl;
  // cout<<"shower start in FH = "<<((float)region_2_classified_events*100.0)/total_events<<"%"<<endl;
  // cout<<"Non-classified events = "<<((float)non_classified_events*100.0)/total_events<<"%"<<endl;
  // //cout<<"Sum = "<<(region_1_classified_events+region_2_classified_events+region_3_classified_events+non_classified_events)<<endl;           
  // cout<<"Sum = "<<total_events<<endl;

    char* temp_EE = new char[500];
    char* temp_FH = new char[500];
    char* temp_all = new char[500];


    cout<<endl<<" Event classification "<<endl;
    sprintf(temp_EE,"TQ showering in EE : %0.1f%%",(TQ_SS_EE*100.0)/total_SS);
    sprintf(temp_FH,"TQ MIPs in EE : %0.1f%%",(TQ_SS_FH*100.0)/total_SS);
    sprintf(temp_all,"TQ Not found : %0.1f%%",(TQ_SS_notFound*100.0)/total_SS);

    cout<<temp_EE<<endl;
    cout<<temp_FH<<endl;
    cout<<temp_all<<endl;

    sprintf(temp_EE,"SP showering in EE : %0.1f%%",(SP_SS_EE*100.0)/total_SS);
    sprintf(temp_FH,"SP MIPs in EE : %0.1f%%",(SP_SS_FH*100.0)/total_SS);
    sprintf(temp_all,"SP Not found : %0.1f%%",(SP_SS_notFound*100.0)/total_SS);


    cout<<temp_EE<<endl;
    cout<<temp_FH<<endl;
    cout<<temp_all<<endl;


    sprintf(temp_EE,"true showering in EE : %0.1f%%",(true_SS_EE*100.0)/total_SS);
    sprintf(temp_FH,"true MIPs in EE : %0.1f%%",(true_SS_FH*100.0)/total_SS);
    sprintf(temp_all,"true Not found : %0.1f%%",(true_SS_notFound*100.0)/total_SS);


    cout<<temp_EE<<endl;
    cout<<temp_FH<<endl;
    cout<<temp_all<<endl;

    cout<<endl<<" Effiecieny "<<endl;
    sprintf(temp_EE,"CE-E effiecieny between  1 <= layer <= 28 within 1 : within 2 : within 3 :: %0.1f%%: %0.1f%%: %0.1f%%",(within1Layer_EE*100.0)/total_hadInt_EE,(within2Layer_EE*100.0)/total_hadInt_EE, (within3Layer_EE*100.0)/total_hadInt_EE);

    sprintf(temp_FH,"CE-E effiecieny between  29 <= layer <= 40 within 1 : within 2 :: %0.1f%% : %0.1f%%",(within1Layer_FH*100.0)/total_hadInt_FH,(within2Layer_FH*100.0)/total_hadInt_FH);

    sprintf(temp_all,"Overall effiecieny between  1 <= layer <= 40 within 1 : within 2 :: %0.1f%% : %0.1f%%",(within1Layer*100.0)/total_hadInt,(within2Layer*100.0)/total_hadInt);

    cout<<temp_EE<<endl;
    cout<<temp_FH<<endl;
    cout<<temp_all<<endl;


  // cout<<"EE effiecieny between  3 <= layer <= 28 witin +/- 1 : witin +/- 2 :: "<<(within1Layer_EE*100.0)/total_hadInt_EE
  //   <<"%% : "<<(within2Layer_EE*100.0)/total_hadInt_EE<<"%%"<<endl;
  
  // cout<<"FH effiecieny between  32 <= layer <= 40 witin +/- 1 : witin +/- 2 :: "<<(within1Layer_FH*100.0)/total_hadInt_FH
  //   <<"%% : "<<(within2Layer_FH*100.0)/total_hadInt_FH<<"%%"<<endl;

  // cout<<"Overall effiecieny between  3 <= layer <= 40 witin +/- 1 : witin +/- 2 :: "<<(within1Layer*100.0)/total_hadInt
  //   <<"%% : "<<(within2Layer*100.0)/total_hadInt<<"%%"<<endl;

  

  // cout<<"Total true events =  "<<hadronic_int<<endl;
  // cout<<"Hadronic interaction not found =  "<<notFound_hadronic_int<<endl;
  // cout<<"Hard interaction = "<<hard_hadronic_int<<endl;



  // cout<<"***CHECK***"<<endl;
  // cout<<"(within2Layer) EE : in : FH : Overall :: "<<within2Layer_EE<<" : "<<within2Layer_interface<< " : "<<within2Layer_FH<< " : "<<within2Layer<<endl;
  // cout<<"(Total) EE : in : FH : Overall :: "<<total_hadInt_EE<<" : "<<total_hadInt_interface<< " : "<<total_hadInt_FH<< " : "<<total_hadInt<<endl;
  // cout<<"***********"<<endl;
    cout<<"Energy Threshold : "<<ratio_thres<<endl;

    sprintf(temp_EE,"CE-E effiecieny between  1 <= layer <= 28 within 1 : within 2 : within 3 :: %0.1f%%: %0.1f%%: %0.1f%%",(within1Layer_cog_EE*100.0)/total_hadInt_cog_EE,(within2Layer_cog_EE*100.0)/total_hadInt_cog_EE, (within3Layer_cog_EE*100.0)/total_hadInt_cog_EE);

    sprintf(temp_FH,"CE-E effiecieny between  29 <= layer <= 40 within 1 : within 2 :: %0.1f%% : %0.1f%%",(within1Layer_cog_FH*100.0)/total_hadInt_cog_FH,(within2Layer_cog_FH*100.0)/total_hadInt_cog_FH);

    sprintf(temp_all,"Overall effiecieny between  1 <= layer <= 40 within 1 : within 2 :: %0.1f%% : %0.1f%%",(within1Layer_cog*100.0)/total_hadInt_cog,(within2Layer_cog*100.0)/total_hadInt_cog);

  cout<<"***New SS***"<<endl;
  cout<<temp_EE<<endl;
  cout<<temp_FH<<endl;
  cout<<temp_all<<endl;
  
  cout<<within1Layer_cog_EE<<":"<<within2Layer_cog_EE<<":"<<within3Layer_cog_EE<<":"<<total_hadInt_cog_EE<<endl;
  cout<<within1Layer_cog_FH<<":"<<within2Layer_cog_FH<<":"<<total_hadInt_cog_FH<<endl;
  cout<<within1Layer_cog<<":"<<within2Layer_cog<<":"<<total_hadInt_cog<<endl;
  // cout<<"CE-E effiecieny between  3 <= layer <= 26 witin +/- 1 : witin +/- 2 :: "<<(within1Layer_cog_EE*100.0)/total_hadInt_cog_EE<<"%% : "<<(within2Layer_cog_EE*100.0)/total_hadInt_cog_EE<<"%%"<<endl;
  // cout<<"CE-H effiecieny between  32 <= layer <= 38 witin +/- 1 : witin +/- 2 :: "<<(within1Layer_cog_FH*100.0)/total_hadInt_cog_FH<<"%% : "<<(within2Layer_cog_FH*100.0)/total_hadInt_cog_FH<<"%%"<<endl;
  // cout<<"Overall effiecieny between  3 <= layer <= 38 witin +/- 1 : witin +/- 2 :: "<<(within1Layer_cog*100.0)/total_hadInt_cog<<"%% : "<<(within2Layer_cog*100.0)/total_hadInt_cog<<"%%"<<endl;
  // cout<<ratio_thres<<" <= Thres; Overall effiecieny between  1 < layer < 29 witin +/- 1 : witin +/- 2 :: "<<(within1Layer_cog)<<" : "<<(within2Layer_cog)<<" : "<<total_hadInt_cog<<endl;
  cout<<"***********"<<endl;


  cout<<"======> leading neutral pion fraction <======="<<endl;
  cout<<"Fraction : "<<(total_leading_pi0*100.0)/hadronic_int<< endl;


  if(E_beam < 0) {
    cout<<"E_beam negative!!!"<<endl;
    return;
  }


  
  

  cout<<endl<<endl;

}


