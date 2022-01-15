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
#include "TStopwatch.h"
#include "../interface/AnalyzeHGCOctTB.h"
#include "Math/SMatrix.h"
#include "Math/SVector.h"
#include <math.h>
#include <TF1.h>

using namespace std;



int main(int argc, char* argv[]) { //main function: takes argument at time the times of execution

  if (argc < 3) {
    cerr << "Please give 4 arguments: runList  outputFileName   dataset   energy" <<endl;

    return -1;
  }
  const char *inputFileList = argv[1];
  const char *outFileName   = argv[2];
  const char *data          = argv[3];
  const char *energy = argv[4];
  const char *config          = "alpha";

  bool shower = true;
  bool longi = false;
  bool trans = false;
  
  bool UseSSBasedCat = true;
  int noise_FH_ = 4;


  AnalyzeHGCOctTB hgcOctTB(inputFileList, outFileName, data, config, energy, shower, longi, trans);   //intialzing object
  
  cout << "dataset " << data << " " << endl;
  cout << "energy " << energy << " " << endl;


  hgcOctTB.EventLoop(data,UseSSBasedCat,noise_FH_);   //Call event loop function
  return 0;
}



void AnalyzeHGCOctTB::EventLoop(const char *data, bool UseSSBasedCat, int noise_FH_) {  //Event loop function


  if (fChain == 0) return;

  Long64_t nentries = fChain->GetEntriesFast();

  cout << "nentries " << nentries << endl;
  cout << "Analyzing dataset " << data << " " << endl;

  cout<<endl;
  if(SHOWER_RECO_HIST) std::cout<<BOLDGREEN<<"Shower energy histograms will be filled"<<RESET<<std::endl;
  if(LONGI_PROFILE_HIST) std::cout<<BOLDGREEN<<"Longitudinal histograms will be filled"<<RESET<<std::endl;
  if(TRANSVERSE_PROFILE_HIST)  {
    std::cout<<BOLDGREEN<<"Transverse histograms will be filled: with track impact as center"<<RESET<<std::endl;
  }
  if(!SHOWER_RECO_HIST && !LONGI_PROFILE_HIST && !TRANSVERSE_PROFILE_HIST) std::cout<<BOLDRED<<"WARNING: NONE OF THE HISTOGRAMS SHALL BE FILLED"<<RESET<<std::endl;
  cout<<endl;


  std::cout<<"Event Categorization based on :"<<BOLDGREEN<<" Shower start Algorithm"<<RESET<<std::endl;

  bool ScaleSim = false;
  float ee_rescaling = 1.035;
  float fh_rescaling = 1.095;
  float ah_rescaling = 1.095;

  if(!strcmp(data,"data")) {
    ScaleSim = false;
    ee_rescaling = 1.0;
    fh_rescaling = 1.0;
    ah_rescaling = 1.0;

  }
  else {
    if(ScaleSim) {
      cout<<"CHECK IF YOU REALLY NEED RESCALING, THEN COME BACK!!!"<<endl;
      return;
      //std::cout<<BOLDGREEN<<"Sim rescaling ON, EE with factor : "<<ee_rescaling<<" ; FH with factor : "<<fh_rescaling<<" ; AH with factor : "<<ah_rescaling<<RESET<<std::endl;
      
    }
    else {
      std::cout<<BOLDRED<<"Sim rescaling OFF"<<RESET<<std::endl;
      std::cout<<BOLDGREEN<<"Sim has already been rescaled!! EE with factor : "<<ee_rescaling<<" ; FH with factor : "<<fh_rescaling<<" ; AH with factor : "<<ah_rescaling<<RESET<<std::endl;
      ee_rescaling = 1.0;
      fh_rescaling = 1.0;
      ah_rescaling = 1.0;
    }
    
  }


  // std::cout<<std::endl<<BOLDRED<<"Doing some debugging for transverse shower shape profiles. Search _check_TS_ for extra cout statements"<<RESET<<std::endl<<std::endl;

  Long64_t nbytes = 0, nb = 0;
  Long64_t region_1_classified_events = 0;
  Long64_t region_2_classified_events = 0;
  Long64_t region_3_classified_events = 0;
  Long64_t region_4_classified_events = 0;
  Long64_t non_classified_events = 0;
  int decade = 0;
  Long64_t ahc_zeroHit = 0;
  Long64_t EE_zeroHit = 0;
  Long64_t FH_zeroHit = 0;
  Long64_t EE_zeroHit_R2 = 0;
  int chi2_method = 3; // to do: remove this

  //To do : Update the function parameters
  //function parateres for fitted weights: energy dependent weights (to be used for reco energy as input)
  // TF1* f_EH_w1 = new TF1("f_EH_w1","[0]*[0]+[1]*[1]/sqrt(x)", 20, 300);
  // f_EH_w1->FixParameter(0,1.025153);
  // f_EH_w1->FixParameter(1,1.624093);
  // TF1* f_EH_w2 = new TF1("f_EH_w2","[0]*[0]+[1]*[1]/sqrt(x)", 20, 300);
  // f_EH_w2->FixParameter(0,0.958412);
  // f_EH_w2->FixParameter(1,1.226055);
  // TF1* f_EH_w3 = new TF1("f_EH_w3","sqrt([0]*[0]+[1]*[1]/x)", 20, 300);
  // f_EH_w3->FixParameter(0,1.019201);
  // f_EH_w3->FixParameter(1,3.796437);

  // TF1* f_H_w1 = new TF1("f_H_w1","[0]*x", 20, 300);
  // f_H_w1->FixParameter(0,0.0);
  // TF1* f_H_w2 = new TF1("f_H_w2","[0]*[0]+[1]*[1]/sqrt(x)", 20, 300);
  // f_H_w2->FixParameter(0,0.908137);
  // f_H_w2->FixParameter(1,0.995139);
  // TF1* f_H_w3 = new TF1("f_H_w3","sqrt([0]*[0]+[1]*[1]/x)", 20, 300);
  // f_H_w3->FixParameter(0,0.977450);
  // f_H_w3->FixParameter(1,2.996701);


    ///////////////////////////////////////////
  // Updated weights, and functional form //
  /////////////////////////////////////////

  TF1* f_EH_w1 = new TF1("f_EH_w1","[0] + [1]/sqrt(x)", 5, 320);
  f_EH_w1->FixParameter(0,1.05);
  f_EH_w1->FixParameter(1,2.65);
  TF1* f_EH_w2 = new TF1("f_EH_w2","[0] + [1]/sqrt(x)", 5, 320);
  f_EH_w2->FixParameter(0,0.93);
  f_EH_w2->FixParameter(1,1.49);
  TF1* f_EH_w3 = new TF1("f_EH_w3","[0] + [1]/sqrt(x)", 5, 320);
  f_EH_w3->FixParameter(0,0.96);
  f_EH_w3->FixParameter(1,1.34);

  TF1* f_H_w1 = new TF1("f_H_w1","[0]*x", 5, 320);
  f_H_w1->FixParameter(0,0.0);
  TF1* f_H_w2 = new TF1("f_H_w2","[0] + [1]/sqrt(x)", 5, 320);
  f_H_w2->FixParameter(0,0.83);
  f_H_w2->FixParameter(1,1.01);
  TF1* f_H_w3 = new TF1("f_H_w3","[0] + [1]/sqrt(x)", 5, 320);
  f_H_w3->FixParameter(0,0.94);
  f_H_w3->FixParameter(1,0.92);

  //////////////////  

  //For Track window selection cut
  std::map<int,vector<double>> track_window;
  vector<double> temp_window;

  /////////////////////////////////
  //    For track window @ L1  ///
  ////////////////////////////////
  
  temp_window.clear();
  temp_window.push_back(-2.60);temp_window.push_back(1.80);temp_window.push_back(-1.20);temp_window.push_back(3.20);
  track_window[20] = temp_window;

  temp_window.clear();
  temp_window.push_back(-2.65);temp_window.push_back(1.80);temp_window.push_back(-1.20);temp_window.push_back(3.20);
  track_window[50] = temp_window;

  temp_window.clear();
  temp_window.push_back(-2.65);temp_window.push_back(1.80);temp_window.push_back(-1.20);temp_window.push_back(3.20);
  track_window[80] = temp_window;

  temp_window.clear();
  temp_window.push_back(-2.65);temp_window.push_back(1.80);temp_window.push_back(-1.20);temp_window.push_back(3.20);
  track_window[100] = temp_window;

  temp_window.clear();
  temp_window.push_back(-2.65);temp_window.push_back(1.60);temp_window.push_back(-1.00);temp_window.push_back(3.00);
  track_window[120] = temp_window;


  temp_window.clear();
  temp_window.push_back(-3.20);temp_window.push_back(-0.90);temp_window.push_back(0.00);temp_window.push_back(2.50);
  track_window[200] = temp_window;

  temp_window.clear();
  temp_window.push_back(-2.90);temp_window.push_back(0.00);temp_window.push_back(-0.30);temp_window.push_back(2.10);
  track_window[250] = temp_window;

  temp_window.clear();
  temp_window.push_back(-3.10);temp_window.push_back(-0.90);temp_window.push_back(-0.40);temp_window.push_back(2.00);
  track_window[300] = temp_window;

  ///////////////////////////////////////////////////////////////////

  //Following code snippet takes value from initialized layer-map text file : "config1_lengths.txt"
  float lambda[79];
  for(int i = 0; i < 79; i++) {
    lambda[i] = layer_positions[i+1].at(2); //lambda_int
    

  }





  bool DEBUG = false;  //IMPORTANT FOR DEBUGGING
  int debug_count = 0;
  Long64_t event_count[10];
  for (int i = 0; i < 10; i++) event_count[i] = 0;
  Long64_t nEvents = 0;
  Long64_t MIP_pions = 0;
  int TOTAL_ACTIVE_LAYER = 40; //in HGCAL prototype
  int EE_LAYER = 28;
  int FH_LAYER = 12;

  float FH_AH_relative_scale = 0.4;
  float alpha_ = FH_AH_relative_scale;
  float EE_scale = 94.624; //MIPs per GeV
  // float FH_AH_scale = 12.788; //MIPs per GeV
  float FH_AH_scale = 12.66; //MIPs per GeV


  double E_beam = -1.0;


  if(DEBUG) cout<<"DEBUG: Configuration = "<<conf_<<endl;
  if(DEBUG) cout<<"DEBUG: TOTAL_ACTIVE_LAYER = "<<TOTAL_ACTIVE_LAYER<<endl;
  if(DEBUG) cout<<"DEBUG: EE_LAYER = "<<EE_LAYER<<endl;
  if(DEBUG) cout<<"DEBUG: FH_LAYER = "<<FH_LAYER<<endl;

  if(DEBUG) cout<<"DEBUG: Entering event Loop"<<endl;

  Long64_t jentry;

  TStopwatch sw; // to clock the execution run time
  sw.Start();

  double time_hgc = 0.0;
  double time_ahc = 0.0;
  double time_CompartCalc = 0.0;
  double time_LongiCalc = 0.0;
  double time_EnMeasureCalc = 0.0;
  
  for (jentry=0; jentry<nentries;jentry++) {
    // TStopwatch sw;
    // sw.Start();

    // ==============print number of events done == == == == == == == =
    double progress = 10.0 * jentry / (1.0 * nentries);
    int k = int (progress);
    if (k > decade)
      cout << 10 * k << " %" << endl;
    decade = k;


    // ===============read this entry == == == == == == == == == == ==

    Long64_t ientry = LoadTree(jentry);
    // cout<<" jentry : ientry :: "<<jentry<<" : "<<ientry<<endl; 
    if (ientry < 0) {  cout<<"Breaking"<<endl; break;}
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    
    h_nTracks->Fill(ntracks);
    event_count[0]++;
    
    if(!isGoodTrack) continue;  // track quality cut
    event_count[1]++;

    E_beam = beamEnergy;
    h_beamEnergy->Fill(E_beam);
    h_particleID->Fill(pdgID);
    h_runNumber->Fill(run);

    if(pdgID == 11) {
      double trackx = TrackImpactX_layer->at(2);
      double tracky = TrackImpactY_layer->at(2);
      trackx = -1*trackx;
      tracky = -1*tracky;
      
      std::pair<float, float> dxy_al = dxy_alignment(3);
      float dx_corr = dxy_al.first;
      float dy_corr = dxy_al.second; 
   
      trackx = trackx + dx_corr;
      tracky = tracky + dy_corr;
      if(abs(trackx) > 2.0 || abs(tracky) > 2.0) continue;

      if(isHGC_AHC_sync == 0) continue;
      if(MuonVeto == 1) continue;
      if(rechit_shower_start_layer > 3) continue;

    }

    
    
    if(isFHNoisy) continue;  // skip events where there are abnormally high hits in one of the modules at FH layer 8 (< 1% events)
    event_count[2]++;


    ////// Track information , pickup from tree
    double track_x[79];
    double track_y[79];

    std::vector<RecHit*> rechit_collection[40];  //rechit collection for HGC
    std::vector<RecHit_AHCAL*> rechit_AH_collection[39]; //rechit collection for AHC
    Long_t Nrechit_FH_module[12][7];  
    Double_t RechitsEn_FH_module[12][7];

    if(TrackImpactX_layer->size() != 79 || TrackImpactX_layer->size() != 79) {  //santiy check
      cout<<"Trackvector size is NOT 79 for jentry "<<jentry<<"  Breaking !!!"<<endl;
      break;
      
    }

    Long_t Nrechit_EE = 0;
    Long_t Nrechit_FH = 0;
    Double_t rechitEnergySum_EE = 0.0;
    Double_t rechitEnergySum_FH = 0.0;
    double cogX_AH[39];  // will be used to store  COG of AHCAL
    double cogY_AH[39]; // will be used to  store COG of AHCAL
    double rechitEnergyLayer_EE[28];
    double rechitEnergyLayer_FH[12];
    double rechitEnergyLayer_AH[39];

    
    for(int ii = 0; ii < 79; ii++) { // load event info from tree + initialize variables
      track_x[ii] = TrackImpactX_layer->at(ii);
      track_y[ii] = TrackImpactY_layer->at(ii);
      
      if(ii < 40) {

	if(ii < 28) {
	  rechit_energyPerLayer->at(ii) = rechit_energyPerLayer->at(ii)/ee_rescaling;
	  Nrechit_EE += rechit_nHitsPerLayer->at(ii);
	  rechitEnergySum_EE += rechit_energyPerLayer->at(ii);
	  rechitEnergyLayer_EE[ii] = 0.0;
	}
	else {
	  rechit_energyPerLayer->at(ii) = rechit_energyPerLayer->at(ii)/fh_rescaling;
	  Nrechit_FH += rechit_nHitsPerLayer->at(ii);
	  if(ii+1 != 36)
	    rechitEnergySum_FH += rechit_energyPerLayer->at(ii);
	  
	}
	
	rechit_collection[ii].clear();
	if(ii < 39) {
	  ahc_energyPerLayer->at(ii) = ahc_energyPerLayer->at(ii)/ah_rescaling;
	  if(ii+1 == 38) {
	    ahc_nHitsPerLayer->at(ii) = 0.0;
	    ahc_energyPerLayer->at(ii) = 0.0;
	  }
	  rechit_AH_collection[ii].clear();
	}
	if(ii < 12) {
	  rechitEnergyLayer_FH[ii] = 0.0;
	  for(int j = 0; j < 7; j++) { 
	    Nrechit_FH_module[ii][j] = 0; 
	    RechitsEn_FH_module[ii][j] = 0; 
	  }
	}
	
      }

      if(ii < 39) {
	cogX_AH[ii] = 0.0;
	cogY_AH[ii] = 0.0;
	rechitEnergyLayer_AH[ii] = 0.0;
      }
    }

    ////////////////////////////////////////////////////////
    //////// Updated track window cut at CE-E layer 1 //////
    ////////////////////////////////////////////////////////

    double trackX_al_1 = track_x[0];
    double trackY_al_1 = track_y[0];
    if(!strcmp(data,"data")) { // Track alignment for data [Not needed for sim]
	trackX_al_1 = -1.0*track_x[0];
	trackY_al_1 = -1.0*track_y[0];
	trackX_al_1 += dxy_alignment(1).first;
	trackY_al_1 += dxy_alignment(1).second;;
      }
    vector<double> track_cut = track_window[(int)E_beam];
    if((trackX_al_1 < track_cut.at(0) || trackX_al_1 > track_cut.at(1) || trackY_al_1 < track_cut.at(2) || trackY_al_1 > track_cut.at(3))) {
      isInTrackWindow = false;
    }

      ////////////////////////////////////////////////////////////
    if(strcmp(data,"data") != 0 && Nrechit_EE == 0) continue;  // For protection (was needed in earlier sim samples)

    ////////////////////////////////////////////////////////////

    if(DEBUG) { cout << " DEBUG :: SW >>>>>> Starting of HGCAL loop, Took: "; sw.Print("u");  }

    TStopwatch sw_hgc;
    sw_hgc.Start();

    double energyInFHLayer8 = 0.0;
    
    ////////////////////////////////////////////
    //            HGCAL Part                  //
    ////////////////////////////////////////////

    /// Read Tree
    if(DEBUG) cout<<"DEBUG: Start Analylizing RecHits!!"<<endl;
    
    std::vector<int> temp_moduleID;
    

    double cog_z_EE = 0.0;
    double cog_z_FH = 0.0;
    double cog_z_HGCAL = 0.0;

    double cog_layer_EE = 0.0;
    double cog_layer_FH = 0.0;
    double cog_layer_HGCAL = 0.0;


    Double_t rechitEnergySum = 0.0;
    int module_part_ = -1;
    int module_layer_ = -1;
    int module_position_ = -1;

    double noise_sigma_ = -1.0; 

    for(int i = 0 ; i < NRechits; i++){
      temp_moduleID.clear();
      
      int temp_layer = rechit_layer->at(i);
      int temp_chip = rechit_chip->at(i);
      int temp_channel = rechit_channel->at(i);
      int en_chan = temp_chip*1000+temp_channel;
      
      if(temp_layer <= 28) rechit_energy->at(i) = rechit_energy->at(i)/ee_rescaling;
      else                 rechit_energy->at(i) = rechit_energy->at(i)/fh_rescaling;
      //channel masking
      if(hgc_channel_mask->at(i)) continue;

      //Noise cut ( pass_noise_thres == 1 => indicates it is above noise cut threshold => not noisy
      if(!pass_noise_thres->at(i)) continue;


      //to do: no alignment needed?
      
      // double trackx = track_x[temp_layer-1];
      // double tracky = track_y[temp_layer-1];
      
      // if(!strcmp(data,"data")) {
      // 	trackx = -1*trackx;
      // 	tracky = -1*tracky;

      // 	std::pair<float, float> dxy_al = dxy_alignment(temp_layer);
      // 	float dx_corr = dxy_al.first;
      // 	float dy_corr = dxy_al.second; 
      // 	recx -= dx_corr;
      // 	recy -= dy_corr;
      // }


      temp_moduleID = getModuleLocation(rechit_module->at(i)); //get module location from module map file of config-1: "../txt_maps/config_maps/moduleMAP_config1.txt"

      if(!temp_moduleID.size() || temp_moduleID.size()<3) { // protection
	cout<<"ERROR: Could NOT locate MODULE location for module "<<rechit_module->at(i)<<endl;
	cout<<"\t more info temp_moduleID.size() = "<<temp_moduleID.size()<<endl;
	return;
      }

      module_part_ = temp_moduleID.at(0); // EE = 0 or FH = 1
      module_layer_ = temp_moduleID.at(1); 
      module_position_ = temp_moduleID.at(2); // 0 for EE, 1 to 7 for FH

      // bool iscentral = false;
      // if(module_position_ == 0 || module_position_ == 4) iscentral = true;
      // if(!iscentral) continue;


      // cout<<module_position_<<endl;
      // if(module_position_ > 0 && module_position_ != 4) {cout << "module posistion not 4 or 0."<<endl;}
      h_moduleID->Fill(module_part_);

            
      //////////////////////////////////////////////////////////////////////////

      if(rechit_layer->at(i) > 28) { //FH
	Nrechit_FH_module[rechit_layer->at(i)-29][module_position_-1]++;
	RechitsEn_FH_module[rechit_layer->at(i)-29][module_position_-1] += rechit_energy->at(i);
	rechitEnergyLayer_FH[rechit_layer->at(i)-29] += rechit_energy->at(i);
	if(rechit_layer->at(i) == 36) {
	  if(!(rechit_module->at(i) == 38 && (en_chan == 16 || en_chan == 3004)))
	    energyInFHLayer8 += rechit_energy->at(i)/fh_rescaling;
	  
	}
      }
      else { //EE
	rechitEnergyLayer_EE[rechit_layer->at(i)-1] += rechit_energy->at(i);
      }
      if(module_part_ == 0) {
	cog_z_EE += rechit_energy->at(i)*rechit_z->at(i);
	cog_layer_EE += rechit_energy->at(i)*rechit_layer->at(i);
	cog_z_HGCAL += rechit_energy->at(i)*rechit_z->at(i);
	cog_layer_HGCAL += rechit_energy->at(i)*rechit_layer->at(i);

      }
      else if(module_part_ == 1) {
	cog_z_FH += rechit_energy->at(i)*rechit_z->at(i);
	cog_layer_FH += rechit_energy->at(i)*rechit_layer->at(i);
	cog_z_HGCAL += rechit_energy->at(i)*rechit_z->at(i);
	cog_layer_HGCAL += rechit_energy->at(i)*rechit_layer->at(i);
	
      }
      else { // Protection
	cout<<"ERROR: Unknown Module Part detected!!!!"<<endl;
	return;
	
      }

      rechitEnergySum+=rechit_energy->at(i);



      RecHit* rechit_temp = new RecHit(temp_layer, temp_chip, temp_channel, en_chan, module_part_, module_position_, rechit_module->at(i), (int)rechit_iu->at(i), (int)rechit_iv->at(i), rechit_x->at(i), rechit_y->at(i), rechit_z->at(i), rechit_energy->at(i));       // for HGCAL rechit collection
      rechit_collection[temp_layer-1].push_back(rechit_temp);       // for HGCAL rechit collection
      
    }

    //correcting rechitHit_energy per layer for FH8_, HGCAL L36
    rechit_energyPerLayer->at(36-1) = energyInFHLayer8;
    rechitEnergySum_FH += energyInFHLayer8;
    
    ///////////////////////////////////////////
    ////////// End of HGCAL nRechit Loop /////
    //////////////////////////////////////////

    // cout<<" *************** "<<endl;
    // for(int i = 0; i < 40; i++) {
    //   if(i+1 != 36) continue;
    //   if(i < 28)
    // 	cout<<"jentry : layer : rechit_energyPerLayer : rechitEnergyLayer_EE :: "<<jentry<<" : "<<i+1<<" : "<<rechit_energyPerLayer->at(i)<<" : "<<rechitEnergyLayer_EE[i]<<endl;
    //   else
    // 	cout<<"jentry : layer : rechit_energyPerLayer : rechitEnergyLayer_FH :: "<<jentry<<" : "<<i+1<<" : "<<rechit_energyPerLayer->at(i)<<" : "<<rechitEnergyLayer_FH[i-28]<<endl;
	
    // }
    
    
    time_hgc += sw_hgc.RealTime();
    sw_hgc.Stop();

    if(DEBUG) { cout << "SW >>>>>>END of HGCAL loop; Took: "; sw.Print();  }


    TStopwatch sw_ahc;
    sw_ahc.Start();

    ////////////////////////////////////////////
    //            AHCAL Part                  //
    ////////////////////////////////////////////
    
    Double_t rechitEnergySum_AH = 0.0;
    Double_t cog_layer_AH = 0.0;
    Long_t Nrechit_AH = ahc_nHits;
    for(int i = 0 ; i < ahc_nHits; i++) {
      int temp_layer = ahc_hitK->at(i);
      
      ahc_hitEnergy->at(i) = ahc_hitEnergy->at(i)/ah_rescaling;
      
      rechitEnergyLayer_AH[temp_layer-1] += ahc_hitEnergy->at(i);
      //Skip AHCAL layer 38
      if(temp_layer == 38) continue;
      // if(ahc_channel_mask->at(i)) continue;  //to do: Why not using masking here?
      
      rechitEnergySum_AH += ahc_hitEnergy->at(i);
      cog_layer_AH += ahc_hitEnergy->at(i)*temp_layer;
      cogX_AH[temp_layer-1] += (ahc_hitX->at(i)*ahc_hitEnergy->at(i));  // calculating AHCAL COG X
      cogY_AH[temp_layer-1] += (ahc_hitY->at(i)*ahc_hitEnergy->at(i)); // calculating AHCAL COG Y
      
      RecHit_AHCAL* rechit_temp = new RecHit_AHCAL(temp_layer, ahc_hitI->at(i), ahc_hitJ->at(i), ahc_hitX->at(i), ahc_hitY->at(i), ahc_hitZ->at(i), ahc_hitEnergy->at(i)); // for AHCAL rechit collection
      rechit_AH_collection[temp_layer-1].push_back(rechit_temp); // for AHCAL rechit collection
    }


    //////////////////////////////////////////////
    ////////// End of AHCAL nRechit Loop  ////////
    //////////////////////////////////////////////

    time_ahc += sw_ahc.RealTime();
    sw_ahc.Stop();

    
    if(DEBUG) { cout << "SW >>>>>>END of AHCAL loop; Took: "; sw_ahc.Print();  }
    
    /////////////////////////////////////////////
    //// Shower start finder information   /////
    ////////////////////////////////////////////

    bool MIP = true;
    int shower_start_index = -1;
    float shower_lambda_ = -1.0;
    
    if(rechit_shower_start_layer < 0 ) { // loading shower start layer from the tree
      shower_start_index = -1;
      shower_lambda_ = -1;
      MIP = true;

    }
    else {
      shower_start_index = rechit_shower_start_layer -1 ;
      shower_lambda_ = lambda[rechit_shower_start_layer - 1];
      MIP = false;
    }


    ///////////////////////////////////////////////////////


    TStopwatch sw_compart;
    sw_compart.Start();

    
    ////////////////////////////////////////////////////////////////
    //                Comapartment wise calculation              //
    ///////////////////////////////////////////////////////////////

    // This section assigns values to layer-variables in a layerwise
    // loop after the calculation is done in rechit-loop.

    // All the further calculations related to layers are done here, such
    // as COG caluculation, lateral shower variables, seeds etc.

    float seed_energy_incl = 0.0;    //to do : can be deleted?
    float r1_energy_incl = 0.0;     //to do : can be deleted?
    float r2_energy_incl = 0.0;     //to do : can be deleted?
    float E1_inc = 0.0; float E7_inc = 0.0; float E19_inc = 0.0;
    float E1_EE = 0.0; float E7_EE = 0.0; float E19_EE = 0.0;
    float E1_FH = 0.0; float E7_FH = 0.0; float E19_FH = 0.0;
    float S1_AH = 0.0; float S9_AH = 0.0; float S25_AH = 0.0;


    double total_energy_EE = rechitEnergySum_EE; //to do: can be moved down?
    double total_energy_FH = rechitEnergySum_FH; //to do: can be moved down?
    double total_energy_AH = rechitEnergySum_AH; //to do: can be moved down?


    vector< std::pair<double,double> > seed_xy;
    seed_xy.clear();

    
    for(int iL = 0; iL < TOTAL_ACTIVE_LAYER; iL++){


      // if(iL > EE_LAYER -1) { // for FH DEBUG 
      // 	for(int j = 0; j < 7; j++) {
      // 	  if(iL+1 >= 38 && j+1 != 4) continue;
      // 	  h_nrechit_FH[iL-EE_LAYER][j]->Fill(Nrechit_FH_module[iL-EE_LAYER][j]);
      // 	  h_energy_FH[iL-EE_LAYER][j]->Fill(RechitsEn_FH_module[iL-EE_LAYER][j]);
      // 	}

      // }

      // For layerwise distributions
      h_ld_energy[iL]->Fill(rechit_energyPerLayer->at(iL));
      h_ld_nrechit[iL]->Fill(rechit_nHitsPerLayer->at(iL));

      h_cogX[iL]->Fill(rechit_cogX->at(iL));
      h_cogY[iL]->Fill(rechit_cogY->at(iL));
      h_cogX_cogY[iL]->Fill(rechit_cogX->at(iL),rechit_cogY->at(iL));

      h_trackX[iL]->Fill(track_x[iL]);
      h_trackY[iL]->Fill(track_y[iL]);
      h_trackX_trackY[iL]->Fill(track_x[iL],track_y[iL]);



      // For lateral shower variable calculation


      // cout<<"size at layer "<<iL+1<<" : "<<rechit_collection[iL].size()<<endl;
      LateralRings* temp_LR = new LateralRings(rechit_collection[iL]); // see LateralRings.h file
      if(!rechit_collection[iL].size()) { 
	temp_LR->setLayer(iL+1);
      }
      
      float seed = temp_LR->getSeedEnergy();   // get seed-energy for this layer
      int seed_index = temp_LR->getSeedIndex();
      float r1 = temp_LR->getEnergyInRing1(); // get energy in 1st ring around the seed
      float r2 = temp_LR->getEnergyInRing2();  // get energy in 2nd ring around the seed
      float E1 = seed; float E7 = E1+r1; float E19 = E7+r2;  // sum them up to define E1, E7 etc variables

      if( seed_index < 0) seed_xy.push_back(std::make_pair(-999.0,-999.0));  // if no seed found (i.e. zero energy in this layer), then assign dummy seed coordinates
      else {
	vector<float> temp = ((rechit_collection[iL]).at(seed_index))->getCellCoordinate_xyz();
	if((int)temp.size() !=3) { cout<<"halting at point tango"<<endl; return;} // protection
	seed_xy.push_back(std::make_pair(temp[0],temp[1]));  // seed coordinates
      }
      
      seed_energy_incl += seed;  //to do: can be deleted?
      r1_energy_incl += r1; //to do: can be deleted?
      r2_energy_incl += r2; //to do: can be deleted?
      E1_inc += seed; //to do: can be deleted?
      E7_inc += E7; //to do: can be deleted?
      E19_inc += E19; //to do: can be deleted?
      
      if(iL < 28) {  // sum lateral variables for whole EE
	E1_EE += seed;
	E7_EE += E7;
	E19_EE += E19;

      }
      else {   // sum lateral variables for whole FH
	E1_FH += seed;
	E7_FH += E7;
	E19_FH += E19;
      }

      if(iL < 39) {  // caluculation for AHCAL
	cogX_AH[iL] = cogX_AH[iL]/rechitEnergyLayer_AH[iL];
	cogY_AH[iL] = cogY_AH[iL]/rechitEnergyLayer_AH[iL];
	
	LateralSquares* temp_LS = new LateralSquares(rechit_AH_collection[iL]); // see LateralSquares.h file

	
	float seed_AH = temp_LS->getSeedEnergy(); // ==> 1 cell
	float s1 = temp_LS->getEnergyInSquare1(); // ==> 8 cells
	float s2 = temp_LS->getEnergyInSquare2(); // ==> 16 cells
	float S1 = seed_AH; float S9 = S1+s1; float S25 = S9+s2;
	S1_AH += seed_AH;
	S9_AH += S9;
	S25_AH += S25;
	// cout<<iL+1<<" seed_AH: "<<seed_AH<<" ,s1: "<<s1<<" ,s2: "<<s2<<" :: S1:"<<S1<<" ,S9: "<<S9<<" ,S25: "<<S25<<endl;
	
	// For Layerwise distribution
	h_ld_energy[iL+40]->Fill(ahc_energyPerLayer->at(iL));
	h_ld_nrechit[iL+40]->Fill(ahc_nHitsPerLayer->at(iL));

      }
      

      ///////  END of lateral information calculation ///////////




    }

    
    ///////////////////////////////////////////////////////////////
    ///////     END of compartment wise calculation       /////////
    ///////////////////////////////////////////////////////////////

    

    time_CompartCalc += sw_compart.RealTime();
    sw_compart.Stop();

    if(DEBUG) { cout << "SW >>>>>>END of compartment wise calculation loop; Took: "; sw_compart.Print();  }

    cog_z_EE = cog_z_EE/rechitEnergySum_EE;
    cog_layer_EE = cog_layer_EE/rechitEnergySum_EE;

    cog_z_FH = cog_z_FH/rechitEnergySum_FH;
    cog_layer_FH = cog_layer_FH/rechitEnergySum_FH;
    cog_layer_FH = cog_layer_FH - 28;

    cog_z_HGCAL = cog_z_HGCAL/rechitEnergySum;
    cog_layer_HGCAL = cog_layer_HGCAL/rechitEnergySum;

    cog_layer_AH = cog_layer_AH/rechitEnergySum_AH;

    // Fill histograms with calculated variables
    
    h_cog_layer_EE->Fill(cog_layer_EE);
    h_cog_layer_FH->Fill(cog_layer_FH);
    h_cog_layer_HGCAL->Fill(cog_layer_HGCAL);
    h_cog_layer_AH->Fill(cog_layer_AH);

    
    if(E1_EE !=0 ) {
      h_E1_E7_EE->Fill(E1_EE/E7_EE);
      h_E1_E19_EE->Fill(E1_EE/E19_EE);
      h_E7_E19_EE->Fill(E7_EE/E19_EE);
    }
    if(E1_FH != 0) {
      h_E1_E7_FH->Fill(E1_FH/E7_FH);
      h_E1_E19_FH->Fill(E1_FH/E19_FH);
      h_E7_E19_FH->Fill(E7_FH/E19_FH);
    }
    if(E1_inc != 0) {
      h_E1_E7_HGCAL->Fill(E1_inc/E7_inc);
      h_E1_E19_HGCAL->Fill(E1_inc/E19_inc);
      h_E7_E19_HGCAL->Fill(E7_inc/E19_inc);
    }

    if(!MuonVeto) {
      h_S1_S9_AH->Fill(S1_AH/S9_AH);
      h_S1_S25_AH->Fill(S1_AH/S25_AH);
      h_S9_S25_AH->Fill(S9_AH/S25_AH);
    }


    if(Nrechit_EE == 0) EE_zeroHit++;
    if(Nrechit_FH == 0) FH_zeroHit++;
    if(Nrechit_AH == 0) ahc_zeroHit++;



    // if(shower_start_index+1 == 15) h_run_SS15->Fill(run);  //to do:  Not needed anymore
    h_shower_start->Fill(shower_lambda_);

    

    if(shower_start_index+1 > 0.0 && shower_start_index+1 <= 28) {
      if(E1_EE != 0) {
	h_E1_E7_SS_EE->Fill(E1_EE/E7_EE);
	h_E1_E19_SS_EE->Fill(E1_EE/E19_EE);
	h_E7_E19_SS_EE->Fill(E7_EE/E19_EE);
      }
    }
    else if(shower_start_index+1 > 28) {
      if(E1_EE != 0) {
	h_E1_E7_MIPs_in_EE->Fill(E1_EE/E7_EE);
	h_E1_E19_MIPs_in_EE->Fill(E1_EE/E19_EE);
	h_E7_E19_MIPs_in_EE->Fill(E7_EE/E19_EE);
      }
      if(E1_FH != 0) {
	h_E1_E7_SS_FH->Fill(E1_FH/E7_FH);
	h_E1_E19_SS_FH->Fill(E1_FH/E19_FH);
	h_E7_E19_SS_FH->Fill(E7_FH/E19_FH);
      }
    }
    else {
      if(E1_FH != 0) {
	h_E1_E7_MIPs_in_FH->Fill(E1_FH/E7_FH);
	h_E1_E19_MIPs_in_FH->Fill(E1_FH/E19_FH);
	h_E7_E19_MIPs_in_FH->Fill(E7_FH/E19_FH);
      }
    }


    ///////////////////////////////////////////////
    //                Collapsed EE               //
    ///////////////////////////////////////////////

    // This section is important for shower start algo fits and plots
    // Sums up back-to-back shower start layer in EE to make it more presentable
    
    if(shower_start_index+1 <= 28 && (shower_start_index+1) > 0) {
      h_shower_start_full_collapsed_EE->Fill(lambda[0]);
      if((shower_start_index+1)%2 == 0) {
	h_shower_start_part_collapsed_EE->Fill(lambda[shower_start_index-1]);
	// cout<<shower_start_index<<" "<<(shower_start_index+1)%2<<" "<<shower_start_index-1<<" "<<lambda[(shower_start_index+1)/2]<<endl;
      }
      else {
	h_shower_start_part_collapsed_EE->Fill(lambda[shower_start_index]);
      }
    }
    else {
      h_shower_start_full_collapsed_EE->Fill(shower_lambda_);
      h_shower_start_part_collapsed_EE->Fill(shower_lambda_);
    }

    if(MuonVeto) MIP = true; //to do: WHY NOT USING muon veto cleaning cut on the top?

    TStopwatch sw_longi;
    sw_longi.Start();

    
    
    ////////////////////////////////////////////////////////////
    ////           S H O W E R       P R O F I L E          ////
    ///////////////////////////////////////////////////////////

    // This section calculates longitudinal shower variables and fills it in histograms
    
    bool DoWeight = false;
    double beta_ = -1.0;
    double ee_mip2gev = 0.0106; // GeV-per-MIP fixed scale for EE 
    double fh_mip2gev = 0.0789; // GeV-per-MIP fixed scale for FH
    double ah_mip2gev = 0.4*fh_mip2gev; // GeV-per-MIP fixed scale for AH
    
    float w1 = 0; float w2 = 0; float w3 = 0;

    w1 = ee_mip2gev; w2 = fh_mip2gev; w3 = ah_mip2gev;


    double total_energy = total_energy_EE + total_energy_FH + total_energy_AH;
    double total_energy_gev = w1*total_energy_EE + w2*total_energy_FH + w3*total_energy_AH;

    ////// Energy dependent weighting, with boot-strap energy determination
    // to do: use beam energy as input !!!!
    if(shower_start_index+1 >= 1 && shower_start_index+1 <= 28) {
      // float recoE_w1 = f_EH_w1->Eval(total_energy_gev);
      // float recoE_w2 = f_EH_w2->Eval(total_energy_gev);
      // float recoE_w3 = f_EH_w3->Eval(total_energy_gev);
      //w1 = recoE_w1; w2 = recoE_w2; w3 = recoE_w3;
      w1 = getChi2Weights_EH(beamEnergy).at(0);
      w2 = getChi2Weights_EH(beamEnergy).at(1);
      w3 = getChi2Weights_EH(beamEnergy).at(2);

    }
    else if(shower_start_index+1 >= 29) {
      // float recoE_w1 = f_H_w1->Eval(total_energy_gev);
      // float recoE_w2 = f_H_w2->Eval(total_energy_gev);
      // float recoE_w3 = f_H_w3->Eval(total_energy_gev);
      //w1 = recoE_w1; w2 = recoE_w2; w3 = recoE_w3;
      w1 = getChi2Weights_H(beamEnergy).at(0);
      w2 = getChi2Weights_H(beamEnergy).at(1);
      w3 = getChi2Weights_H(beamEnergy).at(2);

    }
    
    total_energy_gev = w1*(ee_mip2gev*total_energy_EE) + w2*(fh_mip2gev*total_energy_FH) + w3*(ah_mip2gev*total_energy_AH);  //calculating total energy in GeV using detector scale energy 

    
    double EE_frac = 100.0*(w1*(ee_mip2gev*total_energy_EE))/total_energy_gev;
    double FH_frac = 100.0*(w2*(fh_mip2gev*total_energy_FH))/total_energy_gev;
    double AH_frac = 100.0*(w3*(ah_mip2gev*total_energy_AH))/total_energy_gev;
      
    if(shower_start_index+1 == 1 && debug_count < 20 && false) {
      cout<<"jentry : EE : FH : AH : total :: "<<EE_frac<<" : "<<FH_frac<<" : "<<AH_frac<<" : "<<(EE_frac+FH_frac+AH_frac)<<endl;
      debug_count++;
    }


    // compartment-wise variables init for transverse-shower-profile
    double sum_dR1_EE = 0.0;
    double sum_dR2_EE = 0.0;
    double sum_dR3_EE = 0.0;
    double sum_dR5_EE = 0.0;
    double sum_dR8_EE = 0.0;
    double sum_dR10_EE = 0.0;
    double sum_dR12_EE = 0.0;
    double sum_dR15_EE = 0.0;
    double sum_dR18_EE = 0.0;
    double sum_dR20_EE = 0.0;

    double sum_dR1_FH = 0.0;
    double sum_dR2_FH = 0.0;
    double sum_dR3_FH = 0.0;
    double sum_dR5_FH = 0.0;
    double sum_dR8_FH = 0.0;
    double sum_dR10_FH = 0.0;
    double sum_dR12_FH = 0.0;
    double sum_dR15_FH = 0.0;
    double sum_dR18_FH = 0.0;
    double sum_dR20_FH = 0.0;


    /// For summed up transverse variables inclusive in 7 EE layers //
    // 10 elements correspond to 10 different radii for energy-sum at i_th layer
    // it makes life easier
    double sum_dRi_01_07[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double sum_dRi_08_14[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double sum_dRi_15_21[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    double sum_dRi_22_28[10] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};


    // float dRs[] = {1.0,2.0,3.0,5.0,8.0,10.0,12.0,15.0,18.0,20.0};

    bool debug_flag = true;
    bool new_event = true;

    // TH1F* htemp_EE;
    // TH1F* htemp_FH;
    // TH1F* htemp_AH;
    // if(!htemp_EE) {
    //   // delete htemp_EE;
    //   htemp_EE = new TH1F("htemp_EE","htemp_EE",500,0,50);
    // }
    // if(!htemp_FH) {
    //   // delete htemp_FH;
    //   htemp_FH = new TH1F("htemp_FH","htemp_FH",500,0,50);
    // }
    // if(!htemp_AH) {
    //   // delete htemp_AH;
    //   htemp_AH = new TH1F("htemp_AH","htemp_AH",500,0,100);
    // }

    TH1F* htemp_EE = new TH1F("htemp_EE","htemp_EE",500,0,50);
    TH1F* htemp_FH = new TH1F("htemp_FH","htemp_FH",500,0,50);
    TH1F* htemp_AH = new TH1F("htemp_AH","htemp_AH",500,0,100);
    
    for(int iL = 0; iL < TOTAL_ACTIVE_LAYER; iL++) { // H E R E    W E    G O
      if(!LONGI_PROFILE_HIST && !TRANSVERSE_PROFILE_HIST) break;
      if(MIP || !isInTrackWindow) continue; //to do: use these cleaning cuts on top
      if(shower_start_index+1 <= 2) continue; // Rejection for pre-showering events


      if(LONGI_PROFILE_HIST) {
	//////////////////////////////////////////////////
	/////////// Longitudinal shower profile part /////
	//////////////////////////////////////////////////
      
	// First we shall fill hists for longi shower profiles

      
	if(iL < EE_LAYER) { //EE part
	
	
	  float EE_energy = rechit_energyPerLayer->at(iL);
	  float EE_detscale = ee_mip2gev*EE_energy;

	  if(shower_start_index+1 == 10  && rechitEnergySum_EE > 0) {  // to check shower profile when shower start at EE-L10
	    h_longi_profile_raw_SS10_check->Fill(iL+1,EE_energy/rechitEnergySum_EE);
	    if(total_energy_gev > 0) h_longi_profile_gev_SS10_check->Fill(iL+1,w1*EE_detscale/total_energy_gev);
	  }
	
	  if(shower_start_index+1 <= 28) { // Showering in CE-E
	    h_longi_profile_ShowerInEE->Fill(lambda[iL],EE_energy);
	    h_longi_profile_ShowerInEE_gev->Fill(lambda[iL],EE_energy*ee_mip2gev);
	
	    h_longi_profile_raw_ShowerInEE_layer->Fill(iL+1,EE_energy);
	    h_longi_profile_gev_ShowerInEE_layer->Fill(iL+1,EE_energy*ee_mip2gev);

	  

	  }
	  else { // MIPs in CE-E
	    h_longi_profile_MipsInEE->Fill(lambda[iL],EE_energy);
	    h_longi_profile_MipsInEE_gev->Fill(lambda[iL],EE_energy*ee_mip2gev);
	    h_longi_profile_MipsInEE_SS_ref->Fill(lambda[iL]-lambda[shower_start_index],EE_energy);

	    h_longi_profile_raw_MipsInEE_layer->Fill(iL+1,EE_energy);
	    h_longi_profile_gev_MipsInEE_layer->Fill(iL+1,EE_energy*ee_mip2gev);
	    h_longi_profile_raw_MipsInEE_SS_ref_layer->Fill(iL+1-shower_start_index,EE_energy);
  
	  }
	  // spandey --> Protection against zero energy in GeV in EE for pions MIPs in EE
	  if(shower_start_index+1 <= 28) {
	    h_longi_profile_gev[shower_start_index]->Fill(lambda[iL],EE_detscale*w1);
	    h_longi_profile_gev_layer[shower_start_index]->Fill(iL+1,EE_detscale*w1);
	    h_longi_profile_gev_fraction[shower_start_index]->Fill(lambda[iL],EE_detscale*w1/total_energy_gev);
	    h_longi_profile_gev_layer_fraction[shower_start_index]->Fill(iL+1,EE_detscale*w1/total_energy_gev);
	  }
	  else {
	    //double mip_energy_per_layer = 0.4/28; // 0.4 GeV in EE compartment divided by 28 EE layers
	    double mip_energy_per_layer = EE_energy*ee_mip2gev; 
	    h_longi_profile_gev[shower_start_index]->Fill(lambda[iL],mip_energy_per_layer);
	    h_longi_profile_gev_layer[shower_start_index]->Fill(iL+1,mip_energy_per_layer);
	    h_longi_profile_gev_fraction[shower_start_index]->Fill(lambda[iL],mip_energy_per_layer/total_energy_gev);
	    h_longi_profile_gev_layer_fraction[shower_start_index]->Fill(iL+1,mip_energy_per_layer/total_energy_gev);

	  }
	  h_longi_profile_raw[shower_start_index]->Fill(lambda[iL],EE_energy);
	  h_longi_profile_raw_layer[shower_start_index]->Fill(iL+1,EE_energy);
	  h_longi_profile_raw_fraction[shower_start_index]->Fill(lambda[iL],EE_energy/total_energy);
	  h_longi_profile_raw_layer_fraction[shower_start_index]->Fill(iL+1,EE_energy/total_energy);
	
	
	  h_longi_profile_inclusive->Fill(lambda[iL],EE_energy);
	  h_longi_profile_inclusive_frac->Fill(lambda[iL],EE_energy/total_energy);

	  h_longi_profile_raw_inclusive_layer->Fill(iL+1,EE_energy);
	  h_longi_profile_raw_inclusive_frac_layer->Fill(iL+1,EE_energy/total_energy);

	  
	  h_Rechits_nrec_SS[shower_start_index][iL]->Fill(rechit_nHitsPerLayer->at(iL));
	  h_Rechits_En_SS[shower_start_index][iL]->Fill(EE_energy);

	  if(shower_start_index+1 <= 7) {  // for summed up shower profile
	    h_longi_profile_raw_inc[0]->Fill(lambda[iL],EE_energy);
	    h_longi_profile_gev_inc[0]->Fill(lambda[iL],EE_detscale*w1);
	    h_longi_profile_raw_layer_inc[0]->Fill(iL+1,EE_energy);
	    h_longi_profile_gev_layer_inc[0]->Fill(iL+1,EE_detscale*w1);
	  }
	  else if(shower_start_index+1 <= 14) {  // for summed up shower profile
	    h_longi_profile_raw_inc[1]->Fill(lambda[iL],EE_energy);
	    h_longi_profile_gev_inc[1]->Fill(lambda[iL],EE_detscale*w1);
	    h_longi_profile_raw_layer_inc[1]->Fill(iL+1,EE_energy);
	    h_longi_profile_gev_layer_inc[1]->Fill(iL+1,EE_detscale*w1);
	  }
	  else if(shower_start_index+1 <= 21) {  // for summed up shower profile
	    h_longi_profile_raw_inc[2]->Fill(lambda[iL],EE_energy);
	    h_longi_profile_gev_inc[2]->Fill(lambda[iL],EE_detscale*w1);
	    h_longi_profile_raw_layer_inc[2]->Fill(iL+1,EE_energy);
	    h_longi_profile_gev_layer_inc[2]->Fill(iL+1,EE_detscale*w1);
	  }
	  else if(shower_start_index+1 <= 28) {  // for summed up shower profile
	    h_longi_profile_raw_inc[3]->Fill(lambda[iL],EE_energy);
	    h_longi_profile_gev_inc[3]->Fill(lambda[iL],EE_detscale*w1);
	    h_longi_profile_raw_layer_inc[3]->Fill(iL+1,EE_energy);
	    h_longi_profile_gev_layer_inc[3]->Fill(iL+1,EE_detscale*w1);
	  }
	  else {  // for summed up shower profile
	    h_longi_profile_raw_inc[shower_start_index-24]->Fill(lambda[iL],EE_energy);
	    h_longi_profile_gev_inc[shower_start_index-24]->Fill(lambda[iL],EE_detscale*w1);
	    h_longi_profile_raw_layer_inc[shower_start_index-24]->Fill(iL+1,EE_energy);
	    h_longi_profile_gev_layer_inc[shower_start_index-24]->Fill(iL+1,EE_detscale*w1);
	  }


	}  // End of EE part
	else { // start of FH part
	  float FH_energy = rechit_energyPerLayer->at(iL);
	  float FH_detscale = fh_mip2gev*FH_energy;
	  if(shower_start_index+1 == 10  && rechitEnergySum_FH > 0) {   // to check shower profile when shower start at EE-L10
	    h_longi_profile_raw_SS10_check->Fill(iL+1,FH_energy/rechitEnergySum_FH);
	    if(total_energy_gev > 0) h_longi_profile_gev_SS10_check->Fill(iL+1,w2*FH_detscale/total_energy_gev);

	  }
	
	  if(shower_start_index+1 <= 28) { // Showering in CE-E
	    h_longi_profile_ShowerInEE->Fill(lambda[iL],FH_energy);
	    h_longi_profile_ShowerInEE_gev->Fill(lambda[iL],FH_energy*fh_mip2gev);

	    h_longi_profile_raw_ShowerInEE_layer->Fill(iL+1,FH_energy);
	    h_longi_profile_gev_ShowerInEE_layer->Fill(iL+1,FH_energy*fh_mip2gev);

	  }
	  else { // MIPs in CE-E
	    h_longi_profile_MipsInEE->Fill(lambda[iL],FH_energy);
	    h_longi_profile_MipsInEE_gev->Fill(lambda[iL],FH_energy*fh_mip2gev);
	    h_longi_profile_MipsInEE_SS_ref->Fill(lambda[iL]-lambda[shower_start_index],FH_energy);

	    h_longi_profile_raw_MipsInEE_layer->Fill(iL+1,FH_energy);
	    h_longi_profile_gev_MipsInEE_layer->Fill(iL+1,FH_energy*fh_mip2gev);
	    h_longi_profile_raw_MipsInEE_SS_ref_layer->Fill(iL+1-shower_start_index,FH_energy);


	  }

	  h_longi_profile_raw[shower_start_index]->Fill(lambda[iL],FH_energy);
	  h_longi_profile_gev[shower_start_index]->Fill(lambda[iL],FH_detscale*w2);

	  h_longi_profile_raw_layer[shower_start_index]->Fill(iL+1,FH_energy);
	  h_longi_profile_gev_layer[shower_start_index]->Fill(iL+1,FH_detscale*w2);

	  h_longi_profile_raw_fraction[shower_start_index]->Fill(lambda[iL],FH_energy/total_energy);
	  h_longi_profile_gev_fraction[shower_start_index]->Fill(lambda[iL],FH_detscale*w2/total_energy_gev);

	  h_longi_profile_raw_layer_fraction[shower_start_index]->Fill(iL+1,FH_energy/total_energy);
	  h_longi_profile_gev_layer_fraction[shower_start_index]->Fill(iL+1,FH_detscale*w2/total_energy_gev);

	  h_longi_profile_inclusive->Fill(lambda[iL],FH_energy);
	  h_longi_profile_inclusive_frac->Fill(lambda[iL],FH_energy/total_energy);


	  h_longi_profile_raw_inclusive_layer->Fill(iL+1,FH_energy);
	  h_longi_profile_raw_inclusive_frac_layer->Fill(iL+1,FH_energy/total_energy);

	  h_Rechits_nrec_SS[shower_start_index][iL]->Fill(rechit_nHitsPerLayer->at(iL));
	  h_Rechits_En_SS[shower_start_index][iL]->Fill(FH_energy);

	  if(shower_start_index+1 <= 7) {
	    h_longi_profile_raw_inc[0]->Fill(lambda[iL],FH_energy);
	    h_longi_profile_gev_inc[0]->Fill(lambda[iL],FH_detscale*w2);
	    h_longi_profile_raw_layer_inc[0]->Fill(iL+1,FH_energy);
	    h_longi_profile_gev_layer_inc[0]->Fill(iL+1,FH_detscale*w2);
	  }
	  else if(shower_start_index+1 <= 14) {
	    h_longi_profile_raw_inc[1]->Fill(lambda[iL],FH_energy);
	    h_longi_profile_gev_inc[1]->Fill(lambda[iL],FH_detscale*w2);
	    h_longi_profile_raw_layer_inc[1]->Fill(iL+1,FH_energy);
	    h_longi_profile_gev_layer_inc[1]->Fill(iL+1,FH_detscale*w2);
	  }
	  else if(shower_start_index+1 <= 21) {
	    h_longi_profile_raw_inc[2]->Fill(lambda[iL],FH_energy);
	    h_longi_profile_gev_inc[2]->Fill(lambda[iL],FH_detscale*w2);
	    h_longi_profile_raw_layer_inc[2]->Fill(iL+1,FH_energy);
	    h_longi_profile_gev_layer_inc[2]->Fill(iL+1,FH_detscale*w2);
	  }
	  else if(shower_start_index+1 <= 28) {
	    h_longi_profile_raw_inc[3]->Fill(lambda[iL],FH_energy);
	    h_longi_profile_gev_inc[3]->Fill(lambda[iL],FH_detscale*w2);
	    h_longi_profile_raw_layer_inc[3]->Fill(iL+1,FH_energy);
	    h_longi_profile_gev_layer_inc[3]->Fill(iL+1,FH_detscale*w2);
	  }
	  else {
	    h_longi_profile_raw_inc[shower_start_index-24]->Fill(lambda[iL],FH_energy);
	    h_longi_profile_gev_inc[shower_start_index-24]->Fill(lambda[iL],FH_detscale*w2);
	    h_longi_profile_raw_layer_inc[shower_start_index-24]->Fill(iL+1,FH_energy);
	    h_longi_profile_gev_layer_inc[shower_start_index-24]->Fill(iL+1,FH_detscale*w2);
	  }

	
	} // End of FH part
	if(iL < 39) { // Use same loop to also fill AHCAL part
	  float AH_energy = ahc_energyPerLayer->at(iL);
	  float AH_detscale = ah_mip2gev*AH_energy;
	  if(shower_start_index+1 == 10 && rechitEnergySum_AH > 0)  {  // to check shower profile when shower start at EE-L10
	    h_longi_profile_raw_SS10_check->Fill(iL+1+40,AH_energy/rechitEnergySum_AH);
	    if(total_energy_gev > 0) h_longi_profile_gev_SS10_check->Fill(iL+1+40,w3*AH_detscale/total_energy_gev);

	  }
	  if(shower_start_index+1 <= 28) { // Showering in CE-E
	    h_longi_profile_ShowerInEE->Fill(lambda[iL+40],AH_energy);
	    h_longi_profile_ShowerInEE_gev->Fill(lambda[iL+40],AH_energy*ah_mip2gev);

	    h_longi_profile_raw_ShowerInEE_layer->Fill(iL+1+40,AH_energy);
	    h_longi_profile_gev_ShowerInEE_layer->Fill(iL+1+40,AH_energy*ah_mip2gev);


	  }
	  else {  // MIPs in CE-E
	    h_longi_profile_MipsInEE->Fill(lambda[iL+40],AH_energy);
	    h_longi_profile_MipsInEE_gev->Fill(lambda[iL+40],AH_energy*ah_mip2gev);
	    h_longi_profile_MipsInEE_SS_ref->Fill(lambda[iL+40]-lambda[shower_start_index],AH_energy);

	    h_longi_profile_raw_MipsInEE_layer->Fill(iL+1+40,AH_energy);
	    h_longi_profile_gev_MipsInEE_layer->Fill(iL+1+40,AH_energy*ah_mip2gev);
	    h_longi_profile_raw_MipsInEE_SS_ref_layer->Fill(iL+1+40-shower_start_index,AH_energy);

	  }
	
	  h_longi_profile_raw[shower_start_index]->Fill(lambda[iL+40],AH_energy);
	  h_longi_profile_gev[shower_start_index]->Fill(lambda[iL+40],AH_detscale*w3);

	  h_longi_profile_raw_layer[shower_start_index]->Fill(iL+1+40,AH_energy);
	  h_longi_profile_gev_layer[shower_start_index]->Fill(iL+1+40,AH_detscale*w3);

	  h_longi_profile_raw_fraction[shower_start_index]->Fill(lambda[iL+40],AH_energy/total_energy);
	  h_longi_profile_gev_fraction[shower_start_index]->Fill(lambda[iL+40],AH_detscale*w3/total_energy_gev);

	  h_longi_profile_raw_layer_fraction[shower_start_index]->Fill(iL+1+40,AH_energy/total_energy);
	  h_longi_profile_gev_layer_fraction[shower_start_index]->Fill(iL+1+40,AH_detscale*w3/total_energy_gev);

	  h_longi_profile_inclusive->Fill(lambda[iL+40],AH_energy);
	  h_longi_profile_inclusive_frac->Fill(lambda[iL+40],AH_energy/total_energy);

	  h_longi_profile_raw_inclusive_layer->Fill(iL+1+40,AH_energy);
	  h_longi_profile_raw_inclusive_frac_layer->Fill(iL+1+40,AH_energy/total_energy);

	  h_Rechits_nrec_SS[shower_start_index][iL+40]->Fill( ahc_nHitsPerLayer->at(iL) );
	  h_Rechits_En_SS[shower_start_index][iL+40]->Fill(AH_energy);

	  if(shower_start_index+1 <= 7) {
	    h_longi_profile_raw_inc[0]->Fill(lambda[iL+40],AH_energy);
	    h_longi_profile_gev_inc[0]->Fill(lambda[iL+40],AH_detscale*w3);
	    h_longi_profile_raw_layer_inc[0]->Fill(iL+40+1,AH_energy);
	    h_longi_profile_gev_layer_inc[0]->Fill(iL+40+1,AH_detscale*w3);
	  }
	  else if(shower_start_index+1 <= 14) {
	    h_longi_profile_raw_inc[1]->Fill(lambda[iL+40],AH_energy);
	    h_longi_profile_gev_inc[1]->Fill(lambda[iL+40],AH_detscale*w3);
	    h_longi_profile_raw_layer_inc[1]->Fill(iL+40+1,AH_energy);
	    h_longi_profile_gev_layer_inc[1]->Fill(iL+40+1,AH_detscale*w3);
	  }
	  else if(shower_start_index+1 <= 21) {
	    h_longi_profile_raw_inc[2]->Fill(lambda[iL+40],AH_energy);
	    h_longi_profile_gev_inc[2]->Fill(lambda[iL+40],AH_detscale*w3);
	    h_longi_profile_raw_layer_inc[2]->Fill(iL+40+1,AH_energy);
	    h_longi_profile_gev_layer_inc[2]->Fill(iL+40+1,AH_detscale*w3);
	  }
	  else if(shower_start_index+1 <= 28) {
	    h_longi_profile_raw_inc[3]->Fill(lambda[iL+40],AH_energy);
	    h_longi_profile_gev_inc[3]->Fill(lambda[iL+40],AH_detscale*w3);
	    h_longi_profile_raw_layer_inc[3]->Fill(iL+40+1,AH_energy);
	    h_longi_profile_gev_layer_inc[3]->Fill(iL+40+1,AH_detscale*w3);
	  }
	  else {
	    h_longi_profile_raw_inc[shower_start_index-24]->Fill(lambda[iL+40],AH_energy);
	    h_longi_profile_gev_inc[shower_start_index-24]->Fill(lambda[iL+40],AH_detscale*w3);
	    h_longi_profile_raw_layer_inc[shower_start_index-24]->Fill(iL+40+1,AH_energy);
	    h_longi_profile_gev_layer_inc[shower_start_index-24]->Fill(iL+40+1,AH_detscale*w3);
	  }

	}

      } // E N D    L O N G I    S H O W E R     S H A P E S

      if(TRANSVERSE_PROFILE_HIST) {

	//////////////////////////////////////////////////
	/////////// Transverse shower profile        /////
	//////////////////////////////////////////////////
      

	//extraxt seed and cog coordinates for this HGCAL layer
      
	float seed_x = seed_xy[iL].first;  
	float seed_y = seed_xy[iL].second;

	float cog_x = rechit_cogX->at(iL);
	float cog_y = rechit_cogY->at(iL);
      
	if(seed_x > -999) { // if energy in this layer is non-zero
	  double trackx = track_x[iL];
	  double tracky = track_y[iL];
      
	  if(!strcmp(data,"data")) {  // do alignment for data events (not needed for sim)
	    trackx = -1*trackx;
	    tracky = -1*tracky;
      
	    std::pair<float, float> dxy_al = dxy_alignment(iL+1);
	    float dx_corr = dxy_al.first;
	    float dy_corr = dxy_al.second; 
	    seed_x -= dx_corr;
	    seed_y -= dy_corr;

	    cog_x -= dx_corr;
	    cog_y -= dy_corr;
	  
	  }

	  //calculate and fill differences b/w seed,cog AND trackImpact position
	  double dX = (seed_x - trackx);
	  double dY = (seed_y - tracky);
	  // cout<<"layer : seedX : seedY : dX : dY ::  "<<iL+1<<" : "<<seed_xy[iL].first<<" : "<<seed_xy[iL].second<<" : "<<dX<<" : "<<dY<<endl;
	  h_track_seed_diff_x[iL]->Fill(dX);
	  h_track_seed_diff_y[iL]->Fill(dY);
	  h_track_seed_diff_dR[iL]->Fill(deltaR(seed_x,seed_y,trackx,tracky));

	  dX = (cog_x - trackx);
	  dY = (cog_y - tracky);
	  h_track_cog_diff_x[iL]->Fill(dX);
	  h_track_cog_diff_y[iL]->Fill(dY);
	  h_track_cog_diff_dR[iL]->Fill(deltaR(cog_x,cog_y,trackx,tracky));
	
	
	}

	// layer-wise variables init for transverse-shower-profile
	double sum_dR0p56 = 0.0;
	double sum_dR1 = 0.0;
	double sum_dR2 = 0.0;
	double sum_dR3 = 0.0;
	double sum_dR5 = 0.0;
	double sum_dR8 = 0.0;
	double sum_dR10 = 0.0;
	double sum_dR12 = 0.0;
	double sum_dR15 = 0.0;
	double sum_dR18 = 0.0;
	double sum_dR20 = 0.0;
	double accumulated_energy = 0.0;
	TH1F* htemp = new TH1F("htemp","htemp",500,0,50);

	for(int i = 0; i < (int)rechit_collection[iL].size(); i++) { // loop over HGCAL rechit collection for this layer
	  RecHit* temp_rechit = rechit_collection[iL].at(i);
	  double	recx = (temp_rechit->getCellCoordinate_xyz())[0];
	  double	recy = (temp_rechit->getCellCoordinate_xyz())[1];
	  double trackx = track_x[iL];
	  double tracky = track_y[iL];

	  double trackx_ = 0.0;
	  double tracky_ = 0.0;

	  double recx_ = recx;
	  double recy_ = recy;
	
	  int en_chan = temp_rechit->getEncodedChannel();
	  if(!strcmp(data,"data")) { //alignment for data (not needed for sim)
	    trackx = -1*trackx;
	    tracky = -1*tracky;
	    // trackx = trackx;
	    // tracky = tracky;
	  
	    std::pair<float, float> dxy_al = dxy_alignment(iL+1);
	    float dx_corr = dxy_al.first;
	    float dy_corr = dxy_al.second; 
	    recx -= dx_corr;
	    recy -= dy_corr;

	    trackx_ = trackx + dx_corr;
	    tracky_ = tracky + dy_corr;
	  }
	  //recx & recy are corrected/aligned rechits
	  //recx_ and recy_ are uncorrected/un-aligned rechits
	  //double dR = deltaR(recx_,recy_,rechit_cogX->at(iL),rechit_cogY->at(iL)); // COG as the center
	  double dR = deltaR(recx,recy,trackx,tracky); // track impact as the center
	  double energy = temp_rechit->getEnergy();

	  //convert energy from mips to gev for pions
	  // if(pdgID == 11) {
	  //   if(iL < 28) energy = ee_mip2gev *  energy;
	  //   else        energy = fh_mip2gev *  energy;

	  // }
	  // else {
	  //   if(iL < 28) energy = ee_mip2gev * w1 * energy;
	  //   else        energy = fh_mip2gev * w2 * energy;
	  // }
	  // Fill energy-weighted distance 1D histograms
	  if(shower_start_index+1 <= 7 || shower_start_index+1 == 29 || shower_start_index+1 == 37) {
	    if(temp_rechit->getModuleNumber() == 38 && (en_chan == 16 || en_chan == 3004)) continue; // FH8 masking noisey channel


	    // h_transverse_distance[0][iL]->Fill(dR,energy);
	    // h_transverse_distance_inclusive[0]->Fill(dR,energy);
	    // h_transverse_distance_prof[0][iL]->Fill(dR,energy);
	    // h_transverse_distance_inclusive_prof[0]->Fill(dR,energy);
	    // h_transverse_distance_zoom[0][iL]->Fill(dR,energy);
	    // h_transverse_distance_inclusive_zoom[0]->Fill(dR,energy);
	    // if(energy > 100) {cout<<"Energy = "<<energy<<endl;return;}

	    //if(iL+1 == 14 || iL+1 == 32) htemp->Fill(dR,energy);
	    htemp->Fill(dR,energy);
	    if(iL+1 <= 28) { htemp_EE->Fill(dR,energy);}
	    else htemp_FH->Fill(dR,energy);
	   
	  }

	  
	  //dR = deltaR(recx_,recy_,rechit_cogX->at(iL),rechit_cogY->at(iL)); // COG as the center
	  dR = deltaR(recx,recy,trackx,tracky); // track impact as the center

	  // Calculate energies for different radii in THIS layer
	
	  if(dR < 0.56) sum_dR0p56 += temp_rechit->getEnergy();
	  // if(dR < 0.65) sum_dR0p56 += temp_rechit->getEnergy();
	  if(dR < 1.0) sum_dR1 += temp_rechit->getEnergy();
	  if(dR < 2.0) sum_dR2 += temp_rechit->getEnergy();
	  if(dR < 3.0) sum_dR3 += temp_rechit->getEnergy();
	  if(dR < 5.0) sum_dR5 += temp_rechit->getEnergy();
	  if(dR < 8.0) sum_dR8 += temp_rechit->getEnergy();
	  if(dR < 10.0) sum_dR10 += temp_rechit->getEnergy();
	  if(dR < 12.0) sum_dR12 += temp_rechit->getEnergy();
	  if(dR < 15.0) sum_dR15 += temp_rechit->getEnergy();
	  if(dR < 18.0) sum_dR18 += temp_rechit->getEnergy();

	  if(dR < 20.0) sum_dR20 += temp_rechit->getEnergy();

		
	} // E N D    O F   R E C H I T    C O L L E C T I O N    L O O P

	if(shower_start_index+1 <= 7 ) {
	  if(rechit_energyPerLayer->at(iL) != 0)  {
	    // cout<<"layer : htemp_integral : rechitEnergyLayer :: "<<iL+1<<" : "<<htemp->Integral() <<" : "<<rechit_energyPerLayer->at(iL)<<endl;
	    // return;
	    getProfile(h_transverse_prof_frac[0][iL],htemp,rechit_energyPerLayer->at(iL));

	  }
	}
	else if(shower_start_index+1 == 29) {
	  if(rechit_energyPerLayer->at(iL) != 0) 
	    getProfile(h_transverse_prof_frac[1][iL],htemp,rechit_energyPerLayer->at(iL));
	
	}
	else if(shower_start_index+1 == 37) {
	  if(rechit_energyPerLayer->at(iL) != 0) 
	    getProfile(h_transverse_prof_frac[2][iL],htemp,rechit_energyPerLayer->at(iL));
	
	}
	delete htemp;

	// assign summed up transverse variables inclusive in 7 EE layers
	if(iL+1 <= 7) {
	  sum_dRi_01_07[0] += sum_dR0p56;
	  sum_dRi_01_07[1] += sum_dR1;
	  sum_dRi_01_07[2] += sum_dR2;
	  sum_dRi_01_07[3] += sum_dR3;
	  sum_dRi_01_07[4] += sum_dR5;
	  sum_dRi_01_07[5] += sum_dR8;
	  sum_dRi_01_07[6] += sum_dR10;
	  sum_dRi_01_07[7] += sum_dR12;
	  sum_dRi_01_07[8] += sum_dR15;
	  sum_dRi_01_07[9] += sum_dR18;
	  sum_dRi_01_07[10] += sum_dR20;
	}
	else if(iL+1 <= 14) {
	  sum_dRi_08_14[0] += sum_dR0p56;
	  sum_dRi_08_14[1] += sum_dR1;
	  sum_dRi_08_14[2] += sum_dR2;
	  sum_dRi_08_14[3] += sum_dR3;
	  sum_dRi_08_14[4] += sum_dR5;
	  sum_dRi_08_14[5] += sum_dR8;
	  sum_dRi_08_14[6] += sum_dR10;
	  sum_dRi_08_14[7] += sum_dR12;
	  sum_dRi_08_14[8] += sum_dR15;
	  sum_dRi_08_14[9] += sum_dR18;
	  sum_dRi_08_14[10] += sum_dR20;
	}
	else if(iL+1 <= 21) {
	  sum_dRi_15_21[0] += sum_dR0p56;
	  sum_dRi_15_21[1] += sum_dR1;
	  sum_dRi_15_21[2] += sum_dR2;
	  sum_dRi_15_21[3] += sum_dR3;
	  sum_dRi_15_21[4] += sum_dR5;
	  sum_dRi_15_21[5] += sum_dR8;
	  sum_dRi_15_21[6] += sum_dR10;
	  sum_dRi_15_21[7] += sum_dR12;
	  sum_dRi_15_21[8] += sum_dR15;
	  sum_dRi_15_21[9] += sum_dR18;
	  sum_dRi_15_21[10] += sum_dR20;
	  
	}
	else if(iL+1 <= 28) {
	  sum_dRi_22_28[0] += sum_dR0p56;
	  sum_dRi_22_28[1] += sum_dR1;
	  sum_dRi_22_28[2] += sum_dR2;
	  sum_dRi_22_28[3] += sum_dR3;
	  sum_dRi_22_28[4] += sum_dR5;
	  sum_dRi_22_28[5] += sum_dR8;
	  sum_dRi_22_28[6] += sum_dR10;
	  sum_dRi_22_28[7] += sum_dR12;
	  sum_dRi_22_28[8] += sum_dR15;
	  sum_dRi_22_28[9] += sum_dR18;
	  sum_dRi_22_28[10] += sum_dR20;
	  
	}



	// w_t variable shall be used to convert corresponding energy from MIP to GeV
	double w_t = 0.0;
	if(iL < 28) w_t = w1 * ee_mip2gev;
	else w_t = w2 * fh_mip2gev;



      
	if(iL < 28) {
	  sum_dR1_EE += sum_dR1;
	  sum_dR2_EE += sum_dR2;
	  sum_dR3_EE += sum_dR3;
	  sum_dR5_EE += sum_dR5;
	  sum_dR8_EE += sum_dR8;
	  sum_dR10_EE += sum_dR10;
	  sum_dR12_EE += sum_dR12;
	  sum_dR15_EE += sum_dR15;
	  sum_dR18_EE += sum_dR18;
	  sum_dR20_EE += sum_dR20;
	
	}
	else {
	  sum_dR1_FH += sum_dR1;
	  sum_dR2_FH += sum_dR2;
	  sum_dR3_FH += sum_dR3;
	  sum_dR5_FH += sum_dR5;
	  sum_dR8_FH += sum_dR8;
	  sum_dR10_FH += sum_dR10;
	  sum_dR12_FH += sum_dR12;
	  sum_dR15_FH += sum_dR15;
	  sum_dR18_FH += sum_dR18;
	  sum_dR20_FH += sum_dR20;

	}

	/// For AHCAL
	if(iL < 39) {
	  // if(jentry < 10) {
	  //   cout<<"iL+1 : cogX_AH[iL] : cogY_AH[iL] : rechitEnergyLayer_AH[iL] :: "<<iL+1<<" : "<<cogX_AH[iL]<<" : "<<cogY_AH[iL]<<" : "<<rechitEnergyLayer_AH[iL]<<endl;
	  // }
	  if(rechitEnergyLayer_AH[iL] == 0.0) continue;
	  TH1F* htemp1 = new TH1F("htemp1","htemp1",500,0,50);
	  for(int i = 0; i < (int)rechit_AH_collection[iL].size(); i++) {
	    RecHit_AHCAL* temp_rechit = rechit_AH_collection[iL].at(i);
	    double recx = (temp_rechit->getCellCoordinate_xyz())[0];
	    double recy = (temp_rechit->getCellCoordinate_xyz())[1];
	    double trackx = track_x[iL];
	    double tracky = track_y[iL];

	    double recx_ = recx;
	    double recy_ = recy;


	    if(!strcmp(data,"data")) { //alignment for data (not needed for sim)
	      trackx = -1*trackx;
	      tracky = -1*tracky;
	      // trackx = trackx;
	      // tracky = tracky;
	      
	      std::pair<float, float> dxy_al = dxy_alignment(iL+1);
	      float dx_corr = dxy_al.first;
	      float dy_corr = dxy_al.second; 
	      recx -= dx_corr;
	      recy -= dy_corr;
	      
	    }

	    // double dR = deltaR(recx_,recy_,cogX_AH[iL],cogY_AH[iL]); // COG as the center
	    double dR = deltaR(recx,recy,trackx,tracky); // track impact as the center
	  
	    double energy = ah_mip2gev * w3 * temp_rechit->getEnergy();


	    if(shower_start_index+1 <= 7 || shower_start_index+1 == 29 || shower_start_index+1 == 37) {
	      htemp1->Fill(dR,temp_rechit->getEnergy());
	      htemp_AH->Fill(dR,temp_rechit->getEnergy());
	    } 
	    // if(shower_start_index+1 <= 7) {
	    //   h_transverse_distance[0][iL+40]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive[0]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive_compart[0][2]->Fill(dR,energy);
	    //   h_transverse_distance_prof[0][iL+40]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive_prof[0]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive_compart_prof[0][2]->Fill(dR,energy);
	    //   h_transverse_distance_zoom[0][iL+40]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive_zoom[0]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive_compart_zoom[0][2]->Fill(dR,energy);


	    // }
	    // else if(shower_start_index+1 <= 14) {
	    //   h_transverse_distance[1][iL+40]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive[1]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive_compart[1][2]->Fill(dR,energy);
	    //   h_transverse_distance_prof[1][iL+40]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive_prof[1]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive_compart_prof[1][2]->Fill(dR,energy);
	    //   h_transverse_distance_zoom[1][iL+40]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive_zoom[1]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive_compart_zoom[1][2]->Fill(dR,energy);


	    // }
	    // else if(shower_start_index+1 <= 21) {
	    //   h_transverse_distance[2][iL+40]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive[2]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive_compart[2][2]->Fill(dR,energy);
	    //   h_transverse_distance_prof[2][iL+40]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive_prof[2]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive_compart_prof[2][2]->Fill(dR,energy);
	    //   h_transverse_distance_zoom[2][iL+40]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive_zoom[2]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive_compart_zoom[2][2]->Fill(dR,energy);


	    // }
	    // else if(shower_start_index+1 <= 28) {
	    //   h_transverse_distance[3][iL+40]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive[3]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive_compart[3][2]->Fill(dR,energy);
	    //   h_transverse_distance_prof[3][iL+40]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive_prof[3]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive_compart_prof[3][2]->Fill(dR,energy);
	    //   h_transverse_distance_zoom[3][iL+40]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive_zoom[3]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive_compart_zoom[3][2]->Fill(dR,energy);


	    // }
	    // else {
	    //   h_transverse_distance[shower_start_index-24][iL+40]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive[shower_start_index-24]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive_compart[shower_start_index-24][2]->Fill(dR,energy);
	    //   h_transverse_distance_prof[shower_start_index-24][iL+40]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive_prof[shower_start_index-24]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive_compart_prof[shower_start_index-24][2]->Fill(dR,energy);
	    //   h_transverse_distance_zoom[shower_start_index-24][iL+40]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive_zoom[shower_start_index-24]->Fill(dR,energy);
	    //   h_transverse_distance_inclusive_compart_zoom[shower_start_index-24][2]->Fill(dR,energy);


	    // }
	  } // E N D    O F   R E C H I T    C O L L E C T I O N    L O O P for AHCAL
	  if(shower_start_index+1 <= 7 ) {
	    if(rechitEnergyLayer_AH[iL] != 0) 
	      getProfile(h_transverse_prof_frac[0][iL+40],htemp1,rechitEnergyLayer_AH[iL]);
	  }
	  else if(shower_start_index+1 == 29) {
	    if(rechitEnergyLayer_AH[iL] != 0) 
	      getProfile(h_transverse_prof_frac[1][iL+40],htemp1,rechitEnergyLayer_AH[iL]);
	    
	  }
	  else if(shower_start_index+1 == 37) {
	    if(rechitEnergyLayer_AH[iL] != 0) 
	      getProfile(h_transverse_prof_frac[2][iL+40],htemp1,rechitEnergyLayer_AH[iL]);
	    
	  }

	  delete htemp1;
	}
	//// End for AHCAL
	
      }
      
    }


    if(TRANSVERSE_PROFILE_HIST && !MIP && isInTrackWindow && shower_start_index+1 >= 3) {
      // if(MIP || !isInTrackWindow) continue; //to do: use these cleaning cuts on top
      // if(shower_start_index+1 <= 2) continue; // Rejection for pre-showering events

      if(shower_start_index+1 <= 7) {
	// if(rechitEnergySum_EE !=0) {
	//   // cout<<"jentry : htemp_EE_integral : Entries : rechitEnergyLayer_EE :: "<<jentry<<" : "<<htemp_EE->Integral() <<" : "<<htemp_EE->GetEntries()<<" : "<<rechitEnergySum_EE<<endl;
	//   // for(int i = 0 ; i < 28; i++) {
	//   //   cout<<"layer : layer_energy :: "<<i+1<<" : "<<rechit_energyPerLayer->at(i)<<endl;
	//   // }
	//   return;
	// }
      	if(rechitEnergySum_EE !=0) getProfile(h_transverse_prof_frac_Inc[0][0],htemp_EE,rechitEnergySum_EE);
      	if(rechitEnergySum_FH !=0) getProfile(h_transverse_prof_frac_Inc[0][1],htemp_FH,rechitEnergySum_FH);
      	if(rechitEnergySum_AH !=0) getProfile(h_transverse_prof_frac_Inc[0][2],htemp_AH,rechitEnergySum_AH);
      }
      else if(shower_start_index+1 == 29) {
      	if(rechitEnergySum_EE !=0) getProfile(h_transverse_prof_frac_Inc[1][0],htemp_EE,rechitEnergySum_EE);
      	if(rechitEnergySum_FH !=0) getProfile(h_transverse_prof_frac_Inc[1][1],htemp_FH,rechitEnergySum_FH);
      	if(rechitEnergySum_AH !=0) getProfile(h_transverse_prof_frac_Inc[1][2],htemp_AH,rechitEnergySum_AH);

      }
      else if(shower_start_index+1 == 37) {
      	if(rechitEnergySum_EE !=0) getProfile(h_transverse_prof_frac_Inc[2][0],htemp_EE,rechitEnergySum_EE);
      	if(rechitEnergySum_FH !=0) getProfile(h_transverse_prof_frac_Inc[2][1],htemp_FH,rechitEnergySum_FH);
      	if(rechitEnergySum_AH !=0) getProfile(h_transverse_prof_frac_Inc[2][2],htemp_AH,rechitEnergySum_AH);

      }
      
    }
    delete htemp_EE;
    delete htemp_FH;
    delete htemp_AH;
    // if(TRANSVERSE_PROFILE_HIST) {
    //   h_transverse_prof_EE[shower_start_index]->Fill(1.0,sum_dR1_EE);
    //   h_transverse_prof_EE[shower_start_index]->Fill(2.0,sum_dR2_EE);
    //   h_transverse_prof_EE[shower_start_index]->Fill(3.0,sum_dR3_EE);
    //   h_transverse_prof_EE[shower_start_index]->Fill(5.0,sum_dR5_EE);
    //   h_transverse_prof_EE[shower_start_index]->Fill(8.0,sum_dR8_EE);
    //   h_transverse_prof_EE[shower_start_index]->Fill(10.0,sum_dR10_EE);
    //   h_transverse_prof_EE[shower_start_index]->Fill(12.0,sum_dR12_EE);
    //   h_transverse_prof_EE[shower_start_index]->Fill(15.0,sum_dR15_EE);
    //   h_transverse_prof_EE[shower_start_index]->Fill(18.0,sum_dR18_EE);
    //   h_transverse_prof_EE[shower_start_index]->Fill(20.0,sum_dR20_EE);

    //   h_transverse_prof_FH[shower_start_index]->Fill(1.0,sum_dR1_FH);
    //   h_transverse_prof_FH[shower_start_index]->Fill(2.0,sum_dR2_FH);
    //   h_transverse_prof_FH[shower_start_index]->Fill(3.0,sum_dR3_FH);
    //   h_transverse_prof_FH[shower_start_index]->Fill(5.0,sum_dR5_FH);
    //   h_transverse_prof_FH[shower_start_index]->Fill(8.0,sum_dR8_FH);
    //   h_transverse_prof_FH[shower_start_index]->Fill(10.0,sum_dR10_FH);
    //   h_transverse_prof_FH[shower_start_index]->Fill(12.0,sum_dR12_FH);
    //   h_transverse_prof_FH[shower_start_index]->Fill(15.0,sum_dR15_FH);
    //   h_transverse_prof_FH[shower_start_index]->Fill(18.0,sum_dR18_FH);
    //   h_transverse_prof_FH[shower_start_index]->Fill(20.0,sum_dR20_FH);
    // }

    

    ////////////////////////////////////////////////////////////
    /////////    End of shower shape calculation    ///////////
    ///////////////////////////////////////////////////////////
    
    time_LongiCalc += sw_longi.RealTime();
    sw_longi.Stop();
    
    if(DEBUG) { cout << "SW >>>>>>END of shower shape loop; Took: "; sw_longi.Print();  }

    double isMuonLike = false;
    bool isRegion1 = false; //MIP like
    bool isRegion2 = false; //showering in CE-H
    bool isRegion3 = false;  //showering in CE-E
    bool isRegion4 = false; //showering in AHCAL

    TStopwatch sw_energy;
    sw_energy.Start();

    

    if(SHOWER_RECO_HIST) {
      ///////////////////////////////////////////////////////////
      //// E V E N T     C A T E G O R I Z A T I O N       /////
      ///////////////////////////////////////////////////////////


      
      if(!MIP && isInTrackWindow) { // to do: use these two cleaning cuts at the top !!!!
        
      	h_EE_fraction[shower_start_index]->Fill(EE_frac/100.0);
      	h_FH_fraction[shower_start_index]->Fill(FH_frac/100.0);
      	h_AH_fraction[shower_start_index]->Fill(AH_frac/100.0);
        
      	if(shower_start_index+1 <= 28) {
      	  h_ShowerInEE_EE_fraction->Fill(EE_frac/100.0);
      	  h_ShowerInEE_FH_fraction->Fill(FH_frac/100.0);
      	  h_ShowerInEE_AH_fraction->Fill(AH_frac/100.0);
      	}
      }

    
      /////// BASED ON SHOWER START
      //MUON VETO
      if(MuonVeto) {
	isRegion1 = true;
	isRegion2 = false;
	isRegion3 = false;
      }
      if(shower_start_index == -1) {  //MIP LIKE
	//if(shower_start_index < 2) {  //MIP LIKE
	isRegion1 = true;
	isRegion2 = false;
	isRegion3 = false;
      }
      else if(shower_start_index > 27) { //H hadrons
	isRegion1 = false;
	isRegion2 = true;
	isRegion3 = false;
      }
      else {  //EH hadron
	isRegion1 = false;
	isRegion2 = false;
	isRegion3 = true;
      }

      



      h_Nrechit_EE->Fill(Nrechit_EE);
      h_Nrechit_FH->Fill(Nrechit_FH);
      h_Nrechit_AH->Fill(ahc_nHits);

      h_Nrechit_low_EE->Fill(Nrechit_EE);
      h_Nrechit_low_FH->Fill(Nrechit_FH);
      h_Nrechit_low_AH->Fill(ahc_nHits);

    
      h_rechit_energy_raw_EE->Fill(rechitEnergySum_EE);
      h_rechit_energy_raw_FH->Fill(rechitEnergySum_FH);
      h_rechit_energy_raw_AH->Fill(rechitEnergySum_AH);
      h_rechit_energy_raw_all->Fill(rechitEnergySum_EE+rechitEnergySum_FH+rechitEnergySum_AH);


      h_rechit_energy_raw_low_EE->Fill(rechitEnergySum_EE);
      h_rechit_energy_raw_low_FH->Fill(rechitEnergySum_FH);
      h_rechit_energy_raw_low_AH->Fill(rechitEnergySum_AH);
      h_rechit_energy_raw_low_all->Fill(rechitEnergySum_EE+rechitEnergySum_FH+rechitEnergySum_AH);

      h_rechit_energy_raw_verylow_EE->Fill(rechitEnergySum_EE);
      h_rechit_energy_raw_verylow_FH->Fill(rechitEnergySum_FH);
      h_rechit_energy_raw_verylow_AH->Fill(rechitEnergySum_AH);
      h_rechit_energy_raw_verylow_all->Fill(rechitEnergySum_EE+rechitEnergySum_FH+rechitEnergySum_AH);

      h_rechit_energy_raw_mid_EE->Fill(rechitEnergySum_EE);
      h_rechit_energy_raw_mid_FH->Fill(rechitEnergySum_FH);
      h_rechit_energy_raw_mid_all->Fill(rechitEnergySum_EE+rechitEnergySum_FH+rechitEnergySum_AH);

      h_rechit_energy_raw_high_EE->Fill(rechitEnergySum_EE);
      h_rechit_energy_raw_high_FH->Fill(rechitEnergySum_FH);
      h_rechit_energy_raw_high_all->Fill(rechitEnergySum_EE+rechitEnergySum_FH+rechitEnergySum_AH);


      h_Nrechit_EE_vs_FH->Fill(Nrechit_EE,Nrechit_FH);




      // MIPs in CE-E and MIPs in CE-H
      if(isRegion1) {
	h_shower_start_reg1->Fill(shower_lambda_);
	region_1_classified_events++; 
      
	h_rechit_energy_raw_EE_MipsInEEFH->Fill(rechitEnergySum_EE);
	h_rechit_energy_raw_FH_MipsInEEFH->Fill(rechitEnergySum_FH);
	h_rechit_energy_raw_EE_vs_FH_MipsInEEFH->Fill(rechitEnergySum_EE,rechitEnergySum_FH);
	h_rechit_energy_raw_FH_vs_AH_MipsInEEFH->Fill(rechitEnergySum_FH,rechitEnergySum_AH);
	h_rechit_energy_raw_AH_MipsInEEFH->Fill(rechitEnergySum_AH);
	h_rechit_energy_raw_EE_MipsInEEFH_extended->Fill(rechitEnergySum_EE);
	h_rechit_energy_raw_FH_MipsInEEFH_extended->Fill(rechitEnergySum_FH);
	h_rechit_energy_raw_AH_MipsInEEFH_extended->Fill(rechitEnergySum_AH);
	h_rechit_energy_raw_AH_MipsInEEFH_extended_v1->Fill(rechitEnergySum_AH);

	h_rechit_energy_raw_all_MipsInEEFH->Fill(rechitEnergySum_EE+rechitEnergySum_FH);
      
	h_rechit_energy_raw_EE_MipsInFH_extended->Fill(rechitEnergySum_EE);
	h_rechit_energy_raw_FH_MipsInFH_extended->Fill(rechitEnergySum_FH);
	h_rechit_energy_raw_AH_MipsInFH_extended->Fill(rechitEnergySum_AH);
	h_rechit_energy_raw_AH_MipsInFH_extended_v1->Fill(rechitEnergySum_AH);
      
      
      
      }
    
      // MIPs in CE-E
      else if(isRegion2) {
	h_shower_start_reg2->Fill(shower_lambda_);
	region_2_classified_events++; 
      

	if(Nrechit_EE == 0) EE_zeroHit_R2++;
      
	bool trackwindow = isInTrackWindow;
	bool SS1_reject = (shower_start_index+1 > 1);
	bool SS2_reject = (shower_start_index+1 > 2);
    

	//////////////////////////////////////////////////////////////////////////////////////
	///////////////   Based on Electron/pion scale in CE-E & CE-H-Si ////////////////////
	/////////////////////////////////////////////////////////////////////////////////////



	float total_energy = (rechitEnergySum_EE/EE_scale) + (rechitEnergySum_FH+ (alpha_ * rechitEnergySum_AH))/FH_AH_scale;
	float total_energy_withoutAH = (rechitEnergySum_EE/EE_scale) + (rechitEnergySum_FH)/FH_AH_scale;
      
      

	/////////////////////////////////////////////////////////
	/////     H hadrons ; Chi2 matrix initialzation     /////
	/////////////////////////////////////////////////////////
	double EE_detscale = (rechitEnergySum_EE/EE_scale);
	double FH_detscale = rechitEnergySum_FH/FH_AH_scale;
	double AH_detscale = (alpha_ * rechitEnergySum_AH)/FH_AH_scale;

	//double O = 0.4;
	double O = EE_detscale;

	float scale_ = scaling_factor_H[(int)beamEnergy];

	float w1 = getChi2Weights_H(beamEnergy).at(0);
	float w2 = getChi2Weights_H(beamEnergy).at(1);
	float w3 = getChi2Weights_H(beamEnergy).at(2);

	// float w1 = f_H_w1->Eval(beamEnergy);
	// float w2 = f_H_w2->Eval(beamEnergy);
	// float w3 = f_H_w3->Eval(beamEnergy);

	float recoE_w1 = f_H_w1->Eval(total_energy);
	float recoE_w2 = f_H_w2->Eval(total_energy);
	float recoE_w3 = f_H_w3->Eval(total_energy);

	float chi2_energy_recoE = recoE_w1*1.0 + recoE_w2*FH_detscale + recoE_w3*AH_detscale + O;
      
	float chi2_energy = w1*1.0 + w2*FH_detscale + w3*AH_detscale + O;
	float chi2_energy_withoutAHCAL = w1*1.0 + w2*FH_detscale + O;


      

	/////////////////////////////////////////////////////
	//////////////// selection cut check  ////////////////
	/////////////////////////////////////////////////////

	//here, there should be no checks!! make it consistent with the latest energy measurement code
	h_baseline_FH->Fill(total_energy);
	h_baseline_chi2_FH->Fill(chi2_energy);

	if(trackwindow)  {
	  h_trackwindow_FH->Fill(total_energy);
	  h_trackwindow_chi2_FH->Fill(chi2_energy);
	  h_trackwindow_chi2_FH_recoE->Fill(chi2_energy_recoE);
	
	}
	if(SS1_reject) {
	  h_SS1_reject_FH->Fill(total_energy);
	}
	if(SS2_reject) {
	  h_SS2_reject_FH->Fill(total_energy);
	}
	if(SS1_reject && trackwindow) {
	  h_trackwindow_SS1_reject_FH->Fill(total_energy);
	}
	if(SS2_reject && trackwindow) {
	  h_trackwindow_SS2_reject_FH->Fill(total_energy);

	  h_trackwindow_SS2_reject_chi2_FH->Fill(chi2_energy);
	  h_trackwindow_SS2_reject_chi2_FH_recoE->Fill(chi2_energy_recoE);

	  h_trackwindow_SS2_reject_chi2_all->Fill(chi2_energy);
	  h_trackwindow_SS2_reject_chi2_all_recoE->Fill(chi2_energy_recoE);

	
	  h_rechit_energy_raw_FH_MipsInEE_extended->Fill(rechitEnergySum_FH);
	  h_rechit_energy_raw_FH_MipsInEE_extended_v1->Fill(rechitEnergySum_FH);
	  h_rechit_energy_raw_AH_MipsInEE_extended->Fill(rechitEnergySum_AH);
	  h_rechit_energy_raw_AH_MipsInEE_extended_v1->Fill(rechitEnergySum_AH);
	
	  h_rechit_energy_raw_EE_vs_FH_MipsInEE->Fill(rechitEnergySum_EE,rechitEnergySum_FH);
	  h_rechit_energy_raw_FH_vs_AH_MipsInEE->Fill(rechitEnergySum_FH,rechitEnergySum_AH);
	  h_rechit_energy_raw_EE_MipsInEE->Fill(rechitEnergySum_EE);
	  h_rechit_energy_raw_EE_MipsInEE_extended->Fill(rechitEnergySum_EE);
	  h_rechit_energy_raw_FH_MipsInEE->Fill(rechitEnergySum_FH);
	  h_rechit_energy_raw_AH_MipsInEE->Fill(rechitEnergySum_AH);
      
	  h_rechit_energy_raw_all_MipsInEE->Fill(rechitEnergySum_EE+rechitEnergySum_FH);


	  h_EE_MIPsinEE->Fill(rechitEnergySum_EE);
	  h_FH_MIPsinEE->Fill(rechitEnergySum_FH);
	  h_AH_MIPsinEE->Fill(rechitEnergySum_AH);
	  h_all_MIPsinEE->Fill(rechitEnergySum_EE+rechitEnergySum_FH+rechitEnergySum_AH);
	
	  h_MIPsinEE_elecpi_scale->Fill(total_energy);
	  h_MIPsinEE_elecpi_scale_withoutAH->Fill(total_energy_withoutAH);


	  h_EE_inclusive->Fill(rechitEnergySum_EE);
	  h_FH_inclusive->Fill(rechitEnergySum_FH);
	  h_AH_inclusive->Fill(rechitEnergySum_AH);
	  h_all_inclusive->Fill(rechitEnergySum_EE+rechitEnergySum_FH+rechitEnergySum_AH);
	
	  h_inclusive_elecpi_scale->Fill(total_energy);


	
	}

	/////////////////////////////////////////////////////
	/////////////////////////////////////////////////////


      }

      else if(isRegion3) {  //Showering in CE-E
	h_shower_start_reg3->Fill(shower_lambda_);
	region_3_classified_events++; 
      

	bool trackwindow = isInTrackWindow;
	bool SS1_reject = (shower_start_index+1 > 1);
	bool SS2_reject = (shower_start_index+1 > 2);

	//////////////////////////////////////////////////////////////////////////////////////
	///////////////   Based on Electron/pion scale in CE-E & CE-H-Si //////////////////////
	//////////////////////////////////////////////////////////////////////////////////////

	float total_energy = (rechitEnergySum_EE/EE_scale) + (rechitEnergySum_FH+ (alpha_ * rechitEnergySum_AH))/FH_AH_scale;
	float total_energy_withoutAH = (rechitEnergySum_EE/EE_scale) + (rechitEnergySum_FH)/FH_AH_scale;
      


	///////////////////////////////////////////////////////////
	/////     EH hadrons ; Chi2 matrix initialzation     /////
	//////////////////////////////////////////////////////////
      
	double EE_detscale = (rechitEnergySum_EE/EE_scale);
	double FH_detscale = (rechitEnergySum_FH/FH_AH_scale);
	double AH_detscale = (alpha_*rechitEnergySum_AH)/FH_AH_scale;

	float scale_ = scaling_factor_EH[(int)beamEnergy];
      
	float w1 = getChi2Weights_EH(beamEnergy).at(0);
	float w2 = getChi2Weights_EH(beamEnergy).at(1);
	float w3 = getChi2Weights_EH(beamEnergy).at(2);
      
	// float w1 = f_EH_w1->Eval(beamEnergy);
	// float w2 = f_EH_w2->Eval(beamEnergy);
	// float w3 = f_EH_w3->Eval(beamEnergy);

	float recoE_w1 = f_EH_w1->Eval(total_energy);
	float recoE_w2 = f_EH_w2->Eval(total_energy);
	float recoE_w3 = f_EH_w3->Eval(total_energy);

	if(DEBUG) { cout<<"DEBUG :: beamE : w1 : w2 : w3 : total_E : Rw1 : Rw2 : Rw3 :: "<<beamEnergy<< " : " <<" : "<<w1 <<" : "<<w2<<" : "<<w3<<" : "<<total_energy<<" : "<<recoE_w1<< " : " << recoE_w2 << " : "<< recoE_w3 <<endl;}
      
	float chi2_energy_recoE = recoE_w1*EE_detscale + recoE_w2*FH_detscale + recoE_w3*AH_detscale;

	double O = 0.0;

	float chi2_energy = w1*EE_detscale + w2*FH_detscale + w3*AH_detscale + O;
	float chi2_energy_withoutAHCAL = w1*EE_detscale + w2*FH_detscale + O;



	/////////////////////////////////////////////////////
	//////////////// selection cut check  ////////////////
	/////////////////////////////////////////////////////

      
	h_baseline->Fill(total_energy);
	h_baseline_chi2_EE->Fill(chi2_energy);

	if(trackwindow)  {
	  h_trackwindow->Fill(total_energy);
	  h_trackwindow_chi2_EE->Fill(chi2_energy);
	
	}
	if(SS1_reject) {
	  h_SS1_reject->Fill(total_energy);
	  h_SS1_reject_chi2_EE->Fill(chi2_energy);
	}
	if(SS2_reject) {
	  h_SS2_reject->Fill(total_energy);
	  h_SS2_reject_chi2_EE->Fill(chi2_energy);
	}
	if(SS1_reject && trackwindow) {
	  h_trackwindow_SS1_reject->Fill(total_energy);
	  h_trackwindow_SS1_reject_chi2_EE->Fill(chi2_energy);
	}
	if(SS2_reject && trackwindow) {
	  h_trackwindow_SS2_reject->Fill(total_energy);

	  h_trackwindow_SS2_reject_chi2_EE->Fill(chi2_energy);
	  h_trackwindow_SS2_reject_chi2_EE_recoE->Fill(chi2_energy_recoE);
	  h_trackwindow_SS2_reject_chi2_all->Fill(chi2_energy);
	  h_trackwindow_SS2_reject_chi2_all_recoE->Fill(chi2_energy_recoE);

	
	  h_EE_showerinEE->Fill(rechitEnergySum_EE);
	  h_FH_showerinEE->Fill(rechitEnergySum_FH);
	  h_AH_showerinEE->Fill(rechitEnergySum_AH);
	  h_all_showerinEE->Fill(rechitEnergySum_EE+rechitEnergySum_FH+rechitEnergySum_AH);

	  h_showerinEE_elecpi_scale->Fill(total_energy);
	  h_showerinEE_elecpi_scale_withoutAH->Fill(total_energy_withoutAH);


	  h_EE_inclusive->Fill(rechitEnergySum_EE);
	  h_FH_inclusive->Fill(rechitEnergySum_FH);
	  h_AH_inclusive->Fill(rechitEnergySum_AH);
	  h_all_inclusive->Fill(rechitEnergySum_EE+rechitEnergySum_FH+rechitEnergySum_AH);
	
	  h_inclusive_elecpi_scale->Fill(total_energy);

	  h_rechit_energy_raw_EE_ShowerInEE->Fill(rechitEnergySum_EE);
	  h_rechit_energy_raw_FH_ShowerInEE->Fill(rechitEnergySum_FH);
	  h_rechit_energy_raw_EE_vs_FH_ShowerInEE->Fill(rechitEnergySum_EE,rechitEnergySum_FH);
	  h_rechit_energy_raw_FH_vs_AH_ShowerInEE->Fill(rechitEnergySum_FH,rechitEnergySum_AH);
	  h_rechit_energy_raw_AH_ShowerInEE->Fill(rechitEnergySum_AH);
	  h_rechit_energy_raw_all_ShowerInEE->Fill(rechitEnergySum_EE+rechitEnergySum_FH+rechitEnergySum_AH);

	}

	/////////////////////////////////////////////////////
	/////////////////////////////////////////////////////


      }

      //Showering in AHCAL ????
      if(isRegion4) {
	region_4_classified_events++;
      }
    }
    time_EnMeasureCalc += sw_energy.RealTime();
    sw_energy.Stop();
    
    if(DEBUG) { cout << "SW >>>>>>END of shower energy loop; Took: "; sw_energy.Print();  }

    if(Nrechit_EE == 0 && false) {
      cout<<" **** Nrechit_EE is zero for jentry = "<<jentry<<endl;
      cout<<" isRegion1 "<<isRegion1
	  <<" isRegion2 "<<isRegion2
	  <<" isRegion3 "<<isRegion3
	  <<" isRegion4 "<<isRegion4<<endl;
      cout<<" isGoodTrack : isInTrackWindow :: "<<isGoodTrack<<" : "<<isInTrackWindow<<endl;
      cout<<"*************"<<endl;
    }
    
    // if(MuonVeto) continue;
    // event_count[3]++;
    

    if(DEBUG) { cout << "SW >>>>>>END of event, Took: "; sw.Print(); }
    if(DEBUG) cout<<"DEBUG: End of Event = "<<jentry+1<<endl;
    if(DEBUG) cout<<"DEBUG: ****************** "<<endl;

    if(DEBUG && jentry > 5) break;
    
    // if(jentry > 1000) break;
    
  }
  //////////////////////////////////////////////
  //    E N D   O F   E N T R Y    L O O P    /        
  /////////////////////////////////////////////
  if(E_beam < 0) {
    cout<<"E_beam negative!!!"<<endl;
    return;
  }


  cout<<"***** Timing budget summmary ******"<<endl;
  cout << " Central clocktime; Took: "; sw.Print("m"); 
  sw.Stop();

  
  cout<<"Total calulated time(sec) : "<<(time_hgc+time_ahc+time_CompartCalc+time_LongiCalc+time_EnMeasureCalc)<<endl;
  cout<<"HGC took = "<< time_hgc <<endl;
  cout<<"AHC took = "<< time_ahc <<endl;
  cout<<"Compart Calc took = "<< time_CompartCalc <<endl;
  cout<<"Longi Calc took = "<< time_LongiCalc <<endl;
  cout<<"Energy measured took = "<< time_EnMeasureCalc <<endl;
  cout<<"***************************"<<endl<<endl;
  //cout<<"Got Out, total events = "<<jentry<<" ; selected events = "<<selected_events<<endl;

  cout<<"***** Event summary ******"<<endl;
  cout<<"Total events = "<<event_count[0]<<endl;
  cout<<"Events with good track = "<<event_count[1]<<endl;
  cout<<"Events with Nrechits in FH[8][3] < 80 and FH[8][4] < 80 = "<<event_count[2]<<endl;
  cout<<" Events passing Muon veto = "<<event_count[3]<<endl;
  cout<<"***************************"<<endl<<endl;
  

  Long64_t total_events_temp = (region_1_classified_events+region_2_classified_events+region_3_classified_events+region_4_classified_events+non_classified_events);
  cout<<"Events with zero hits in EE = "<<EE_zeroHit<<endl;
  cout<<"Events with zero hits in FH = "<<FH_zeroHit<<endl;
  cout<<"Events with zero hits in AHCAL = "<<ahc_zeroHit<<endl;
  cout<<"MIP like events = "<<((float)region_1_classified_events*100.0)/total_events_temp<<"%"<<endl;
  cout<<"shower start in EE = "<<((float)region_3_classified_events*100.0)/total_events_temp<<"%"<<endl;
  cout<<"shower start in FH = "<<((float)region_2_classified_events*100.0)/total_events_temp<<"%"<<endl;
  cout<<"shower start in AH = "<<((float)region_4_classified_events*100.0)/total_events_temp<<"%"<<endl;
  cout<<"Non-classified events = "<<((float)non_classified_events*100.0)/total_events_temp<<"%"<<endl;
  
  cout<<"Number of events with (Nrechit_EE == 0) = "<<EE_zeroHit_R2<<endl;
  cout<<"Sum = "<<total_events_temp<<endl;

  cout<<endl<<endl;

  cout<<"Total number of events = "<<debug_count<<endl;
}

