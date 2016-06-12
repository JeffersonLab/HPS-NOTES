#include <TH1D.h>
#include <TH2D.h>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TStyle.h>
#include <TChain.h>
#include <TSystem.h>
#include <TCanvas.h>
#include <TLorentzVector.h>

#include <HpsEvent.h>
#include <GblTrack.h>
#include <SvtTrack.h>
#include <MCParticle.h>

using namespace std;

int main()
{
  
  const double radian = TMath::RadToDeg();
  const int pair = 2;
  
  const HpsEvent::collection_t fs_part_type = HpsEvent::FINAL_STATE_PARTICLES;  // Collection name for FINAL_STATE_PARTICLES
  const HpsEvent::collection_t uc_part_type = HpsEvent::UC_VTX_PARTICLES;  // Collection name for unconstrained particles


  TChain *tr1 = new TChain("HPS_Event", "HPS_Event");

  //tr1->Add("/path/to/file");    // Add files to the chain
//tr1->Add("/cache/mss/hallb/hps/production/dst/wab-beam-tri/1pt05/wabv1-egsv3-triv2-g4v1_s2d6_HPS-EngRun2015-Nominal-v1_3.4.0-20150710_singles1_*.root"); // WAB-BEAM-TRI
  //tr1->Add("/u/home/byale/mc/egsv3-triv2-g4v1_s2d6_HPS-EngRun2015-Nominal-v1_3.4.0-20150710_pairs1_18.root");
tr1->Add("/work/hallb/hps/byale/moller/1pt05/molv1_s2d6_HPS-EngRun2015-Nominal-v1_3.4.0-20150710_singles1_18.root"); // PURE MOLLERS
//tr1->Add("/work/hallb/hps/data/engrun2015/pass1/dst/hps_005772.*_dst_R3321.root"); // DATA

  HpsEvent *ev1 = new HpsEvent(); // pointer to HPS Events, all tracks are accessed from HPSEvent
  HpsParticle *fs_part, *uc_part, *fs_part_upper, *dau_particles, *dau_part1, *dau_part2;
  EcalCluster *ec_clust, *ec_clust_upper, *dau_clust1, *dau_clust2;
  int charge1, charge2;

  tr1->SetBranchAddress("Event", &ev1);

  TLorentzVector L_em, L_ep;      // Lorentz vectors, that to be assigned to e- and e+
  
  TFile *file_out = new TFile("moller.root", "Recreate"); // Histograms, and other necessary objects will be saved in this file

  TH1D *h_n_fs_particles = new TH1D("h_n_fs_particles", "", 30, 0., 3);
  TH1D *h_n_uc_particles = new TH1D("h_n_uc_particles", "", 30, 0., 3);
  TH1D *h_fs_Minv1 = new TH1D("h_fs_Minv1", "", 200, 0., .3);
  TH1D *h_uc_Minv1 = new TH1D("h_uc_Minv1", "", 200, 0., .3);
  
  TH2D *moller_thetaE = new TH2D("moller_thetaE","Moller Theta vs. E (singles1 recon)",100,0.0,1.056,100,0.0,0.25);
  TH1D *moller_E = new TH1D("moller_E", "", 200, 0., 1.06);
  TH1D *moller_Theta = new TH1D("moller_Theta", "", 200, 0.0, 0.25);
  TH1D *moller_ESum = new TH1D("moller_E1+E2", "", 200, 0., 1.06);

  TH2D *electron_thetaE = new TH2D("electron_thetaE","Electron Theta vs. E (singles1 recon)",100,0.0,1.056,100,0.0,0.25);
  TH1D *fs_electron_E = new TH1D("fs_electron_E", "", 200, 0., 1.06);
  TH1D *fs_electron_Theta = new TH1D("fs_electron_Theta", "", 200, 0.0, 0.25);


  int nev = tr1->GetEntries();
  cout<<"Entries in the Tree = "<<nev<<endl;

  int n_neg1 = 0;
  int n_neg2 = 0;
  double E,theta;
  char outputname[200];

  // Loop over events
  for( int i = 0; i < nev; i++ )
    {

      if( i%50000 == 0 )
	{
	  cout.flush()<<"Processed "<<i<<" events:          about "<<(100.*double(i)/double(nev))<<" % \r";
	}

      // Read the current (i-th) entry from the tree
      tr1->GetEntry(i);
      
      // Number of tracks in the collection "FINAL_STATE_PARTICLES"
      int n_fs_particles = ev1->getNumberOfParticles(fs_part_type);
      
      // Fill the histogram with the number of final_state_particles
      h_n_fs_particles->Fill(n_fs_particles);


      n_neg1 = 0; 
      n_neg2 = 0;
      
      for( int j = 0; j < n_fs_particles; j++ )  // Now loop over particles in the collection FINAL_STATE_PARTICLES
	{
	  fs_part = ev1->getParticle(fs_part_type, j); // get the j-th particle from that collection
	  
	  TRefArray *ec_clusters = fs_part->getClusters();   // Get ECcal clusters assoicated with that particle
	  int n_clust = ec_clusters->GetEntries();           // NUmber of clusters (Should always be 0 or 1)
	  
	  TRefArray *svt_tracks = fs_part->getTracks();      // Get SVT tracks associated with this particle
	  int n_tracks = svt_tracks->GetEntries();           // Should always be (0 or 1)
	  
	  if( n_clust == 1 && n_tracks == 1 ) // neglect events with no tracks and with no clusters (e- and e+ should always produce hits on SVT and ECal)
	    {
	      
	      ec_clust = (EcalCluster*)ec_clusters->At(0);
	      vector<double>  mom = fs_part->getMomentum(); double px = mom[0]; double py = mom[1]; double pz = mom[2];
	      double P = sqrt( px*px + py*py + pz*pz );
	      double E = ec_clust->getEnergy();
              double theta = atan2(sqrt(px*px + py*py),pz*pz);
              double time = ec_clust->getSeed()->getTime();

	      int q = fs_part->getCharge();

              // Select electrons
	      if( q < 0 )
		{
		  L_em.SetPxPyPzE(px, py, pz, P);  // Just Assume that negative particle is electron
		  n_neg1 = n_neg1 + 1;
                  electron_thetaE->Fill(E,theta);
                  fs_electron_E->Fill(E);
                  fs_electron_Theta->Fill(theta);

              // 
                  //if(E>0.47-0.08 && E<0.47+0.08 && theta>0.0875-0.03 && theta<0.0875+0.03){
                  if(py<0 && E>=0.2 && E<=0.7 && theta>=0.01 && theta<=0.2){
                  //if(py<0)
//                    { 
                       vector<double> vtx_lower = fs_part->getVertexPosition();
                       for( int k = 0; k < n_fs_particles; k++ )
                         {                          
                           fs_part_upper = ev1->getParticle(fs_part_type, k);
                           TRefArray *ec_clusters_upper = fs_part_upper->getClusters();
                           int n_clust_upper = ec_clusters_upper->GetEntries();

                           TRefArray *svt_tracks_upper = fs_part_upper->getTracks();      // Get SVT tracks associated with this particle
                           int n_tracks_upper = svt_tracks_upper->GetEntries();

                           if( n_clust_upper == 1 && n_tracks_upper == 1 ) // neglect events with no tracks and with no clusters (e- and e+ should always produce hits on SVT and ECal)
            		     {
                           	ec_clust_upper = (EcalCluster*)ec_clusters_upper->At(0);

                           	vector<double> mom_upper = fs_part_upper->getMomentum(); double px_upper = mom_upper[0]; double py_upper = mom_upper[1]; double pz_upper = mom_upper[2];
                           	double P_upper = sqrt( px_upper*px_upper + py_upper*py_upper + pz_upper*pz_upper );
                           	double E_upper = ec_clust_upper->getEnergy();
                           	double theta_upper = atan2(sqrt(px_upper*px_upper + py_upper*py_upper),pz_upper*pz_upper);
                                double time_upper = ec_clust_upper->getSeed()->getTime();

                           	vector<double> vtx_upper = fs_part_upper->getVertexPosition();

                                // MOLLER CUTS
                                if(py_upper>=0 && vtx_upper == vtx_lower && E_upper>=0.2 && E_upper<=0.7 && theta_upper>=0.01 && theta_upper<=0.2 && E+E_upper >= 0.82 && abs(time-time_upper)<=0.3) // within 1% from before pass these cuts
                           	//if(py_upper>0 && vtx_upper == vtx_lower && abs(time-time_upper)<=0.3) // within 1% (89%) from below pass this timecut (300 ps), essentially all
                           	//if(py_upper>0 && vtx_upper == vtx_lower)// 89% pass this cut => 100% the mollers from the cut below have the same vertex
                           	//if(py_upper>0) // 89% of mollers pass this cut (mollers on opposite sides of the ECal) (This is a requirement for pairs1 anyway) Generator only 93% efficient anyway, so within 4% of that.
                             	{
                               		moller_thetaE->Fill(E,theta); moller_thetaE->Fill(E_upper,theta_upper);
                               		moller_E->Fill(E); moller_E->Fill(E_upper);
                               		moller_Theta->Fill(theta); moller_Theta->Fill(theta_upper);
                                        moller_ESum->Fill(E+E_upper);
                           	}
             		    }		
                          
                        }


  
                      // moller_thetaE->Fill(E,theta);
                      // moller_E->Fill(E);
                      // moller_Theta->Fill(theta);
//                  } // if(py<0)
                  }  // Moller cuts

		}
	     // else if( q > 0 )
	//	{
		  //L_ep.SetPxPyPzE(px, py, pz, P);
		  //n_pos = n_pos + 1;
		  // Do something with positrons;
	//	}
		
	    }
	  
	}
      // cout<<"Nneg = "<<n_neg<<"      n_pos = "<<n_pos<<endl;
//      if( n_neg1 == 1 && n_neg2 == 1 )  // Events with 1 negative and 1 positive tracks
//	{
//	  TLorentzVector L_emep = L_em + L_ep;
//	  double Minv = L_emep.M();
	  //h_fs_Minv1->Fill(E);
          //moller_thetaE->Fill(E,theta);
//	}
      
      

//      int n_uc_particles = ev1->getNumberOfParticles(uc_part_type);   // Number of tracks in the collection "UC_VTX_PARTICLES"
      // These are particles which actually contain 2 sub particles, i.e. two tracks were fitted to have a same vertex
      
//      h_n_uc_particles->Fill(n_uc_particles);
      

  // EVERYTHING YOU WANT TO DO WITH UC VERTEXED PARTICLES GOES HERE!!!!!!!!!! (MOLLERS ARE PROBABLY NOT VERTEXED THOUGH!!!)
//      for( int jj = 0; jj < n_uc_particles; jj++ )
//	{
//	  uc_part = ev1->getParticle(uc_part_type, jj); // get the jj-th particle from the UC_VTX_PARTICLES collection
//	  double m_uc_part = uc_part->Mass();   // Since this is a parrent of two particles (neg and pos), we can directly get the mass of it, which is the
//	                                        // invariant mass of it's child particles
//	  h_uc_Minv1->Fill(m_uc_part);
//                  
//          // Get the daughter particles from the uc_part
//          TRefArray *dau_particles = uc_part->getParticles();
//          int n_dau = dau_particles->GetEntries();
//          dau_part1 = (HpsParticle*)dau_particles->At(0);
//          dau_part2 = (HpsParticle*)dau_particles->At(1);
//
//          TRefArray *dau_clusters1 = dau_part1->getClusters();   // Get ECcal clusters assoicated with that particle
//          TRefArray *dau_clusters2 = dau_part2->getClusters();
//          //int n_dau_clust = dau_clusters->GetEntries();           // NUmber of clusters (Should always be 0 or 1)
//
//          TRefArray *dau_tracks1 = dau_part1->getTracks();      // Get SVT tracks associated with this particle
//          TRefArray *dau_tracks2 = dau_part2->getTracks();
//          //int n_dau_tracks = dau_tracks->GetEntries();           // Should always be (0 or 1)
//
//          //if( n_dau_clust == 1 && n_tracks == 1 ) // neglect events with no tracks and with no clusters (e- and e+ should always produce hits on SVT and ECal)
//            //{
//                   
//
//          dau_clust1 = (EcalCluster*)dau_clusters1->At(0);
//          vector<double>  mom1 = dau_part1->getMomentum(); double px1 = mom1[0]; double py1 = mom1[1]; double pz1 = mom1[2];
//          double P1 = sqrt( px1*px1 + py1*py1 + pz1*pz1 );
//          double E1 = dau_clust1->getEnergy();
//          double theta1 = atan2(sqrt(px1*px1 + py1*py1),pz1*pz1);
//
//
//          dau_clust2 = (EcalCluster*)dau_clusters2->At(0);
//          vector<double>  mom2 = dau_part2->getMomentum(); double px2 = mom2[0]; double py2 = mom2[1]; double pz2 = mom2[2];
//          double P2 = sqrt( px2*px2 + py2*py2 + pz2*pz2 );
//          double E2 = dau_clust2->getEnergy();
//          double theta2 = atan2(sqrt(px2*px2 + py2*py2),pz2*pz2);
//
//
//          // Get the charges of the daughter particles
//          charge1 = dau_part1->getCharge();
//          charge2 = dau_part2->getCharge();
//          //TRefArray *dau_particles = uc_part->getParticles();
//          
//          // Select vertexed particles with 2 negative daughters
//          //if(n_dau == 2 && charge1<0 && charge2<0){
//
//              // moller_thetaE->Fill(E1,theta1);
//              // moller_thetaE->Fill(E2,theta2);
//
//          // }
//
//
//
//	}
//     

    }
      
      TCanvas *c1 = new TCanvas("c1", "c1",10,10,1000,750);
      c1->SetLogz();

      h_n_fs_particles->Write();
      h_fs_Minv1->Write();
      h_n_uc_particles->Write();
      h_uc_Minv1->Write();
      electron_thetaE->Write();
      moller_thetaE->Write();      

      electron_thetaE->Draw("COLZ");
      strcpy(outputname,"electron_ThetaE");
      c1->SaveAs(strcat(outputname,"_electron_thetaE_singles1.png"));      

      fs_electron_E->Draw("COLZ");
      strcpy(outputname,"electron_E");
      c1->SaveAs(strcat(outputname,"_electron_E_singles1.png"));

      fs_electron_Theta->Draw("COLZ");
      strcpy(outputname,"electron_Theta");
      c1->SaveAs(strcat(outputname,"_electron_Theta_singles1.png"));

      moller_thetaE->Draw("COLZ");
      strcpy(outputname,"moller_ThetaE");
      c1->SaveAs(strcat(outputname,"_moller_thetaE_singles1.png"));

      moller_E->Draw("COLZ");
      strcpy(outputname,"moller_E");
      c1->SaveAs(strcat(outputname,"_moller_E_singles1.png"));

      moller_Theta->Draw("COLZ");
      strcpy(outputname,"moller_Theta");
      c1->SaveAs(strcat(outputname,"_moller_Theta_singles1.png"));

      moller_ESum->Draw("COLZ");
      strcpy(outputname,"moller_ESum");
      c1->SaveAs(strcat(outputname,"_moller_ESum_singles1.png"));

}
