
const double m_proton=0.938272046; //GeV/c2
const int binsize=5; // cm

double dedx_min=30.;

//x-y cut ---------------------------//
double mean_x=-26.58; //prod4
double mean_y=423.5; //prod4
double dev_x=1.5*3.753; //prod4
double dev_y=1.5*4.354; //prod4

double mean_x_mc=-29.25; //prod4 mc
double mean_y_mc=422.1; //prod4 mc
double dev_x_mc=1.5*4.456; //prod4 mc
double dev_y_mc=1.5*3.933; //prod4 mc


//double mean_x=-26.41;
//double mean_y=423.6;
//double dev_x=1.5*3.687;
//double dev_y=1.5*4.349;
//double dev_x=100000.;
//double dev_y=100000.;
//x-y cut ---------------------------//

//z0 offset (before SCE corr. ------------------------------------------//
double z0_before_sce=32.4672;
double sigma0_before_sce=0.872971;


	//mu of beam ---------------------------------------------------------------//
	double mu_data=4.45451e+02;
	double err_mu_data=7.21974e-01;
	double sigma_data=5.09380e+01;
	double err_sigma_data=5.22759e-01;

	double mu_stop_data=4.05686e+02; //range-based calc
	double err_mu_stop_data=1.23011e+00;
	double sigma_stop_data=5.25443e+01;
	double err_sigma_stop_data=1.20796e+00;

	double mu_calo_stop_data=3.84030e+02; //calo-based calc
	double err_mu_calo_stop_data=1.15843e+00;
	double sigma_calo_stop_data=5.65971e+01;
	double err_sigma_calo_stop_data=9.02569e-01;

	//double de_upstream_data=mu_data-mu_stop_data;
	//double err_de_upstream_data=sqrt(pow(err_mu_data,2)+pow(err_mu_stop_data,2));
	//double de_upstream_data=mu_data-mu_calo_stop_data;
	double de_upstream_data=mu_data-mu_stop_data;
	//double err_de_upstream_data=sqrt(pow(err_mu_data,2)+pow(err_mu_calo_stop_data,2));
	//--------------------------------------------------------------------------//

		//track selection cuts -------------------------------------------------------------------------------------------------------------------------------------------//
		//stopping protons
		//double mean_norm_trklen_csda=8.37593e-01;
		//double sigma_norm_trklen_csda=7.09110e-02;
		//double mean_norm_trklen_csda=8.49529e-01; //new
		//double sigma_norm_trklen_csda=7.32112e-02; //new
		//double mean_norm_trklen_csda=8.82596e-01; //new2
		//double sigma_norm_trklen_csda=7.06573e-02; //new2
		double mean_norm_trklen_csda=8.84707e-01; //prod4
		double sigma_norm_trklen_csda=7.02962e-02; //prod4

		//double min_norm_trklen_csda=0.8; //trk length cut
		double min_norm_trklen_csda=mean_norm_trklen_csda-2.*sigma_norm_trklen_csda;
		double max_norm_trklen_csda=mean_norm_trklen_csda+3.*sigma_norm_trklen_csda;
		//double min_norm_trklen_csda=mean_norm_trklen_csda-1.*sigma_norm_trklen_csda;
		//double max_norm_trklen_csda=mean_norm_trklen_csda+1.*sigma_norm_trklen_csda;
		//double min_cosine=9.78573e-01-4.*1.10895e-02; //cosine cut(0.934), new2
		double min_cosine=9.78550e-01-4.*1.10031e-02; //cosine cut(0.934), prod4

		//position cuts (diff)
		double min1_dz=3.67683e+00-3.*1.01896e+00; //new2
		double min2_dz=3.67683e+00+3.*1.01896e+00; //new2
		double min1_dy=1.85320e+00-3.*1.87173e+00;
		double min2_dy=1.85320e+00+3.*1.87173e+00;
		double min1_dx=3.17667e+00-3.*1.22454e+00; //new2
		double min2_dx=3.17667e+00+3.*1.22454e+00; //new2

		//position cuts
		//double min1_z=3.67892e+00-3.*1.02152e+00; //new2
		//double min2_z=3.67892e+00+3.*1.02152e+00; //new2
		//double min1_y=4.24055e+02-3.*4.54622e+00; //new2
		//double min2_y=4.24055e+02+3.*4.54622e+00; //new2
		//double min1_x=-2.82470e+01-3.*3.83924e+00; //new2
		//double min2_x=-2.82470e+01+3.*3.83924e+00; //new2

		double min1_x=-2.82351e+01-3.*3.98756e+00; //prod4
		double min2_x=-2.82351e+01+3.*3.98756e+00; //prod4
		double min1_y=4.24140e+02-3.*4.68599e+00; //prod4
		double min2_y=4.24140e+02+3.*4.68599e+00; //prod4
		double min1_z=3.79565e+00-3.*1.01964e+00; //prod4
		double min2_z=3.79565e+00+3.*1.01964e+00; //prod4



