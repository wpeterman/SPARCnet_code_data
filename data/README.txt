File List
DateMat_osu.csv
egg_model_normal_posteriors.csv
fitted_model_df.rds
Lotter_Data.csv
RSmat.csv
SalDatGrowth.csv
sCJS_Data.rds
SexMat.csv

*****************************

File: 
	DateMat_osu.csv

Contents: 
	Dates that salamander surveys were conducted

Fields:
	Unit: 1–4 identfier
	Plot: Plot identifier
	Survey_x: Date that survey (1–6) was conducted
	Site_num: Numeric identifer of ridge (1) and slope (2)
	Plot_Num: Numeric indicator of plot

*****************************

File: 
	egg_model_normal_posteriors.csv

Contents: 
	20,0000 Posterior estimates of linear regression model parameters from model predicting salamander clutch size in relation to snout-vent-length

Fields:
	b_Intercept: Intercept term
	b_svl: beta slope term
	sigma: sigma term

*****************************

File: 
	Lotter_Data.csv

Contents: 
	Data extracted from Lotter 1978
	Lotter, F. 1978. Reproductive ecology of the salamander Plethodon cinereus (Amphibia, Urodela, Plethodontidae) in Connecticut. Journal of Herpetology 12:231–236.


Fields:
	svl: Snout-vent-length, measured in mm
	eggs: Number of eggs counted in females

*****************************

File: 
	RSmat.csv

Contents: 
	Matrix indicating whether observed salamanders are from Ridge or Slope plots


Fields:
	UniqueID: Unique code, including plot identifier
	Sx_is_R: Six columns indicating whether or not an individual was observed on a Ridge plot (1) or a Slope plot (0) during each of the six surveys
	
*****************************

File: 
	SalDatGrowth.csv

Contents: 
	Table containing capture-recapture information for each salamander


Fields:
	[blank]: Nameless column with numeric row identfier
	season: Character column indicating early vs late sample survey within each year
	UniqueID: Unique identification code assigned to each salamander
	SVL: snount-vent-length measured in mm
	Sex: Male (M), Female (F), or Unknown (U)
	RS: Character indicating whether a salamander was observed on a Ridge plot (R) or Slope plot (S)
	PosID: Numeric 1 = slope; 2 = ridge
	Plot: Alphanumeric plot identifier
	date: date that a salamander was observed
	SurveyNum: Numeric survey number (1–6) that a salamander was observed

*****************************

File: 
	SexMat.csv

Contents: 
	Matrix indicating the sex of each salamander


Fields:	
	UniqueID: Unique identification code assigned to each salamander
	Sx_is_M: Indicator variable as to whether a salamander was identified as being male (1) in surveys 1–6
	Sx_is_F: Indicator variable as to whether a salamander was identified as being female (1) in surveys 1–6
	
*****************************

File: 
	fitted_model_df.rds

Contents: 
	R data frame containing 175,000 posterior draws from fitted growth model


Fields:	
	size_male: Estimated asymptotic snout-vent-length for males
	size_female: Estimated asymptotic snout-vent-length for males
	K_slope: Growth parameter estimate for slope
	K_ridge: Growth parameter estimate for ridge
	K_male: Growth parameter estimate for males
	K_ridge_m: Growth parameter estimate for males on ridges
	K_ridge_f: Growth parameter estimate for females on ridges
	K_slope_m: Growth parameter estimate for males on slopes
	K_slope_f: Growth parameter estimate for females on slopes
	sigma.SVL: Uncertainty in estimated snout-vent-length
	deviance: Deviance estimate from fitted Bayesian model
	
*****************************

File: 
	sCJS_Data.rds

Contents: 
	List of 20 R data objects, including all elements needed to run sCJS model