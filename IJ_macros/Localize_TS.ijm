path 	= "/Volumes/CSieben_SSD/4_data_PALM/2020-07-01_ACE2_HEK/"; 
respath = "/Volumes/CSieben_SSD/4_data_PALM/2020-07-01_ACE2_HEK/locResults/";
name 	= "HEK_gain300_20ms_37deg_PCglas_Cell3";

run("Bio-Formats Importer", "open="+path+name+".nd2");

run("Z Project...", "projection=[Standard Deviation]");
run("16-bit");
saveAs("Tiff",respath+"STD_"+name+".tif");
close();

// selectWindow(path+name+".nd2");
selectWindow(name+".nd2");

run("Camera setup", "offset=100.0 quantumefficiency=0.9 isemgain=true photons2adu=6.0"+ 
	"gainem=300.0 pixelsize=160.0");

run("Run analysis", "filter=[Wavelet filter (B-Spline)] scale=2.0 order=3 detector "+
	"detector=[Local maximum] connectivity=8-neighbourhood threshold=std(Wave.F1) "+
	"estimator=[PSF: Integrated Gaussian] sigma=1.6 method=[Weighted Least squares] "+ 
	"fitradius=3 mfaenabled=false renderer=[Averaged shifted histograms] magnification=5.0 "+ 
	"colorizez=true shifts=2 repaint=50 threed=false");

run("Export results", "floatprecision=3 filepath="+respath+name+".csv " + 
	"fileformat=[CSV (comma separated)] sigma=true intensity=true chi2=false "+
	"offset=false saveprotocol=false x=true y=true bkgstd=false id=false "+
	"uncertainty_xy=true frame=true");

run("Close All");
	

