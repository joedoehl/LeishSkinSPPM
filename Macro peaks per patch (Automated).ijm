//#####################################################################################################################################

// Title:        Skin patch analysis

// Script Name:  Macro Complete Skin Analysis (Single).ijm
// Copyright:    Johannes S. P. Doehl (NIAID) & Paul M. Kaye (University of York)

// Purpose:      Analysis of visceralizing Leishmania-infected whole-body skin images for the extraction of the following data:
//				 1 - extraction of x & y coordinates of the skin ouline and ear wholes, if present
//				 2 - detection of skin patches, measurement of patch size and density and nearest enighbour analysis
//				 3 - distance measurements of patches from the right ear

// Author:       Johannes S. P. Doehl
// Date:         March 2020

// Publication: 

// Acknowledgements: Adapted from the protocol from Andrew McNaughton (2013) from the Unviersity of Otago

//#####################################################################################################################################

//SCRIPT

//clear previously used RAM from previous runs
run("Collect Garbage");

//activates Bio-Formats Importer options
run("Bio-Formats Macro Extensions");

//get image directory and extract images
dir1 = getDirectory(" "); //change this according to file location
list = getFileList(dir1);

for (ii=0; ii<list.length; ii++) { //begin of outer loop
	path = dir1+list[ii];

//#####################################################################################################################################

//ADJUSTABLE VARIABLES

//Set Image Scale in µm !!!!!!!!!!!!!!!! change appropriately !!!!!!!!!!!!!!!!
//scale = 465.67; // adjust according to your image resolution by measuring scale bar
scale = 927.02; // adjust according to your image resolution by measuring scale bar

//Marker line thickness
yG = 20; //red line; //full size = 20; half size = 10

//Color Threshold for patch selection
//Hue
Hmin = 0;
Hmax = 255;
//Saturation
Smin = 0;
Smax = 255;
//Brightness
Bmin = 70;
Bmin2 = 170;
Bmax = 255;

//Mimimum size that describes a parasite patch
patchSize = 4148.65; //µm^2
//according to Guertin & Sabatini (2015) [doi: 10.1038/npg.els.0003359], the average macrophage diameter is 20-30 µm;
//if 25 µm is the mean diameter for macrophages, then one macrophage has an average 2D surface of ~490 µm^2
//a pixel is in half size about 21.47 µm x 21.47 µm = 461 µm^2
//Thus, 4 macrophages will require at least 3 x 3 = 9 pixels to be covered; Thus, 64.41 µm x 64.42 µm = 4148.65  µm^2
//Since, the "scale" is adjusted for full sized images, the pixel number will dobule but the dimensions in µm will stay the same

patchMulti = 2; // mutiplier for the patchSize to determine peak size

//#####################################################################################################################################

//SCRIPT

	//opens file without opening Bio-Fomats importer window
	run("Bio-Formats Windowless Importer", "open=[path]");
	
	//set scale
	run("Set Scale...", "distance="+scale+" known=10000 pixel=1 unit=µm global");
	
	//split colors to get red only
	run("Split Channels");
	Blue = getTitle();
	close(Blue);
	Green=getTitle();
	run("Put Behind [tab]");
	Red=getTitle();
	run("Duplicate...", "title=[Area Measurements]");
	image2 = getTitle();
	run("Merge Channels...", "c1=["+Red+"] c2=["+Green+"] create");
	rename("Red Green Merge");
	run("Flatten");
	image3 = getTitle();
	close("Red Green Merge");
	
	//-------------------------------------------------------------------------------------------------------------------------------------
	
	//saving patch and density images and the results before closing
	dir2 = dir1+"/../../OutputPeaks"+File.separator; //this creates the new path; "/../"selects parent folder
	File.makeDirectory(dir2); //this creates the new path
	path = dir2+list[ii];//this creates the new path
	dotIndex = lastIndexOf(path, "."); 
	path = substring(path, 0, dotIndex); // remove ext.
	
	selectWindow(image2);
	setBatchMode(true);
	run("Duplicate...", " ");
	run("RGB Color");	
	//analysis for patches
	run("Colors...", "foreground=white background=black selection=green"); //adjust marker line colour
	run("Line Width... ", "  width="+yG); //modifies boarder line width
	run("Color Threshold...");
	// Color Thresholder 2.0.0-rc-69/1.52n
	min=newArray(3);
	max=newArray(3);
	filter=newArray(3);
	a=getTitle();
	run("HSB Stack");
	run("Convert Stack to Images");
	selectWindow("Hue");
	rename("0");
	selectWindow("Saturation");
	rename("1");
	selectWindow("Brightness");
	rename("2");
	min[0]=Hmin;
	max[0]=Hmax;
	filter[0]="pass";
	min[1]=Smin;
	max[1]=Smax;
	filter[1]="pass";
	min[2]=Bmin;
	max[2]=Bmax;
	filter[2]="pass";
	for (i=0;i<3;i++){
	  selectWindow(""+i);
	  setThreshold(min[i], max[i]);
	  run("Convert to Mask");
	  if (filter[i]=="stop")  run("Invert");
	}
	imageCalculator("AND create", "0","1");
	imageCalculator("AND create", "Result of 0","2");
	for (i=0;i<3;i++) {
	  selectWindow(""+i);
	  close();
	}
	selectWindow("Result of 0");
	close();
	selectWindow("Result of Result of 0");
	rename(a);
	// Colour Thresholding-------------end
	setBatchMode(false);
	//Measuring patch areas mean gray values and centre of mass
	run("Set Measurements...", "area mean standard min centroid center perimeter fit skewness kurtosis stack display redirect=[Area Measurements] decimal=2");
	run("Analyze Particles...", "size="+patchSize+"-Infinity show=[Masks] display exclude record add"); //measure selected areas by color threshold
	rename("NND");
	maxCounts = nResults;
	print("Patch Count");
	print(maxCounts);
	selectWindow("Threshold Color"); //select and close Color Threshold Window
	run("Close");
	close("Area Measurements");
	selectWindow("Area Measurements-1");
	run("Flatten");
	save(path+"_patches.tif");
	close("Area Measurements-1");
	close("Area Measurements-2");
	close("NND");
	selectWindow("Results");
	run("Close");
	run("Collect Garbage");
	
	//-------------------------------------------
	
	peakSize = patchSize * patchMulti;
	selectWindow(image3);
	run("Color Threshold...");
	setBatchMode(true);
	// Color Thresholder 2.0.0-rc-69/1.52p
	// Autogenerated macro, single images only!
	min=newArray(3);
	max=newArray(3);
	filter=newArray(3);
	a=getTitle();
	run("HSB Stack");
	run("Convert Stack to Images");
	selectWindow("Hue");
	rename("0");
	selectWindow("Saturation");
	rename("1");
	selectWindow("Brightness");
	rename("2");
	min[0]=Hmin;
	max[0]=Hmax;
	filter[0]="pass";
	min[1]=Smin;
	max[1]=Smax;
	filter[1]="pass";
	min[2]=Bmin2;
	max[2]=Bmax;
	filter[2]="pass";
	for (i=0;i<3;i++){
	  selectWindow(""+i);
	  setThreshold(min[i], max[i]);
	  run("Convert to Mask");
	  if (filter[i]=="stop")  run("Invert");
	}
	imageCalculator("AND create", "0","1");
	imageCalculator("AND create", "Result of 0","2");
	for (i=0;i<3;i++){
	  selectWindow(""+i);
	  close();
	}
	selectWindow("Result of 0");
	close();
	selectWindow("Result of Result of 0");
	rename(a);
	// Colour Thresholding-------------
	setBatchMode(false);
	run("Set Measurements...", "area center display redirect=None decimal=2");
	run("Analyze Particles...", "size="+peakSize+"-Infinity show=[Masks] display exclude record"); //measure selected areas by color threshold
	run("Invert");
	rename("Peak Counts");
	image5 = getTitle();
	
	maxPeaks = nResults;
	print("Peak Count");
	print(maxPeaks);
	selectWindow("Threshold Color"); //select and close Color Threshold Window
	run("Close");
	selectWindow("Results");
	run("Close");
	close("Red Green Merge (RGB)");
	selectWindow("Log");
	saveAs("Text", path+"_Counts.csv");
	selectWindow("Log");
	run("Close");
	
	//-------------------------------------------
	
	for (i=0;i<maxCounts;i++) {
		selectWindow(image5);
		roiManager("Select", i);
		run("Find Maxima...", "prominence=20 exclude output=Count");
		close("Mask of Peak Counts");
		run("Collect Garbage");
	}
	
	selectWindow("Results");
	saveAs("Results", path+"_Peak_Counts.csv");
	selectWindow("Results");
	run("Close");
	selectWindow(image5);
	save(path+"_peaks.tif");
	close(image5);
	selectWindow("ROI Manager");
	run("Close");
	run("Collect Garbage");
}