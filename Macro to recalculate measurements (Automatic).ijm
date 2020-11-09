//#####################################################################################################################################

// Title:        Skin patch analysis

// Script Name:  Macro Complete Skin Analysis (Automated).ijm
// Copyright:    Johannes S. P. Doehl (NIAID) & Paul M. Kaye (University of York)

// Purpose:      Analysis of visceralizing Leishmania-infected whole-body skin images for the extraction of the following data:
//				 1 - extraction of x & y coordinates of the skin ouline and ear wholes, if present
//				 2 - detection of skin patches, measurement of patch size and density and nearest enighbour analysis
//				 3 - distance measurements of patches from the right ear

// Author:       Johannes S. P. Doehl
// Date:         March 2020

// Publication:

//#####################################################################################################################################

//SCRIPT

//This macro measures patch size and density

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

//Converison of XM and XY to x and y
conv = 21.47444; //half size images
//conv = 10.78725; //full size images

//Marker line thickness
lw1 = 1; // for measurement
yY = 10; //yellow line; //full size = 10; half size = 5
yR = 20; //red line; //full size = 20; half size = 10
yW = 20; //white line; //full size = 20; half size = 10
yM = 10; //magenta line; //full size = 10; half size = 5
lineWidth = 40; //// line thickness of distance measurement line

//label dimesions for label in density image
zoomMag = 12; //full size = 12; half size = 7

//For Distance Measurement
//define earhole size to exclude abitrary skin tears / gaps
earHole = 30000000; //30,000,000 µm^2 = 30 mm^2 = 0.3 cm^2

//Color Threshold for patch selection
//Hue
Hmin = 0;
Hmax = 255;
//Saturation
Smin = 0;
Smax = 255;
//Brightness
Bmin = 70;
Bmax = 255;

//Mimimum size that describes a parasite patch
patchSize = 4148.65; //µm^2
//according to Guertin & Sabatini (2015) [doi: 10.1038/npg.els.0003359], the average macrophage diameter is 20-30 µm;
//if 25 µm is the mean diameter for macrophages, then one macrophage has an average 2D surface of ~490 µm^2
//a pixel is in half size about 21.47 µm x 21.47 µm = 461 µm^2
//Thus, 4 macrophages will require at least 3 x 3 = 9 pixels to be covered; Thus, 64.41 µm x 64.42 µm = 4148.65  µm^2
//Since, the "scale" is adjusted for full sized images, the pixel number will dobule but the dimensions in µm will stay the same

//Max distance allowed between patches for cluster analysis
maxDis = 3900; // µm
//Lutz & Neiva (1912) decibed Lutzomyia longipalpis female adult sand flies as about 2 mm full body length, a wingspan around 2.3 mm and a body with of about 6.5 mm
//Thus, if one assumes that one body length into any direction from where a sand fly lands on the skin can be considered a foraging area, then patches in a cluster should be max about 4 mm apart (edge to edge)
//Since patches will have to reach into this radius and the labrium is about 35 µm in diameter according to Brinson (1993)than from any side the patch needs to reach into the foragin area by at least 50 µm
//Thus, 4000 µm - 2x50 µm = 3900 µm diameter or max patch separation (edge to edge)

//define last line for the loop
loopEnd = 5;

//#####################################################################################################################################

//SCRIPT (continued)

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
dir2 = dir1+"/../../OutputWhole"+File.separator; //this creates the new path; "/../"selects parent folder
File.makeDirectory(dir2); //this creates the new path
path = dir2+list[ii];//this creates the new path
dotIndex = lastIndexOf(path, "."); 
path = substring(path, 0, dotIndex); // remove ext.

selectWindow(image2);
setBatchMode(true);
run("Duplicate...", " ");
run("RGB Color");	

//analysis for patches
run("Colors...", "foreground=white background=black selection=yellow"); //adjust marker line colour
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
for (i=0;i<3;i++){
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
run("Set Measurements...", "area mean standard min centroid center perimeter fit feret's skewness kurtosis stack display redirect=[Area Measurements] decimal=2");
run("Analyze Particles...", "size="+patchSize+"-Infinity show=[Masks] display exclude record"); //measure selected areas by color threshold
rename("NND");
image4 = getTitle();
maxCount = nResults;
print("Patch Counts");
print(maxCount);
selectWindow("Log");
saveAs("Text", path+"_patch_counts.csv");
run("Close");
selectWindow("Threshold Color"); //select and close Color Threshold Window
run("Close");

//-------------------------------------------

//generate image density images
selectWindow(image2);
rename("Density Image");
run("royal");
run("Calibration Bar...", "location=[Upper Right] fill=White label=Black number=5 decimal=0 font=12 zoom="+zoomMag+" bold overlay");
run("Flatten");
image5 = getTitle();
save(path+"_density.tif");
close("Density Image");
close("Area Measurements-1");
run("Collect Garbage");

//-------------------------------------------

//get complete patch coverage
selectWindow(image4);
run("Create Selection");
run("Colors...", "foreground=white background=black selection=yellow"); //adjust marker line colour
run("Set Measurements...", "area mean standard min centroid center perimeter fit skewness kurtosis stack display redirect=None decimal=2");
selectWindow(image3);
run("Duplicate...", "title=[Patch Boarders]");
image6=getTitle(); //loads image with patch boarders into variable
run("Restore Selection");
run("Properties... ", "  width="+yY); //modifies boarder line width
run("ROI Manager...");
roiManager("Add");
roiManager("Select", 0);
roiManager("Rename", "Sum Patch Area");
roiManager("Measure");
			
//create image with marked boarders for measured patches
selectWindow(image6);
run("Flatten");
run("Select None");
run("Restore Selection");
run("Flatten");
image6 = getTitle();
selectWindow(image5);
run("Colors...", "foreground=white background=black selection=magenta"); //adjust marker line colour
run("Restore Selection");
run("Properties... ", "  width="+yM); //modifies boarder line width
run("Flatten");
image5 = getTitle();
selectWindow("ROI Manager");
run("Close");
run("Collect Garbage");

//-------------------------------------------

//get whole skin area and skin boarders
selectWindow(image3);
run("Duplicate...", "title=[Whole Skin]");
image7=getTitle(); //loads image with patch boarders into variable
run("Duplicate...", "title=[Skin Boarder]");
run("Colors...", "foreground=white background=black selection=white"); //adjust marker line colour
run("RGB Color");
run("8-bit");
//setAutoThreshold("Default dark");
run("Threshold...");
setThreshold(1, 255);
run("Create Selection");
run("ROI Manager...");
roiManager("Add");
roiManager("Select", 0);
roiManager("Rename", "Skin Area");
roiManager("Measure");
run("Close");
selectWindow(image5);
run("Restore Selection");
run("Properties... ", "  width="+yW);
run("Flatten");
image5=getTitle(); //loads image with patch boarders into variable
selectWindow(image6);
run("Restore Selection");
run("Properties... ", "  width="+yW);
run("Flatten");
image6=getTitle(); //loads image with patch boarders into variable

//close windows not needed
selectWindow("Threshold");
run("Close");
close("Whole Skin");
close("Skin Boarder");
close("Red Green Merge (RGB)-1");
close("Red Green Merge (RGB)-2");
close("NND");
close("Density Image-1");
close("Density Image-2");
close("Patch Boarders");
close("Patch Boarders-1");
close("Patch Boarders-2");
run("Collect Garbage");

//-------------------------------------------

//save
selectWindow(image5);
save(path+"_densityboarder.tif");
close();
selectWindow(image6);
save(path+"_patches.tif");
selectWindow("Results");
saveAs("Results", path+"_Results.csv");
selectWindow("Results");
run("Close");
close("Red Green Merge (RGB)");
close("Patch Boarders-3");

//-------------------------------------------------------------------------------------------------------------------------------------

//clear previously used RAM from previous runs
run("Collect Garbage");

} //closes the first loop
	
