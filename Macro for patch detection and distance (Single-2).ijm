//This macro measures patch size and density

//all variables in the code to be adjusted

//Adjustable variables
imageName = "Rag14"; //image name as it appears in the folder; name of the mouse here
fileExt = ".psd"; //defines the file type loaded, e.g. .tif, .psd (Photoshop), etc.

//Set Image Scale in µm !!!!!!!!!!!!!!!! change appropriately !!!!!!!!!!!!!!!!
//scale = 465.67; // this is for EXP2 and EXP4 Rag10,12,13,19,20,21,22,23,24,28,29,30,31,32,33
//scale = 927.02; // this is for EXP5 Rag14,15,16,17,18,25,26,27,34,35,36

//Color Threshold for patch selection
//Hue
Hmin = 0;
Hmax = 255;
//Saturation
Smin = 0;
Smax = 255;
//Brightness
Bmin = 75;
Bmax = 255;

//Patch size selection
patchSize = 1960; //according to Guertin & Sabatini (2015) [doi: 10.1038/npg.els.0003359], the average macrophage diameter is 20-30 µm;
//if 25 µm is the mean diameter for macrophages then one macrophage has an average 2D surface of ~490 µm^2, which means 4 macrophages have a combined 2D surface of ~1960 µm^2

//Marker line thickness
yY = 10; //yellow line; //full size = 10; half size = 5
yW = 20; //white line; //full size = 20; half size = 10

//measuring line thickness
lineWidth = 20;

//label dimesions for label in density image
zoomMag = 12; //full size = 12; half size = 7

//For Distance Measurement
//define earhole size to exclude abitrary skin tears / gaps
earHole = 20000000 //20,000,000 µm^2
//Selection of Ear by row; row1 = 0, row2 = 1
mouseEar = 0;

//Converison of XM and XY to x and y
conv = 10.82;

//define last line for the loop
loopEnd = 3;

//--------------------------------------------------------------------------------------------------------------------------------

//clear previously used RAM from previous runs
run("Collect Garbage");

//activates Bio-Formats Importer options
run("Bio-Formats Macro Extensions");

//gets directory from where to load image
dir1 = getDirectory(" "); //select your file location

//opens file without opening Bio-Fomats importer window
run("Bio-Formats Windowless Importer", "open=["+dir1+"/"+imageName+fileExt+"]");

//saving patch and density images and the results before closing
dir2 = dir1+"/../OutputDist"+File.separator; //this creates the new path; "/../"selects parent folder
File.makeDirectory(dir2); //this creates the new path
path = dir2+imageName+fileExt;//this creates the new path
dotIndex = lastIndexOf(path, "."); 
path = substring(path, 0, dotIndex); // remove ext.

//set scale
run("Set Scale...", "distance="+scale+" known=10000 pixel=1 unit=µm global");

//split colors to get red only
run("Split Channels");
Blue = getTitle();
close(Blue);
Green=getTitle();
run("Put Behind [tab]");
Red=getTitle();
run("Duplicate...", " ");
rename("Area Measurements");
image2 = getTitle();
run("Merge Channels...", "c1=["+Red+"] c2=["+Green+"] create");
rename("Red Green Merge");
image3 = getTitle();
run("Duplicate...", "Ear Holes");
run("RGB Color");		

//find the ear holes in the skin
run("8-bit");
//setAutoThreshold("Default dark");
run("Threshold...");
setThreshold(0, 1);
run("Set Measurements...", "center display redirect=None decimal=2");
run("Analyze Particles...", "size="+earHole+"-Infinity show=Masks display exclude add");
selectWindow("Threshold");
run("Close");
close("Ear Holes");

//analysis for patches
selectWindow(image2);
run("Duplicate...", " ");
run("RGB Color");
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

//Measuring patch areas mean gray values and centre of mass
run("Set Measurements...", "area mean standard min center skewness kurtosis display redirect=[Area Measurements] decimal=2");
run("Analyze Particles...", "size="+patchSize+"-Infinity show=[Masks] display"); //measure selected areas by color threshold
run("Create Selection");
run("Properties... ", "  width="+yY); //modifies boarder line width
selectWindow("Threshold Color"); //select and close Color Threshold Window
run("Close");
				
//get complete patch coverage
run("Set Measurements...", "area mean standard min center skewness kurtosis display redirect=None decimal=2");
selectWindow(image3);
run("Restore Selection");
run("ROI Manager...");
roiManager("Add");
roiManager("Select", 0);
roiManager("Rename", "Sum Patch Area");
roiManager("Measure");
selectWindow("ROI Manager");
run("Close");
			
//create image with marked boarders for measured patches
selectWindow(image3);
run("Flatten");
run("Restore Selection");
run("Flatten");

//get whole skin area and skin boarders
rename("Whole Skin");
image5=getTitle(); //loads image with patch boarders into variable
run("Duplicate...", " ");
rename("Skin Boarder");
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

//saving patch image
save(path+"_patches.tif");
image6=getTitle();

//close windows not needed
selectWindow("Threshold");
run("Close");
		
//generate image density images
selectWindow(image2);
rename("Density Image");
run("royal");
run("Calibration Bar...", "location=[Upper Right] fill=White label=Black number=5 decimal=0 font=12 zoom="+zoomMag+" bold overlay");
run("Flatten");

//saving density image
save(path+"_density.tif");

//select image for distance measurement and close all else to recover RAM
selectWindow(image6);
close("\\Others");
run("Collect Garbage");

// open a sample image
maxCount = nResults - loopEnd;
selectWindow(image6);
run("Line Width...", "line="+lineWidth); //adjusts line thickness
run("Colors...", "foreground=white background=black selection=green"); //adjust marker line colour

//run loop without opening new images
setBatchMode(true);

for (i = 1; i < maxCount; i++) {
	//get coordinates from results table
	originX = getResult("XM", mouseEar);
	originY = getResult("YM", mouseEar);
	destinX = getResult("XM", 1 + i);
	destinY = getResult("YM", 1 + i);
	// make an example selection
	;xpoints = newArray(originX/conv, destinX/conv);
	ypoints = newArray(originY/conv, destinY/conv);
	makeSelection("point", xpoints, ypoints);
	//measures the length of the drawn line	
	Roi.getCoordinates(xpoints, ypoints); // store point coordinates
	makeLine(xpoints[0],ypoints[0],xpoints[1],ypoints[1]); // make line
	run("Measure"); // measure the line
	run("Flatten");
	close("\\Others");
	run("Collect Garbage");
}
setBatchMode(false);

//saving distance image and results
save(path+"_distance.tif");
close();
saveAs("Results", path+"_Results.csv");
selectWindow("Results");
run("Close");

//clear previously used RAM from previous runs
run("Collect Garbage");