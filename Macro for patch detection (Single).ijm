//This macro measures patch size and density

//all variables in the code to be adjusted

//Set Image Scale in µm !!!!!!!!!!!!!!!! change appropriately !!!!!!!!!!!!!!!!
//scale = 465.67; // this is for EXP2 and EXP4 Rag10,12,13,19,20,21,22,23,24,28,29,30,31,32,33
scale = 927.02; // this is for EXP5 Rag14,15,16,17,18,25,26,27,34,35,36

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

//Patch size selection
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

//Marker line thickness
yY = 10; //yellow line; //full size = 10; half size = 5
yR = 20; //white line; //full size = 20; half size = 10
yW = 20; //white line; //full size = 20; half size = 10

//label dimesions for label in density image
zoomMag = 12; //full size = 12; half size = 7

//--------------------------------------------------------------------------------------------------------------------------------

//clear previously used RAM from previous runs
run("Collect Garbage");

//activates Bio-Formats Importer options
run("Bio-Formats Macro Extensions");

///gets image name, extention and directory
imagePath = File.openDialog(" "); //Full file path with file name and extention
imageName = File.getName(imagePath); //extracts file name with extention
dotIndex = lastIndexOf(imageName, "."); //index for the . between file name and extention
fileLength = lengthOf(imageName); //gets the total length of the file name.
fileExt = substring(imageName, (dotIndex+1), fileLength); //gets the extention
imageName = substring(imageName, 0, dotIndex); //gets the file name without extention
pathLength = lengthOf(imagePath); //gets the total length of the whole file path
dir1 = substring(imagePath, 0, (pathLength-fileLength)); //gets the directory of the image path

//opens file without opening Bio-Fomats importer window
run("Bio-Formats Windowless Importer", "open=["+dir1+"/"+imageName+"."+fileExt+"]");

//saving patch and density images and the results before closing
dir2 = dir1+"/../../OutputWhole"+File.separator; //this creates the new path; "/../"selects parent folder
File.makeDirectory(dir2); //this creates the new path
path = dir2+imageName+"."+fileExt;//this creates the new path
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
run("Set Measurements...", "area mean standard min centroid center perimeter fit skewness kurtosis display redirect=[Area Measurements] decimal=2");
run("Analyze Particles...", "size="+patchSize+"-Infinity show=[Masks] display record"); //measure selected areas by color threshold
rename("NND");
image4 = getTitle();
selectWindow("Threshold Color"); //select and close Color Threshold Window
run("Close");

if (isOpen("Results")) {//the following two if-conditions are supposed to help to get around the problem when no patches are detected and thus, no selection can be created

	//Nearest Neighbour Distance (NND) and average distance to nearest 6 neichbors measurements
	//Haeri & Haeri (2015), 3, e1
	selectWindow(image4);
	run("ND ");
	selectWindow("Distance Between Neighboring Particles");
	saveAs("Results", path+"_NND.csv");
	run("Close");
	run("Collect Garbage");
	
	//Clustering patches by max distance (Graph)
	//Ben Tupper
	selectWindow(image4);
	run("Create Selection");
	run("Colors...", "foreground=white background=black selection=red"); //adjust marker line colour
	run("Line Width...", "line="+yR); //modifies boarder line width
	run("Graph ", "mask="+image4+" distance=Edges neighbors="+maxDis+" lines label matrix");
	selectWindow("Distance Matrix");
	saveAs("Results", path+"_Distance_Matrix.csv");
	run("Close");
	selectWindow("Log");
	saveAs("Text", path+"_Clusters.csv");
	run("Close");
	selectWindow("Adjacencies");
	run("Flatten");
	close("Adjacencies");
	run("Collect Garbage");

	//get whole skin boarders
	selectWindow(image3);
	run("RGB Color");
	rename("Cluster Boarder");
	run("Colors...", "foreground=white background=black selection=black"); //adjust marker line colour
	run("8-bit");
	//setAutoThreshold("Default dark");
	run("Threshold...");
	setThreshold(1, 255);
	run("Create Selection");
	selectWindow("Adjacencies-1");
	run("Restore Selection");
	run("Properties... ", "  width="+yW);
	run("Flatten");
	save(path+"_distance.tif");
	close("Cluster Boarder");
	close("Adjacencies-1");
	close("Adjacencies-2");
	run("Collect Garbage");

	//get complete patch coverage
	selectWindow(image4);
	run("Create Selection");
	run("Colors...", "foreground=white background=black selection=yellow"); //adjust marker line colour
	run("Set Measurements...", "area mean standard min center skewness kurtosis display redirect=None decimal=2");
	selectWindow(image3);
	run("Restore Selection");
	run("Properties... ", "  width="+yY); //modifies boarder line width
	run("ROI Manager...");
	roiManager("Add");
	roiManager("Select", 0);
	roiManager("Rename", "Sum Patch Area");
	roiManager("Measure");
				
	//create image with marked boarders for measured patches
	selectWindow(image3);
	run("Flatten");
	run("Select None");
	run("Restore Selection");
	run("Flatten");
	image3 = getTitle();
	selectWindow("ROI Manager");
	run("Close");
	run("Collect Garbage");
}

//get whole skin area and skin boarders
selectWindow(image3);
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

//close windows not needed
selectWindow("Threshold");
run("Close");
close("Area Measurements-1");
close("Mask of Area Measurements");
close("Mask of Area Measurements-1");
close("Whole Skin");
close("Skin Boarder");
close("Red Green Merge (RGB)");
close("NND");
run("Collect Garbage");
	
//generate image density images
selectWindow(image2);
rename("Density Image");
run("royal");
run("Calibration Bar...", "location=[Upper Right] fill=White label=Black number=5 decimal=0 font=12 zoom="+zoomMag+" bold overlay");
run("Flatten");
close("Density Image");
run("Collect Garbage");

//save
save(path+"_density.tif");
close();
save(path+"_patches.tif");
close();
saveAs("Results", path+"_Results.csv");
selectWindow("Results");
run("Close");

//clear previously used RAM from previous runs
run("Collect Garbage");

		
