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

//--------------------------------------------------------------------------------------------------------------------------------
//all variables in the code to be adjusted
//!!!!for some reason the macro works more relaibly with the variables within the loop to be redefined every loop around

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
	//according to Guertin & Sabatini (2006) [doi: 10.1038/npg.els.0003359], the average macrophage diameter is 20-30 µm;
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
	yY = 20; //yellow line; //full size = 10; half size = 5
	yR = 20; //white line; //full size = 20; half size = 10
	yW = 20; //white line; //full size = 20; half size = 10

	//label dimesions for label in density image
	zoomMag = 12; //full size = 12; half size = 7

//--------------------------------------------------------------------------------------------------------------------------------

    //opens file without opening Bio-Fomats importer window
    run("Bio-Formats Windowless Importer", "open=[path]");

	//saving patch and density images and the results before closing
	dir2 = dir1+"/../../OutputWhole"+File.separator; //this creates the new path; "/../"selects parent folder
	File.makeDirectory(dir2); //this creates the new path
	path = dir2+list[ii];//this creates the new path
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
	
	//Measuring patch areas mean gray values and centre of mass
	run("Set Measurements...", "area mean standard min centroid center perimeter fit skewness kurtosis display redirect=[Area Measurements] decimal=2");
	run("Analyze Particles...", "size="+patchSize+"-Infinity show=[Masks] display record"); //measure selected areas by color threshold
	rename("NND");
	image4 = getTitle();

	//generate image density images
	selectWindow(image2);
	rename("Density Image");
	run("royal");
	run("Calibration Bar...", "location=[Upper Right] fill=White label=Black number=5 decimal=0 font=12 zoom="+zoomMag+" bold overlay");
	run("Flatten");
	image6 = getTitle();
	run("Collect Garbage");
		
	//create image with marked boarders for measured patches
	selectWindow(image4);
	run("Create Selection");
	run("Colors...", "foreground=white background=black selection=magenta"); //adjust marker line colour
	selectWindow(image6);
	run("Restore Selection");
	run("Properties... ", "  width="+yY);
	run("Flatten");
	image6 = getTitle();
	run("Collect Garbage");
			
	//get whole skin area and skin boarders
	selectWindow(image3);
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
	selectWindow(image6);
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
		
	
	//save
	save(path+"_density.tif");
	close();
	selectWindow("Results");
	run("Close");
	run("Close All");
	
	//clear previously used RAM from previous runs
	run("Collect Garbage");
} //closes the first loop
	
