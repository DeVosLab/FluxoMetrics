/*	
 	**********************************************
 	Fluxometrics.ijm 
 	**********************************************	
 	Authors: 				Winnok H. De Vos & Hugo Steenberghen
	Date Created: 			March 21, 2020
	Date Last Modified:		May 5th, 2022

 	This macro set measures intensity fluctuations (e.g., from the self-quenching dye calcein-AM) in individual cells (ROIs segmented from average ROI)
	The image takes single channel time series stacks of different types (nd2, tiff...) as input. 
	It returns a table containing the cell area and mean instensity as rows, for each cell ROI in a separate column.
	
	1. Install macro set via Plugins > Macros > Install...
	2. Click "S" to set up the segmentation.
		- To use StarDist, first Add "CSBDeep" and "StarDist" to the FIJI update sites.
		- To use presegmented imagej roi-sets, choose "Import_ROIs" for the autothreshold method (e.g. to first segment cells with Cellpose)
	3. Click "T" to test the segmentation on the input stack. 
		- This will average project the stack and segment the cells, including watershed separation by default.
		- The macro will return the segmented image and cell ROIs.
		- Clicking the red "T" will toggle the overlay on and off.
	4. Click "1" to test the whole process (average projection, segmentation, feature extraction).
	5. Click "#" to run the whole process on a folder in batch mode. 
		- All images with the chosen extension will be processed, regardless of the number of timepoints.
		- If checked, a folder with average projected images will be generated, which can be used to check the segmentation on afterwards.
		
 	**********************************************
*/


/*	
 	**********************************************

		Variable initiation

	**********************************************
*/

//	Strings
var dir 				= 	"";											//	directory of single mosaic
var out_dir				=	"";											//	dir for analyze output
var av_im_dir			=   "";											//  dir for AVG projecions
var list 				= 	newArray(0);								//	list of images withn dir
var log_path 			= 	"";											//	path for the log file
var bits 				= 	"16-bit";									//	image bitdepth
var file_type			=	".nd2";										//	image type
var micron				= 	getInfo("micrometer.abbreviation");			//	micron symbol
var unit				=	"pixel";									//	pixel calibration unit
var threshold			= 	"Import_ROIs";									//	autothreshold method for segmentation
var enhance				=	"gauss";									//	enhancement method
		
//	Number variables
var width				= 	1608;										//	image width
var height 				= 	1608;										//	image height
var slices 				= 	1;											//	number of slices
var frames 				= 	1;											//	number of frames	
var pixel_size			= 	0.364;										//	pixel size 20x objective (µm)
var blur_radius			=	1;											// 	radius of Gaussion 
var fixed				= 	700;											// 	value for fixed threshold
var prominence			= 	2000;										// 	prominence of objects for find maxima
var min_cell_size 		= 	50;											//	Min cell size (in pixels)
var max_cell_size 		= 	5000;										//	Max cell size (in pixels)
var roi_nr				= 	0;											// 	number of detected rois
var av_max				= 	0;											//  maximum intensity of average projection image		


//	Arrays
var thresholds 			= 	getList("threshold.methods");				//	Gets all autothreshold methods
var thresholds			= 	Array.concat(thresholds,"Fixed","Import_ROIs");			//	Include a fixed threshold option and manual import (as label images)
var file_types 			=	newArray(".nd2",".tif",".ids",".jpg");		//	file type for finding prefix 

//  Boolean				
var use_stardist		= false;										// checkbox for the use of stardist for segmentation
var av_save 			= true;											// checkbox for saving average images


/*	
 	**********************************************

		Macros

	**********************************************
*/

macro "Autorun"
{
	erase(1);
	setOptions();
}

macro "Settings Action Tool - C438 T5f16S"
{
	setup();
}

macro "Test segmentation settings Action Tool - C555L32c2L33c3L7484L7585L7686L7787L7888L7989L7a8aL7b8bL7c8c"
{
	erase(0);
	setBatchMode(false);
	id = getImageID;
	title = getTitle();
	check_dims(id);
	if(threshold == "Import_ROIs")
	{
		importROIs(id);
		roiManager("Show All");
	}
	else
	{
		if(slices>1)av_project(id);
		segment(id);
	}
	roi_nr = roiManager("count");
	print("number of rois retained: "+roi_nr);
	setBatchMode("exit and display");
}

macro "Analyse Single Image Action Tool - C438 T5f161"
{
	erase(0);
	setBatchMode(false);
	id = getImageID;
	title = getTitle();
	// Measure whole image
	run("Select All");
	roiManager("Add");
	analyze(id);
	selectWindow("Results"); 
	IJ.renameResults("Image Results");
	// Generate rois
	check_dims(id);
	if(threshold == "Import_ROIs")importROIs(id);
	else
	{
		if(slices>1)av_project(id);
		segment(id);
	}
	roi_nr = roiManager("count");
	print("number of rois: "+roi_nr);
	if(roi_nr>0){
		analyze(id);
		roiManager("Deselect");
		//roiManager("Multi Measure");
		selectWindow("Results"); 
		IJ.renameResults("ROI Results");
	}
	else {print("No rois detected");}
	setBatchMode("exit and display");
}

macro "Batch analysis Action Tool - C438 T5f16#"
{	
	erase(2);
	setDirectory();
	files = listFiles(dir);
	setBatchMode(false);
	for(i=0;i<files.length;i++)
	{
		run("Bio-Formats Importer", "open=["+files[i]+"] color_mode=Default view=Hyperstack stack_order=XYCZT");
		id = getImageID;
		title = getTitle;
		prefix = substring(title,0,indexOf(title,file_type));
		print("File",i+1,"/",files.length,":",prefix);
		// Measure whole image
		run("Select All");
		roiManager("Add");
		analyze(id);
		selectWindow("Results");
		saveAs(".txt",out_dir + prefix + "_im_results.txt");
		erase(0);
		// Segment or import rois
		if(threshold == "Import_ROIs")importROIs(id);
		check_dims(id);
		if(threshold == "Import_ROIs")importROIs(id);
		else
		{
			if(slices>1)av_project(id);
			segment(id);
		}
		// Analyse rois (if present)
		if(roi_nr>0){
			analyze(id);
			selectImage(id);close();
			close("*");
			selectWindow("ROI Manager");
			roiManager("save",out_dir + prefix + "_rois.zip");
			selectWindow("Results");
			saveAs(".txt",out_dir + prefix + "_roi_results.txt");
		}
		roiManager("reset");
		erase(0);
	}
	selectWindow("Log");
	saveAs(".txt",log_path);
	erase(2);
	print("Batch Analysis completed");
	setBatchMode("exit and display");		
}

macro "Toggle Overlay Action Tool - Cf88 T5f16T"
{
	 toggleOverlay();
}

macro "[t] Toggle Overlay"
{
	 toggleOverlay();
}

/*	
 	**********************************************

		Functions

	**********************************************
*/

function erase(type)
{
	if(type==1)print("\\Clear");
	if(type==2)
	{
		run("Close All");
		windows = getList("window.titles"); 
     	for (i=0; i<windows.length; i++)
     	{
     		selectWindow(windows[i]);
     		if(windows[i]!="Log")run("Close");
     	} 
	}
	run("Clear Results");
	roiManager("reset");
	//run("Collect Garbage");
}

function setup()
{
	// set segmentation parameters
	erase(1);
	call("ij.ImagePlus.setDefault16bitRange", 16);												// set default bit range to 16 (affects getMinMax(), does not affect 8-bit images)
	setOptions();
	Dialog.create("FluxoMetrics");
	Dialog.addChoice("File Type",file_types,file_type);
	Dialog.addNumber("Pixel Size (ignored if calibrated)",pixel_size,3,7,micron);
	Dialog.addCheckbox("Use Stardist segmentation", use_stardist);
	Dialog.addMessage("----Segmentation parameters if stardist is not used----");
	Dialog.addNumber("Smoothing Radius",blur_radius,1,7,"");
	Dialog.addChoice("Autothreshold Algorithm",thresholds,threshold);
	Dialog.addNumber("Fixed threshold value (-1 = not used)", fixed,0,7,"");
	Dialog.addNumber("Prominence",prominence,0,7,"");
	Dialog.addMessage("---------------------------Other parameters---------------------------");
	Dialog.addNumber("Minimum Cell Area ",min_cell_size,0,7,"pixels");
	Dialog.addNumber("Maximum Cell Area ",max_cell_size,0,7,"pixels");
	Dialog.addCheckbox("Save average projections?", av_save);
	Dialog.show;
	moment = getMoment();
	print(moment+"\n");
	print("*******************************************");
	print("               Settings");
	print("*******************************************");
	file_type 				= Dialog.getChoice();	print("image file type: ",file_type);
	pixel_size				= Dialog.getNumber();	print("image pixel size: ",pixel_size);
	use_stardist			= Dialog.getCheckbox();	print("use stardist: ", use_stardist);
	blur_radius				= Dialog.getNumber();	if(use_stardist== false){print("smoothing radius: ",blur_radius);}
	threshold 				= Dialog.getChoice();	if(use_stardist== false){print("threshold: ",threshold);}
	fixed					= Dialog.getNumber();	if(use_stardist== false){print("fixed threshold: ", fixed);}
	prominence				= Dialog.getNumber();	if(use_stardist== false){print("prominence: ",prominence);}
	min_cell_size 			= Dialog.getNumber();	print("min cell size: ",min_cell_size);
	max_cell_size 			= Dialog.getNumber();	print("max cell size: ",max_cell_size);
	av_save					= Dialog.getCheckbox(); print("save average projections: ", av_save);
	print("*******************************************");	
}

function setOptions()
{
	run("Options...", "iterations = 1 count = 1");
	run("Colors...", "foreground=white background=black selection=yellow");
	run("Appearance...", "  antialiased menu = 0");
	run("Overlay Options...", "stroke = red width = 1 fill = none");
	setOption("BlackBackground", false);
	run("Set Measurements...", "area mean redirect=None decimal=4");
	run("Misc...", "divide=Infinity reverse");
	run("Clear Results");
	setBackgroundColor(0, 0, 0);
	setForegroundColor(255,255,255);
}


function getMoment()
{
     MonthNames = newArray("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
     DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
     getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
     time_string ="Date: "+DayNames[dayOfWeek]+" ";
     if (dayOfMonth<10) {time_string = time_string+"0";}
     time_string = time_string+dayOfMonth+"-"+MonthNames[month]+"-"+year+"\nTime: ";
     if (hour<10) {time_string = time_string+"0";}
     time_string = time_string+hour+":";
     if (minute<10) {time_string = time_string+"0";}
     time_string = time_string+minute+":";
     if (second<10) {time_string = time_string+"0";}
     time_string = time_string+second;
     return time_string;
}


function setDirectory()
{
	dir 					= getDirectory("Choose a Source Directory");
	list 					= getFileList(dir);
	out_dir 				= dir+"Output"+File.separator;
	if(!File.exists(out_dir))File.makeDirectory(out_dir);
	av_im_dir = out_dir + "AVG_images"+File.separator;
	if(!File.exists(av_im_dir))File.makeDirectory(av_im_dir);
	log_path				= out_dir+"log.txt";
}

function listFiles(dir)
{
	files = newArray(0);
	for(i = 0; i < list.length; i++)
	{
		path = dir + list[i];
		if(endsWith(path,file_type))
		{
			files = Array.concat(files,path); 
		}
	}
	return files;
}

function importROIs(id){
	// Import label images from "Labels" subfolder and return ROIs
	selectImage(id);
	Stack.getDimensions(width, height, channels, slices, frames); 
	print(channels,"Channels -",slices,"Slices -",frames,"Frames");
	if(slices == 1 && frames > 1) 
	{
		run("Properties...", "channels=1 slices="+frames+" frames=1");
		print("new slices:"+slices);
	}
	title_temp = getTitle();
	
	// Lookup label image
	Imdir = getDirectory("image");
	open(Imdir+File.separator+"Labels"+File.separator+title_temp);
	label_id = getImageID();
	run("Label image to ROIs", "rm=[RoiManager[visible=true]]");
	selectImage(label_id);close;
	print("Imported labels");
}

function check_dims(id){
	selectImage(id);
	Stack.getDimensions(width, height, channels, slices, frames); 
	print(channels,"Channels -",slices,"Slices -",frames,"Frames");
	if(slices == 1 && frames > 1) 
	{
		run("Properties...", "channels=1 slices="+frames+" frames=1");
		print("new slices:"+slices);
	}
}

function av_project(id){
	print("projecting average image...");
	run("Z Project...", " projection=[Average Intensity]");
	av_id = getImageID;
	av_title = getTitle();
	getMinAndMax(av_min, av_max);
	run("Select None");
	if(av_save == true){
		saveAs("Tiff", av_im_dir+ av_title);
	}
	selectImage(av_id);
}


function segment(id)
{
	// segment cells and return ROIs thereof
	selectImage(id);
	Stack.getDimensions(width, height, channels, slices, frames); 
	if(slices == 1 && frames > 1) 
	{
		run("Properties...", "channels=1 slices="+frames+" frames=1");
	}
	
	
		if(use_stardist==false){
			run("Gaussian Blur...", "sigma="+blur_radius);
			print("fixed = ",fixed);
			if(fixed == -1){
				setAutoThreshold(threshold+" dark");
				//setOption("BlackBackground", false);
				//run("Convert to Mask");
			}
			else{
				print("Fixed threshold");
				setThreshold(fixed, av_max);
			}
			run("Find Maxima...", "prominence="+prominence+" above output=[Segmented Particles]");
			max_id = getImageID; 
			selectImage(max_id);
			setThreshold(1,255);
			run("Analyze Particles...", "size="+min_cell_size+"-"+max_cell_size+" show=Nothing exclude clear add");
			selectImage(max_id); close;
			roi_nr = roiManager("count");
		}
		else 
		{
			print("using stardist...");
			av_title = getTitle();
			print(av_title);
			run("Command From Macro", "command=[de.csbdresden.stardist.StarDist2D], args=['input':'"+av_title+"', 'modelChoice':'Versatile (fluorescent nuclei)', 'normalizeInput':'true', 'percentileBottom':'1.0', 'percentileTop':'99.0', 'probThresh':'0.30000000000000004', 'nmsThresh':'0.1', 'outputType':'ROI Manager', 'nTiles':'1', 'excludeBoundary':'2', 'roiPosition':'Automatic', 'verbose':'false', 'showCsbdeepProgress':'false', 'showProbAndDist':'false'], process=[false]");
			print("2");
			roiManager("Deselect");
			roiManager("Measure");
			roi_nr = roiManager("count");
			//print(roi_nr);
			for(i=roi_nr-1;i>=0;i--)
			{
				area = getResult("Area",i);
				//print(area);
				if(area<min_cell_size || area>max_cell_size)
				{
					roiManager("select",i);
					roiManager("delete");
				}
			}
		}
	
	print("rois detected: "+roi_nr);
	print(slices); //TEMP
	if(slices>1){selectImage(av_id); close;}
	selectImage(id);
	//selectImage(max_id); close;	
}

function analyze(id)
{
	selectImage(id);
	//roiManager("deselect");
	roiManager("multi measure");
	run("Select None");	
}

function toggleOverlay()
{
	run("Select None"); roiManager("deselect");
	roiManager("Show All without labels");
	if(Overlay.size == 0)run("From ROI Manager");
	else run("Remove Overlay");
}