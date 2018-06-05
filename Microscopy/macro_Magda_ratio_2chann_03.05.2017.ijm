// Install the BIOP Library (from PTBIOP update site)
call("BIOP_LibInstaller.installLibrary", "BIOP"+File.separator+"BIOPLib.ijm");
	
run("Set Measurements...", "area limit display redirect=None decimal=3");

chDapi  = 1;
thDapi = "IsoData";

chLipid = 2;
thLipid = "Li";

//processImage();
//waitForUser;

//get directory 
imageFolder = getImageFolder();
saveFolder = getSaveFolder();
imageNumber = getNumberImages();

setBatchMode(true);
for (imageIndex = 0 ; imageIndex < imageNumber ; imageIndex++){
	openImage(imageIndex);					// open the image (an its assiciated roiset.zip)
	processImage();							// process the image
	saveRois("Save");						// save the ROIset with the current image
	saveCurrentImage();						// save the current image
	run("Close All");						// close all
}

if( isOpen("Results") ){
	selectWindow("Results");
	saveAs("Results", saveFolder+"results.txt");// save the results
}

setBatchMode(false);

showMessage("Jobs DONE!");





// required functions 

function toolName() {
	return "macro Magda";
}

function processImage(){
	getDimensions(width, height, channels, slices, frames);
	getVoxelSize(widthPixel, heightPixel, depthPixel, unitPixel);

	// do some processing

	title = getTitle();
	run("Properties...", "unit=pixel pixel_width=1 pixel_height=1 voxel_depth=1");
	run("Split Channels");
	
	selectImage("C"+chLipid+"-"+title);
	run("Gaussian Blur...", "sigma=2");
	setAutoThreshold(thLipid+" dark");
	//waitForUser;
	
	run("Measure");
	setResult("Channel", nResults-1, "Lipid");
	setResult("Label", nResults-1, title);
	
	setOption("BlackBackground", true);
	run("Convert to Mask");
	setMinAndMax(200, 255);
	
	selectImage("C"+chDapi+"-"+title);
	run("Gaussian Blur...", "sigma=3");
	setAutoThreshold(thDapi+" dark");
	//waitForUser;
	
	run("Measure");
	setResult("Channel", nResults-1, "DAPI");
	setResult("Label", nResults-1, title);
	
	setOption("BlackBackground", true);
	run("Convert to Mask");
	setMinAndMax(200, 300);
	
	run("Merge Channels...", "c1=C"+chLipid+"-"+title+" c2=C"+chDapi+"-"+title+" create");
	selectImage("Composite");
	//setSlice(1);
	Stack.setChannel(1);
	run("Yellow");

	//setSlice(2);
	Stack.setChannel(2);
	run("Blue");
	
	run("RGB Color");
	// rename to the original name again
	rename("Processed - "+title);
}