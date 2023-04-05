//open the series file manually 
dir = getDirectory("Choose a Directory ");
setBatchMode(true);
count = 0;
countFiles(dir);
n = 0;
processFiles(dir);
//print(count+" files processed");

function countFiles(dir) {
  list = getFileList(dir);
  for (i=0; i<list.length; i++) {
      if (endsWith(list[i], "/"))
          countFiles(""+dir+list[i]);
      else
          count++;
  }
}

function processFiles(dir) {
  list = getFileList(dir);
  for (i=1; i<list.length; i++) {
      //print("trying to tackle "+list[i]);
      if (File.isDirectory(dir+list[i])){
      	//if it is a directory, go inside and check for tif files 
      	//Note this is not recursive, and therefore will not work with folders more than 1 layer deep 
      	subdir_list = getFileList(dir+list[i]);
      	if (endsWith(subdir_list[0],".tif")){
      		//if there are tifs inside, 
      		print("found tifs in "+list[i]+" , processing as stack");
	    	process_segmented_stack(dir+list[i]);
      	}
      }
  }
}


function process_segmented_stack(source_folder) {
	// close all images
	close("*");
	// empty the ROI manager
	roiManager("reset");
	// empty the results table
	run("Clear Results");
	
	//run("Image Sequence...", "dir=L:/Moritz/2022-03-16_2_photon_for_Ilgin/preliminary_results/Sham_010621_C2/ROI1_series.lsm_C2/ filter=.tif sort");
	source_folder = replace(source_folder,"\\","/");
	target_folder_list = getFileList(source_folder);
	print(target_folder_list[1]);
	//load image sequence
	run("Image Sequence...", "dir=[" + source_folder + "] filter=.tif sort");

	//get the name 
	filename = getTitle();
	filename = replace(filename,File.separator,"");

	//debug
	print("directory is "+source_folder);

	//threshold the probability map (so you get a binary mask)
	setAutoThreshold("Default dark");
	//run("Threshold...");
	setAutoThreshold("Default dark stack");
	setThreshold(118, 255);
	setOption("BlackBackground", false);
	run("Convert to Mask", "method=Default background=Dark");
	run("Watershed", "stack");
	

	//prepare the measurements
	run("Set Measurements...", "area area_fraction redirect=None decimal=3");
	run("Analyze Particles...", "size=100-Infinity clear summarize add stack");
	//saveAs("Results", "L:/Moritz/2022-03-16_2_photon_for_Ilgin/preliminary_results/Sham_010621_C2/ROI1_series.lsm_C2/Summary of ROI1_series.lsm_C2-1_summary_pixel_statisticcs.csv");
	saveAs("Results", source_folder+"Summary of "+filename+"_summary_pixel_statistics.csv");
	
	//split stack
	run("Stack to Images");

	//get the name of the split images
	single_image_name = getTitle();
	single_image_name = substring(single_image_name, 0, lengthOf(single_image_name)-2);

	
	//debug
	print("processing now "+single_image_name);

	//img 0 -> 1
	//selectWindow(single_image_name+"_0");
	//selectWindow(single_image_name+"_1");
	imageCalculator("Subtract create", single_image_name+"_0",single_image_name+"_1");
	selectWindow("Result of "+single_image_name+"_0");
	run("Open");
	run("Close-");
	run("Analyze Particles...", "size=100-Infinity display summarize");
	
	
	//img 1 -> 2
	//selectWindow(single_image_name+"_1");
	//selectWindow(single_image_name+"_2");
	imageCalculator("Subtract create", single_image_name+"_1",single_image_name+"_2");
	selectWindow("Result of "+single_image_name+"_1");
	run("Open");
	run("Close-");
	run("Analyze Particles...", "size=100-Infinity display summarize");
	
	
	//img 2 -> 3
	//selectWindow(single_image_name+"_2");
	//selectWindow(single_image_name+"_3");
	imageCalculator("Subtract create", single_image_name+"_2",single_image_name+"_3");
	selectWindow("Result of "+single_image_name+"_2");
	run("Open");
	run("Close-");
	run("Analyze Particles...", "size=100-Infinity display summarize");
	
	
	//img 3 -> 4
	//selectWindow(single_image_name+"_3");
	//selectWindow(single_image_name+"_4");
	imageCalculator("Subtract create", single_image_name+"_3",single_image_name+"_4");
	selectWindow("Result of "+single_image_name+"_3");
	run("Open");
	run("Close-");
	run("Analyze Particles...", "size=100-Infinity display summarize");
	//saveAs("Results", "L:/Moritz/2022-03-16_2_photon_for_Ilgin/preliminary_results/Sham_010621_C2/ROI1_series.lsm_C2/Summary of ROI1_series.lsm_C2-1_all.csv");
	saveAs("Results", source_folder+"Summary of "+ single_image_name + "_all.csv");
	
	run("Close All");
	close("Summary of "+ single_image_name + "_all.csv");
	close("Summary of "+filename+"_summary_pixel_statistics.csv");
	close("Results");
}


