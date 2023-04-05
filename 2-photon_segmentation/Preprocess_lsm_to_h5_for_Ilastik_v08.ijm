//open the series file manually 
dir = getDirectory("Choose a Directory ");
setBatchMode(true);
//images are saved here. You will need to create this manually with "C1", "C2" and C3" subfolders. 
dir_save = dir+"hdf5"+File.separator;
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
  for (i=0; i<list.length; i++) {
      if (endsWith(list[i], "/"))
          processFiles(""+dir+list[i]);
      else {
         showProgress(n++, count);
         path = dir+list[i];
         processFile(path);
      }
  }
}

function processFile(path) {
   if (endsWith(path, ".lsm")) {
       open(path);
       
       	//get current image title
	filename = getTitle();

	//get parent directory
	parents = split(File.getParent(path),File.separator);
	parent_1 = parents[parents.length -1]; //e.g. 2h
	parent_2 = parents[parents.length -2]; //e.g. mouse name
	parent_3 = parents[parents.length -3]; //e.g. fMCAO
		
	//split channels 
	run("Split Channels");
	wait(500);
	
	//select C3 + save as single h5 
	selectWindow("C3-" + filename);
	run("Stack to Images");
	run("Images to Stack", " title=R4");
	selectWindow("Stack");
	run("Scriptable save HDF5 (new or replace)...", "save=["+dir_save+"C3"+File.separator+parent_3+"_"+parent_2+"_"+parent_1+"_"+filename+"_C3.h5] dsetnametemplate=/channel{c} formattime=%d formatchannel=%d compressionlevel=0");
	selectWindow("Stack");
	close();
	wait(500);
	
	//select C2 + save as single h5 
	selectWindow("C2-" + filename);
	run("Stack to Images");
	run("Images to Stack", " title=R3");
	run("Scriptable save HDF5 (new or replace)...", "save=["+dir_save+"C2"+File.separator+parent_3+"_"+parent_2+"_"+parent_1+"_"+filename+"_C2.h5] dsetnametemplate=/channel{c} formattime=%d formatchannel=%d compressionlevel=0");
	selectWindow("Stack");
	close();
	wait(500);
	
	//select C1 + save as single h5 
	selectWindow("C1-" + filename);
	run("Stack to Images");
	run("Images to Stack", " title=R2");
	run("Scriptable save HDF5 (new or replace)...", "save=["+dir_save+"C1"+File.separator+parent_3+"_"+parent_2+"_"+parent_1+"_"+filename+"_C1.h5] dsetnametemplate=/channel{c} formattime=%d formatchannel=%d compressionlevel=0");
	selectWindow("Stack");
	close();
	wait(500);
		
  }
}
