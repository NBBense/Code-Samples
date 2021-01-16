#####################################################################################
# Title: preProcImages.py
# Author: Nicholas Bense
# Date: 4/26/19
# Description: Python script to pre-process/split SEM .tif images for CNN classification
#####################################################################################

from PIL import Image
import os, sys, argparse, numpy


def pre_proc_images(tif_directory):
	'''Function to pre-process and split SEM .tif images within provided directory to prepare 
	for CNN classification.'''
	
	# Loop through files in provided directory to extract all .tif file names and paths
	for file in os.listdir(tif_directory):
		if file.endswith('.tif'): 
			full_path = os.path.join(tif_directory, file)
			print(file)
			print(full_path)
			img = Image.open(full_path).convert("L")
			
			# Crop images to 1000x1000 to remove SEM instrument panel
			new_img = img.crop((0,0,1000,1000))
			img_arr = numpy.array(new_img)
			
			# Split images into 25 200x200 tiles
			tiles = [img_arr[x:x+200,y:y+200] 
				for x in range(0,img_arr.shape[0],200) 
				for y in range(0,img_arr.shape[1],200)]
			
			
			# Loop through tiles, convert to image and save output files
			os.mkdir(''.join([tif_directory,"/Preprocessed_Images/") 
			i = 0
			for tile in tiles:
				split_image = Image.fromarray(tile)
				split_image.save(''.join([tif_directory,
						"/Preprocessed_Images/",
						os.path.splitext(file)[0],
						"_tile_",
						str(i),
						".tif"]))
				i = i + 1
 
def main(argv):

	# Parse arguments for .tif directory
	parser = argparse.ArgumentParser()
	parser.add_argument("tif_directory", help="path to directory containing .tif files")
	args = parser.parse_args()
	
	# Call pre-processing function
	pre_proc_images(args.tif_directory)

if __name__ == "__main__":
   main(sys.argv[1:])

