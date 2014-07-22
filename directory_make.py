'''Script to generate an index of Software Carpentry
lessons from the novice and intermediate directories 
of the bc repository'''

import glob
import sys

##Takes one argument: the path name for bc

script = sys.argv[0]
directory = sys.argv[1]

##directory names to pull from
##(could automate with * in the for loop below
##but this allows you to exempt non-relevant directories)

lesson_difficulty = ["novice", "intermediate"]
lesson_topic = ["shell","git","python","r","sql","regex","doit","make","extras"]

##output file in markdown
output_file = open("directory.md",'w')

#looping through directories
#loop writes header levels and then links
for difficulty in lesson_difficulty:
	output_file.write("#"+difficulty.capitalize()+"\n")
	for topic in lesson_topic:
		output_file.write("##"+topic+"\n")
		for file in glob.glob(directory+"/"+difficulty+"/"+topic+"/[0-9]*.md"):
			input_file = open(file, 'r')
			for line in input_file:
				if line[0:5] == "title":
					title = line[7:].strip()
					output_file.write("* ["+title+"]("+file+")\n")
				if line
			input_file.close()

output_file.close()

# for difficulty in lesson_difficulty:
	# output_file.write("#"+difficulty.capitalize()+"\n")
	# for topic in lesson_topic:
		# output_file.write("##"+topic+"\n")
		# for file in glob.glob(directory+"/"+difficulty+"/"+topic+"/[0-9]*.md"):
			# input_file = open(file, 'r')
			# for line in input_file:
				# if line[0:3] == "## ":
					# title = line[3:].strip()
					# output_file.write("* ["+title+"]("+file+")\n")
			# input_file.close()

# output_file.close()