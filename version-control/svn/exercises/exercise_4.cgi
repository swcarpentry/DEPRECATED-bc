#!/usr/bin/python
# Accepts a username, and then commits a modification to the "planets.txt" file
# in that users' folder to cause a conflict for exercise 4.
#
# Expects that the repository is structured so that there is a
# /participants folder containing a folder with each username.
#
import cgi
import cgitb
svn_username = "swc_robot"
svn_password = "yggdrasil"
svn_repo_url = 'http://svn.software-carpentry.org/spring2011'
planets_filename = 'planets.txt'

cgitb.enable()
print "Content-Type: text/html"
print
print "<html><title>Version Control Exercise #4</title><body>"

def update_repo(username, tempdir_path):
	import subprocess
	import sys
	import os.path
	print "Using username:", username, "<br/><br/>"

	# fetch list of participants and make sure username is one of them
	svn_path = svn_repo_url + "/participants/"
	svn_ls = subprocess.Popen(
		["svn", "ls", "--username", svn_username, "--password", svn_password, svn_path], stdout=subprocess.PIPE)
	userlist = svn_ls.communicate()[0].split("/\n")

	if username not in userlist:
		raise Exception("Unknown user:" + username)

	# check out a copy of the user's files
	svn_path 		= svn_path + username
	working_copy_path 	= os.path.join(tempdir_path,username)
	planets			= os.path.join(working_copy_path, planets_filename)

	cmd = "svn co " + svn_path + " " + working_copy_path
	cmd = subprocess.Popen(
		["svn", "co", "--username", svn_username, "--password", svn_password, 
		svn_path, working_copy_path], stdout=subprocess.PIPE).communicate()[0]
	

	# modify the planets file 
	try:
		planets_file    = open(planets,"r")
	except: 
		raise Exception("Unable to find or open '%s'! Are you sure it's in your user folder?" % planets_filename)

	contents 		= planets_file.readlines()
	planets_file.close()
	if "miles" not in contents[1].lower():
		raise Exception("Uh oh.  I was expecting to find the second line " + 
			"of %s say '* 10^6 miles', but it doesn't.  I won't continue." % planets_filename)
	contents[1] = "\t\t* 1000000 miles\r\n"
	contents = "".join(contents)
	planets_file = open(planets, "w")
	planets_file.write(contents)
	planets_file.close()

	# check it back in
	cmd = subprocess.Popen(
		["svn", "commit", "--username", svn_username, "--password", svn_password, 
		"-m", "Updating planets file via the web interface",
		working_copy_path], stdout=subprocess.PIPE).communicate()[0]

	print '<div style="color: green;">', "Success! The %s file should have been modified by the user: %s" % (planets_filename, svn_username), '</div>'

if __name__ == '__main__':
	form = cgi.FieldStorage()
	if "username" not in form:
		print 'Enter your username. (Only put in your own username... be nice.)<br/>' 
		print 'When you press submit the planets.txt file in your repository folder will be modified.'
		print '<br/>'
		print "<form>"
		print '   Username: <input type="text" name="username"/>'
		print '<input type="submit" value="Submit"/>'
		print '</form>'
	else:
		import re
		import shutil
		import tempfile
		try:
			username = re.sub("\W","", form.getfirst("username",""))
			tempdir_path = tempfile.mkdtemp()
			update_repo(username, tempdir_path)
		except Exception, e:
			print '<div style="color: red;">Error: ', e, "</div>"
		finally: 
			shutil.rmtree(tempdir_path)

print "</body></html>"
