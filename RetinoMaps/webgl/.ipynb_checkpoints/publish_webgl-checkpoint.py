"""
-----------------------------------------------------------------------------------------
publish_webgl.py
-----------------------------------------------------------------------------------------
Goal of the script:
Publish webgl on webapp
-----------------------------------------------------------------------------------------
Input(s):
sys.argv[1]: main project directory
sys.argv[2]: project name (correspond to directory)
-----------------------------------------------------------------------------------------
Output(s):
sent data
-----------------------------------------------------------------------------------------
To run:
On invibe.nohost.me
1. cd to function
>> cd cd ~/disks/meso_H/projects/[PROJECT]/analysis_code/postproc/prf/webgl/
2. run python command
>> python publish_webgl.py [main directory] [project name]
-----------------------------------------------------------------------------------------
Exemple:
cd ~/disks/meso_H/projects/RetinoMaps/analysis_code/postproc/prf/webgl/
python publish_webgl.py ~/disks/meso_S/data RetinoMaps
-----------------------------------------------------------------------------------------
Written by Martin Szinte (mail@martinszinte.net)
Edited by Uriel Lascombes (uriel.lascombes@laposte.net)
-----------------------------------------------------------------------------------------
"""

# Stop warnings
import warnings
warnings.filterwarnings("ignore")

# Debug import 
import ipdb
deb = ipdb.set_trace

# General imports
import os
import sys
import json

# Inputs
main_dir = sys.argv[1]
project_dir = sys.argv[2]

# Define analysis parameters
with open('../../../settings.json') as f:
    json_s = f.read()
    analysis_info = json.loads(json_s)
webapp_login = analysis_info['webapp_login']
webapp_dir = analysis_info['webapp_dir']

# Define directories
index_dir = "{}".format(os.getcwd())
webgl_dir = "{}/{}/derivatives/webgl".format(main_dir, project_dir)
    
# Copy index.html file in webgl folder
print("Saving {}/index.html".format(webgl_dir))
os.system("rsync -avuz {}/index.html {}/index.html".format(index_dir, webgl_dir))

# send webgl folder
print("Sending {}".format(webgl_dir))
print("Please use server admin password")
os.system("rsync -avuz --progress {}/ {}:{}/www/".format(webgl_dir, webapp_login, webapp_dir))
os.system("ssh {} chmod -Rfv 777 {}/".format(webapp_login, webapp_dir))