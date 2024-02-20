'''
for pythonanywhere
'''
import sys

# add project directory to path
projectHome = u'/home/bio_pam_blosum'
if projectHome not in sys.path:
    sys.path = [projectHome] + sys.path

#pass flask app as "application" for wsgi to work
# for dash, app.server
from app import app
application = app.server