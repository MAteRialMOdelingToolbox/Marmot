import os

# change to doc directory
os.chdir("doc")

# run doxgen
os.system("doxygen Doxyfile")

# run sphinx
os.system("sphinx-build -b html -Dbreathe_projects.Marmot=doc_out/xml . doc_out/sphinx")
