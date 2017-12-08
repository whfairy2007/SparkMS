import os

from fnmatch import fnmatch

def file_name_find(root,pattern):
    FileName = []
    short_name = []
    for path, subdirs, files in os.walk(root):
        for name in files:
            if fnmatch(name, pattern):
                FileName.append(os.path.join(path, name))
                short_name.append(name)
    return(FileName,short_name)