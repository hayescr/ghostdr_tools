import os

def tryMakePath(path):
    try:
        os.mkdir(path)
    except FileExistsError:
        print('Directory {} already exists'.format(path))
    else:
        print('Directory {} created'.format(path))

def makeDirectoryStructure(DATA_parent, obj_name, readmodes):
    """Create's Fletcher's favourite GHOSTDR directory structure for easy trouble-shooting. Does not add your data"""

    dirpath = os.path.join('./', DATA_parent,)
    tryMakePath(dirpath)
    for i in obj_name:
        objpath = os.path.join(dirpath,  i)
        tryMakePath(objpath)
        for j in readmodes:
            readmodepath = os.path.join(objpath,  j)
            tryMakePath(readmodepath)
            rawpath = os.path.join(readmodepath,  'raw')
            tryMakePath(rawpath)
            tryMakePath(os.path.join(readmodepath,  'intermediate'))
            tryMakePath(os.path.join(rawpath,  'packed'))
            tryMakePath(os.path.join(rawpath,  'obj'))    