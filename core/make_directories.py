from os import makedirs, path
import errno

def makeDirectories(dirList):
    
    for directory in dirList:

        try:
            makedirs(directory, 0o777)
        except OSError as exc:
            if exc.errno == errno.EEXIST and path.isdir(directory):
                pass
            else: raise
