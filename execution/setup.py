import pathlib, sys


def setDirectoryTree():
    # List of directories to add to the system path
    mainPath = pathlib.Path(__file__).parent.absolute() / ".."
    dirList = ['', 'cases', 'execution', 'mesh', 'models', 'post', 'solvers', 'util', 'vectors', "fields"]

    # Add the directories
    for d in dirList:
        if d not in sys.path:  # Avoid re adding
            sys.path.append(mainPath + d)

