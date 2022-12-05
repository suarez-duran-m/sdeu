# -*- coding: utf-8 -*-
import ROOT
import glob
import os

FULL  = 0  # read everything
MICRO = 1  # only shower level observables
MINI  = 2  # only traces are turned off


# Locate and load library
for path in \
    ([os.environ[x] + "/lib" for x in ("AUGEROFFLINEROOT", "ADSTROOT") if x in os.environ] +
        os.environ["LD_LIBRARY_PATH"].split(":")):
    path = path.strip()
    if not path:
        continue
    fileNameList = glob.glob(path + "/*RecEvent*")
    if fileNameList:
        ROOT.gSystem.Load(fileNameList[0])
        break
else:
    raise EnvironmentError(
        "Could not find RecEventKG library. "
        "You might have to set the environment variable \"AUGEROFFLINEROOT\"/\"ADSTROOT\" "
        "to your local offline/adst installation.")


# internal parser for input argument
def __ParseArgument(inputArgument):
    fileNameList = []
    if hasattr(inputArgument, "__len__"):
        fileNameList = inputArgument
    else:
        if inputArgument.endswith(".root"):
            fileNameList = [inputArgument]
        elif os.path.exists(inputArgument):
            for line in file(inputArgument):
                line = line.strip()
                if line[0] != "#":
                    fileNameList.append(line)
        else:
            raise IOError("Expecting {0} to be either a path to a ROOT file, \
                           a list of paths or path to a file containing paths"
                          .format(inputArgument))

    if not fileNameList:
        raise IOError("list of filenames is empty")

    fileNameList.sort()
    return fileNameList


def RecEventProvider(inputArgument, mode=FULL):
    """
    Provides an iterator over RecEvent objects saved in ADST files.

    Parameters
    ----------
    inputArgument: string or list of strings
      This can be: A path to a ROOT file, a path to a text file, which contains
      one path of a ROOT file per line, a list of paths to several ROOT files.
    mode: FULL, MICRO, MINI or list of strings
      FULL = read everything (default)
      MICRO = read only high-level observables (station info is turned off!)
      MINI = read everything except traces
      If this is a list of strings it is interpreted as a list of branches to
      be read out, while everything else is ignored (for the power user).

    Examples
    --------
    >>> for event in RecEventProvider("example.root"):
    ...   print event.GetSDEvent().GetEventId()

    Authors
    -------
    Hans Dembinski <hans.dembinski@kit.edu>
    """

    # BEWARE: We don't use SmartOpen, it has a memory leak.
    filenames = __ParseArgument(inputArgument)

    vfilenames = ROOT.std.vector("string")()
    for filename in filenames:
        vfilenames.push_back(filename)
    del filenames

    rf = ROOT.RecEventFile(vfilenames)
    if hasattr(mode, "__iter__"):
        rf.SetBranchStatus("event.*", 0)
        for x in mode:
            rf.SetBranchStatus(x, 1)
    elif mode == MICRO:
        rf.SetMicroADST()
    elif mode == MINI:
        rf.SetMiniADST()

    rev = ROOT.RecEvent()
    rf.SetBuffers(rev)

    detector_geometry = ROOT.DetectorGeometry()
    rf.ReadDetectorGeometry(detector_geometry)

    while rf.ReadNextEvent() == ROOT.RecEventFile.eSuccess:
        yield rev

    rf.Close()
    del rf
    del vfilenames
    del rev


def GetDetectorGeometry(inputArgument):
    """
    Returns the latest DetectorGeometry object from one or several ADST files.

    Parameters
    ----------
    inputArgument: string or list of strings
      Either path of a ROOT file, path of a text file, which contains one
      path of a ROOT file per line, or list of paths to ROOT files.

    Notes
    -----
    If you work on several ADST files, some might have an older
    DetectorGeometry stored than others, where not all SD stations/FD eyes
    where available yet. This function returns the latest DetectorGeometry
    object in this case, whereas "latest" means the one with most SD stations.

    Examples
    --------
    >>> detgeom = GetDetectorGeometry("example.root")

    Authors
    -------
    Hans Dembinski <hans.dembinski@kit.edu>
    """

    # BEWARE: We don't use SmartOpen, it has a memory leak.
    fileNameList = __ParseArgument(inputArgument)

    best = ROOT.DetectorGeometry()
    for fileName in fileNameList:
        rf = ROOT.RecEventFile(fileName)
        geometry = ROOT.DetectorGeometry()
        rf.ReadDetectorGeometry(geometry)
        if geometry.GetNStations() > best.GetNStations():
            best = geometry
    return best


def GetFileInfo(inputFileName):
    """
    Returns the FileInfo-object of an ADST file.
    This can be usefull to access the OffLine Configuration
    used to create an ADST file. (see example)

    Parameters
    ----------
    inputFileName: string
      Path of a ROOT file.

    Examples
    --------
    >>> fileinfo = GetFileInfo("example.root")
    >>> fileinfo.GetOfflineConfiguration()

    Authors
    -------
    Benjamin Fuchs <Benjamin.Fuchs@kit.edu>
    """

    fileinfo = ROOT.FileInfo()
    rf = ROOT.RecEventFile(inputFileName)
    rf.ReadFileInfo(fileinfo)
    return fileinfo
