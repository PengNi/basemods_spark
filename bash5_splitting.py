#!/usr/bin/python
import numpy


def split_holenumbers(holenumbers, folds=1):
    holenumbers_len = len(holenumbers)
    if folds > holenumbers_len:
        folds = holenumbers_len
    onefoldlen_base = holenumbers_len / folds
    fold_yu = holenumbers_len % folds

    hole_splitspots = [onefoldlen_base] * folds
    hole_splitspots[0] += fold_yu
    endOffset = numpy.cumsum(hole_splitspots)
    beginOffset = numpy.hstack(([0], endOffset[0:-1]))
    offsets = zip(beginOffset, endOffset)
    return offsets


def makeOffsetsDataStructure(baxh5obj):
    """
    from pbcore.io.BasH5IO.py
    :param baxh5obj:
    :return:
    """
    numEvent = baxh5obj["/PulseData/BaseCalls/ZMW/NumEvent"].value
    holeNumber = baxh5obj["/PulseData/BaseCalls/ZMW/HoleNumber"].value
    endOffset = numpy.cumsum(numEvent)
    beginOffset = numpy.hstack(([0], endOffset[0:-1]))
    offsets = zip(beginOffset, endOffset)
    return dict(zip(holeNumber, offsets))


def get_movieName(baxh5file):
    """
    from pbcore.io.BasH5IO.py.movieName
    :param baxh5file:
    :return:
    """
    movieNameAttr = baxh5file["/ScanData/RunInfo"].attrs["MovieName"]

    # In old bas.h5 files, attributes of ScanData/RunInfo are stored as
    # strings in arrays of length one.
    if (isinstance(movieNameAttr, (numpy.ndarray, list)) and
                len(movieNameAttr) == 1):
        movieNameString = movieNameAttr[0]
    else:
        movieNameString = movieNameAttr

    if not isinstance(movieNameString, basestring):
        raise TypeError("Unsupported movieName {m} of type {t}."
                        .format(m=movieNameString,
                                t=type(movieNameString)))
    return movieNameString
