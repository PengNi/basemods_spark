#!/usr/bin/python

import xml.etree.ElementTree as ET


def getMoviesChemistry(filepaths):
    """

    :param filepaths: list of filepath
    :return:
    """
    moviestriple = {}
    for filepath in filepaths:
        moviestriple.update(MovieChemistry(filepath).getMovieTriple())
    return moviestriple


class ChemistryLookupError(Exception):
    pass


class MovieChemistry:
    def __init__(self, filepath):
        self._movietriple = {}
        if str(filepath).endswith(".xml"):
            moviename = str(filepath).split("/")[-1].split(".")[0]
            triple = self.tripleFromMetadataXML(filepath)

            tripledict = dict()
            # tripledict['Name'] = moviename
            tripledict['BindingKit'] = triple[0]
            tripledict['SequencingKit'] = triple[1]
            tripledict['SoftwareVersion'] = triple[2]
            self._movietriple[moviename] = tripledict
        else:
            pass

    def getMovieTriple(self):
        return self._movietriple

    def tripleFromMetadataXML(self, metadataXmlPath):
        """
        from pbcore.io.BasH5IO
        Scrape the triple from the metadata.xml, or exception if the file
        or the relevant contents are not found
        """
        nsd = {None: "http://pacificbiosciences.com/PAP/Metadata.xsd",
               "pb": "http://pacificbiosciences.com/PAP/Metadata.xsd"}
        try:
            tree = ET.parse(metadataXmlPath)
            root = tree.getroot()
            bindingKit = root.find("pb:BindingKit/pb:PartNumber", namespaces=nsd).text
            sequencingKit = root.find("pb:SequencingKit/pb:PartNumber", namespaces=nsd).text
            # The instrument version is truncated to the first 3 dot components
            instrumentControlVersion = root.find("pb:InstCtrlVer", namespaces=nsd).text
            verComponents = instrumentControlVersion.split(".")[0:2]
            instrumentControlVersion = ".".join(verComponents)
            return (bindingKit, sequencingKit, instrumentControlVersion)
        except Exception as e:
            raise ChemistryLookupError, \
                ("Could not find, or extract chemistry information from, %s" % (metadataXmlPath,))
