#!/usr/bin/python

COLUMNS = 60


def getRefInfoFromFastaFiles(filepaths):
    """

    :param filepaths: list of filepath
    :return:
    """
    contiginfo_dict = {}
    for fp in filepaths:
        contiginfo_dict.update(_getRefInfoFromFastaFile(fp))
    return contiginfo_dict


def _getRefInfoFromFastaFile(filepath):
    contiginfos = FastaInfo(filepath).getContigs()

    contiginfo_dict = {}
    for contiginfo in contiginfos:
        # contiginfo_dict[contiginfo.getContigName()] = contiginfo
        contiginfo_dict[contiginfo.getContigName()] = {}
        contiginfo_dict[contiginfo.getContigName()]['sequence'] = contiginfo.getSequence()
        contiginfo_dict[contiginfo.getContigName()]['seqLength'] = contiginfo.getSeqLength()
        contiginfo_dict[contiginfo.getContigName()]['md5'] = contiginfo.getMd5()
    return contiginfo_dict


class FastaInfo:

    def __init__(self, filepath):
        self._contigs = []  # list of ContigInfo()

        with open(filepath, mode='r') as rf:
            contigTmp = ContigInfo()
            firstline = next(rf)
            if not str(firstline).startswith(">"):
                print("read fasta file wrong!")
                raise ValueError
            else:
                contigTmp.setContigName(firstline.strip()[1:])

            sequencetmp = ""
            for line in rf:
                if line.startswith(">"):
                    # contigTmp.setSequence(sequencetmp)
                    contigTmp.setSeqLength(len(sequencetmp))
                    contigTmp.setMd5("")
                    self._contigs.append(contigTmp)

                    sequencetmp = ""
                    contigTmp = ContigInfo()
                    contigTmp.setContigName(line.strip()[1:])
                else:
                    sequencetmp += line.strip()

            # the last contig
            contigTmp.setSequence(sequencetmp)
            contigTmp.setSeqLength(len(sequencetmp))
            contigTmp.setMd5("")
            self._contigs.append(contigTmp)

    def getContigs(self):
        return self._contigs


class ContigInfo:
    def __init__(self):
        self._contigName = ""
        self._sequence = ""
        self._seqLength = 0
        self._md5 = ""

    def getContigName(self):
        return self._contigName

    def setContigName(self, contigname):
        self._contigName = contigname

    def getSequence(self):
        return self._sequence

    def setSequence(self, sequence):
        self._sequence = sequence

    def getSeqLength(self):
        return self._seqLength

    def setSeqLength(self, seqlength):
        self._seqLength = seqlength

    def getMd5(self):
        return self._md5

    def setMd5(self, md5):
        self._md5 = md5