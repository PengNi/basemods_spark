#################################################################################
# Copyright (c) 2011-2013, Pacific Biosciences of California, Inc.
#
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
# * Redistributions of source code must retain the above copyright
#   notice, this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright
#   notice, this list of conditions and the following disclaimer in the
#   documentation and/or other materials provided with the distribution.
# * Neither the name of Pacific Biosciences nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE GRANTED BY
# THIS LICENSE.  THIS SOFTWARE IS PROVIDED BY PACIFIC BIOSCIENCES AND ITS
# CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR
# ITS CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
# BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
# IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#################################################################################

class CmpH5Format:
     # def __init__(self, cmpH5):
     def __init__(self):
        # if ('Version' in cmpH5.attrs):
        #     self.VERSION = cmpH5.attrs['Version']

        self.ALN_INFO             = 'AlnInfo'
        self.REF_INFO             = 'RefInfo'
        self.MOVIE_INFO           = 'MovieInfo'
        self.REF_GROUP            = 'RefGroup'
        self.ALN_GROUP            = 'AlnGroup'
        self.ALN_INDEX_NAME       = 'AlnIndex'
        self.FILE_LOG             = 'FileLog'
        self.BARCODE_INFO         = 'BarcodeInfo'

        self.ALN_INDEX            = '/'.join([self.ALN_INFO, self.ALN_INDEX_NAME])
        self.REF_GROUP_ID         = '/'.join([self.REF_GROUP, 'ID'])
        self.REF_GROUP_PATH       = '/'.join([self.REF_GROUP, 'Path'])
        self.REF_GROUP_INFO_ID    = '/'.join([self.REF_GROUP, 'RefInfoID'])

        self.REF_OFFSET_TABLE     = '/'.join([self.REF_GROUP, 'OffsetTable'])
        self.ALN_GROUP_ID         = '/'.join([self.ALN_GROUP, 'ID'])
        self.ALN_GROUP_PATH       = '/'.join([self.ALN_GROUP, 'Path'])

        # Movie Info
        self.MOVIE_INFO_ID         = '/'.join([self.MOVIE_INFO, 'ID'])
        self.MOVIE_INFO_NAME       = '/'.join([self.MOVIE_INFO, 'Name'])
        self.MOVIE_INFO_EXP        = '/'.join([self.MOVIE_INFO, 'Exp'])
        self.MOVIE_INFO_FRAME_RATE = '/'.join([self.MOVIE_INFO, 'FrameRate'])
        self.MOVIE_INFO_RUN        = '/'.join([self.MOVIE_INFO, 'Run'])
        self.MOVIE_INFO_BINDING_KIT      = '/'.join([self.MOVIE_INFO, 'BindingKit'])
        self.MOVIE_INFO_SEQUENCING_KIT   = '/'.join([self.MOVIE_INFO, 'SequencingKit'])
        self.MOVIE_INFO_SOFTWARE_VERSION = '/'.join([self.MOVIE_INFO, 'SoftwareVersion'])

        (self.ID, self.ALN_ID, self.MOVIE_ID, self.REF_ID, self.TARGET_START,
         self.TARGET_END, self.RC_REF_STRAND, self.HOLE_NUMBER, self.SET_NUMBER,
         self.STROBE_NUMBER, self.MOLECULE_ID, self.READ_START, self.READ_END,
         self.MAP_QV, self.N_MATCHES, self.N_MISMATCHES, self.N_INSERTIONS,
         self.N_DELETIONS, self.OFFSET_BEGIN, self.OFFSET_END, self.N_BACK,
         self.N_OVERLAP) = range(0, 22)

        # self.extraTables = ['/'.join([self.ALN_INFO, x]) for x in
        #                     cmpH5[self.ALN_INFO].keys()
        #                     if not x == self.ALN_INDEX_NAME]
        # sorting
        self.INDEX_ATTR = "Index"
        self.INDEX_ELTS = ['REF_ID', 'TARGET_START', 'TARGET_END']
