#本文件使用unittest模块，测试了pysam的常用功能及其是否适用于一些特殊环境
import os
import pysam
import unittest
import json
import collections
import string
import struct
import copy
import array

from TestUtils import (
    checkFieldEqual,
    make_data_files,
    BAM_DATADIR,
    get_temp_filename,
    get_temp_context,
    IS_PYTHON3,
)


if IS_PYTHON3:
    maketrans = str.maketrans
else:
    maketrans = string.maketrans
#不同版本的python指令不同

def setUpModule():
    make_data_files(BAM_DATADIR)


class ReadTest(unittest.TestCase):
    def build_read(self):
        
        """使用pysam包生成测试集数据"""
        header = pysam.AlignmentHeader.from_references(
            ["chr1", "chr2"], [10000000, 10000000]
        )
       
        a = pysam.AlignedSegment(header)
        a.query_name = "read_12345"
        a.query_sequence = "ATGC" * 10
        a.flag = 0
        a.reference_id = 0
        a.reference_start = 20
        a.mapping_quality = 20
        a.cigartuples = ((0, 10), (2, 1), (0, 9), (1, 1), (0, 20))
        a.next_reference_id = 0
        a.next_reference_start = 200
        a.template_length = 167
        a.query_qualities = pysam.qualitystring_to_array("1234") * 10
        return a


class TestAlignedSegment(ReadTest):

    """tests to check if aligned read can be constructed
    and manipulated.
    """
    
    #测试空集的情况
    def testEmpty(self):

        a = pysam.AlignedSegment()
        self.assertEqual(a.query_name, None)
        self.assertEqual(a.query_sequence, None)
        self.assertEqual(pysam.qualities_to_qualitystring(a.query_qualities), None)
        self.assertEqual(a.flag, 0)
        self.assertEqual(a.reference_id, -1)
        self.assertEqual(a.mapping_quality, 0)
        self.assertEqual(a.cigartuples, None)
        self.assertEqual(a.get_tags(), [])
        self.assertEqual(a.next_reference_id, -1)
        self.assertEqual(a.next_reference_start, -1)
        self.assertEqual(a.template_length, 0)

    def testStrOfEmptyRead(self):
        a = pysam.AlignedSegment()
        s = str(a)
        self.assertEqual("None\t0\t*\t0\t0\tNone\t*\t0\t0\tNone\tNone\t[]", s)

    def testSettingTagInEmptyRead(self):
        """see issue 62"""
        a = pysam.AlignedSegment()
        a.tags = (("NM", 1),)
        a.query_qualities = None
        self.assertEqual(a.tags, [("NM", 1),])

    def testCompare(self):
        """check comparison functions."""
        a = self.build_read()
        b = None

        self.assertFalse(a is b)
        self.assertFalse(a == b)
        self.assertEqual(-1, a.compare(b))

        b = self.build_read()

        self.assertEqual(0, a.compare(b))
        self.assertEqual(0, b.compare(a))
        self.assertTrue(a == b)
        self.assertTrue(b == a)
        self.assertFalse(a != b)
        self.assertFalse(b != a)

        b.tid = 1
        self.assertFalse(a == b)
        self.assertFalse(b == a)
        self.assertTrue(a != b)
        self.assertTrue(b != a)


        # check qname
        b.query_name = "read_123"
        checkFieldEqual(self, a, b, "query_name")
        b.query_name = "read_12345678"
        checkFieldEqual(self, a, b, "query_name")
        b.query_name = "read_12345"
        checkFieldEqual(self, a, b)

        # check cigar
        b.cigartuples = ((0, 10),)
        checkFieldEqual(self, a, b, "cigartuples")
        b.cigartuples = ((0, 10), (2, 1), (0, 10))
        checkFieldEqual(self, a, b, "cigartuples")
        b.cigartuples = ((0, 10), (2, 1), (0, 9), (1, 1), (0, 20))
        checkFieldEqual(self, a, b)

        # check seq
        b.query_sequence = "ATGC"
        checkFieldEqual(
            self, a, b, ("query_sequence", "query_qualities", "query_length")
        )
        b.query_sequence = "ATGC" * 3
        checkFieldEqual(
            self, a, b, ("query_sequence", "query_qualities", "query_length")
        )
        b.query_sequence = "ATGC" * 10
        checkFieldEqual(self, a, b, ("query_qualities",))

        # reset qual
        b = self.build_read()

        def dual(name):
            if name.endswith('is_unmapped'): return name.replace('unmapped', 'mapped')
            elif name.endswith('is_mapped'): return name.replace('mapped', 'unmapped')
            elif name.endswith('is_reverse'): return name.replace('reverse', 'forward')
            elif name.endswith('is_forward'): return name.replace('forward', 'reverse')
            else: return name

        # check flags:
        for x in (
            "is_paired",
            "is_proper_pair",
            "is_unmapped",
            "mate_is_unmapped",
            "is_reverse",
            "mate_is_reverse",
            "is_read1",
            "is_read2",
            "is_secondary",
            "is_qcfail",
            "is_duplicate",
            "is_supplementary",
        ):
            setattr(b, x, True)
            self.assertEqual(getattr(b, x), True)
            checkFieldEqual(self, a, b, ("flag", x, dual(x),))
            setattr(b, x, False)
            self.assertEqual(getattr(b, x), False)
            checkFieldEqual(self, a, b)

        for x in (
            "is_mapped",
            "mate_is_mapped",
            "is_forward",
            "mate_is_forward",
        ):
            setattr(b, x, False)
            self.assertEqual(getattr(b, x), False)
            checkFieldEqual(self, a, b, ("flag", x, dual(x),))
            setattr(b, x, True)
            self.assertEqual(getattr(b, x), True)
            checkFieldEqual(self, a, b)

    #测试长的读长的情况
    def testLargeRead(self):
        """build an example read."""

        a = pysam.AlignedSegment()
        a.query_name = "read_12345"
        a.query_sequence = "ATGC" * 200
        a.flag = 0
        a.reference_id = -1
        a.reference_start = 20
        a.mapping_quality = 20
        a.cigartuples = ((0, 4 * 200),)
        a.next_reference_id = 0
        a.next_reference_start = 200
        a.template_length = 167
        a.query_qualities = pysam.qualitystring_to_array("1234") * 200

        self.assertTrue(a)


    def testPositions(self):
        a = self.build_read()
        self.assertEqual(
            a.get_reference_positions(),
            [
                20,
                21,
                22,
                23,
                24,
                25,
                26,
                27,
                28,
                29,
                31,
                32,
                33,
                34,
                35,
                36,
                37,
                38,
                39,
                40,
                41,
                42,
                43,
                44,
                45,
                46,
                47,
                48,
                49,
                50,
                51,
                52,
                53,
                54,
                55,
                56,
                57,
                58,
                59,
            ],
        )

        self.assertEqual(
            a.get_aligned_pairs(),
            [
                (0, 20),
                (1, 21),
                (2, 22),
                (3, 23),
                (4, 24),
                (5, 25),
                (6, 26),
                (7, 27),
                (8, 28),
                (9, 29),
                (None, 30),
                (10, 31),
                (11, 32),
                (12, 33),
                (13, 34),
                (14, 35),
                (15, 36),
                (16, 37),
                (17, 38),
                (18, 39),
                (19, None),
                (20, 40),
                (21, 41),
                (22, 42),
                (23, 43),
                (24, 44),
                (25, 45),
                (26, 46),
                (27, 47),
                (28, 48),
                (29, 49),
                (30, 50),
                (31, 51),
                (32, 52),
                (33, 53),
                (34, 54),
                (35, 55),
                (36, 56),
                (37, 57),
                (38, 58),
                (39, 59),
            ],
        )

        self.assertEqual(
            a.get_reference_positions(),
            [
                x[1]
                for x in a.get_aligned_pairs()
                if x[0] is not None and x[1] is not None
            ],
        )
        # alen is the length of the aligned read in genome
        self.assertEqual(a.reference_length, a.get_aligned_pairs()[-1][0] + 1)
        # aend points to one beyond last aligned base in ref
        self.assertEqual(a.get_reference_positions()[-1], a.reference_end - 1)

#测试是否能检测不同cigar值
    def test_infer_query_length(self):
        """Test infer_query_length on M|=|X|I|D|H|S cigar ops"""
        a = self.build_read()
        a.cigarstring = "40M"
        self.assertEqual(a.infer_query_length(), 40)
        a.cigarstring = "40="
        self.assertEqual(a.infer_query_length(), 40)
        a.cigarstring = "40X"
        self.assertEqual(a.infer_query_length(), 40)
        a.cigarstring = "20M5I20M"
        self.assertEqual(a.infer_query_length(), 45)
        a.cigarstring = "20M5D20M"
        self.assertEqual(a.infer_query_length(), 40)
        a.cigarstring = "5H35M"
        self.assertEqual(a.infer_query_length(), 35)
        a.cigarstring = "5S35M"
        self.assertEqual(a.infer_query_length(), 40)
        a.cigarstring = "35M5H"
        self.assertEqual(a.infer_query_length(), 35)
        a.cigarstring = "35M5S"
        self.assertEqual(a.infer_query_length(), 40)
        a.cigarstring = None
        self.assertEqual(a.infer_query_length(), None)

    def test_infer_read_length(self):
        """Test infer_read_length on M|=|X|I|D|H|S cigar ops"""
        a = self.build_read()
        a.cigarstring = "40M"
        self.assertEqual(a.infer_read_length(), 40)
        a.cigarstring = "40="
        self.assertEqual(a.infer_read_length(), 40)
        a.cigarstring = "40X"
        self.assertEqual(a.infer_read_length(), 40)
        a.cigarstring = "20M5I20M"
        self.assertEqual(a.infer_read_length(), 45)
        a.cigarstring = "20M5D20M"
        self.assertEqual(a.infer_read_length(), 40)
        a.cigarstring = "5H35M"
        self.assertEqual(a.infer_read_length(), 40)
        a.cigarstring = "5S35M"
        self.assertEqual(a.infer_read_length(), 40)
        a.cigarstring = "35M5H"
        self.assertEqual(a.infer_read_length(), 40)
        a.cigarstring = "35M5S"
        self.assertEqual(a.infer_read_length(), 40)
        a.cigarstring = None
        self.assertEqual(a.infer_read_length(), None)

    def test_get_aligned_pairs_soft_clipping(self):
        a = self.build_read()
        a.cigartuples = ((4, 2), (0, 35), (4, 3))
        self.assertEqual(
            a.get_aligned_pairs(),
            [(0, None), (1, None)]
            + [
                (qpos, refpos)
                for (qpos, refpos) in zip(range(2, 2 + 35), range(20, 20 + 35))
            ]
            + [(37, None), (38, None), (39, None)],
        )
        self.assertEqual(
            a.get_aligned_pairs(True),
            # [(0, None), (1, None)] +
            [
                (qpos, refpos)
                for (qpos, refpos) in zip(range(2, 2 + 35), range(20, 20 + 35))
            ]
            # [(37, None), (38, None), (39, None)]
        )

   
    def test_get_aligned_pairs_match_mismatch(self):
        a = self.build_read()
        a.cigartuples = ((7, 20), (8, 20))
        self.assertEqual(
            a.get_aligned_pairs(),
            [
                (qpos, refpos)
                for (qpos, refpos) in zip(range(0, 0 + 40), range(20, 20 + 40))
            ],
        )
        self.assertEqual(
            a.get_aligned_pairs(True),
            [
                (qpos, refpos)
                for (qpos, refpos) in zip(range(0, 0 + 40), range(20, 20 + 40))
            ],
        )

    def test_get_aligned_pairs_padding(self):
        a = self.build_read()
        a.cigartuples = ((0, 1), (6, 1), (0, 1))
        # The padding operation is like an insertion into the reference.
        # See comment in test_get_aligned_pairs_padding_with_seq (below).
        self.assertEqual(a.get_aligned_pairs(),
                         [(0, 20), (1, None), (2, 21)])


    def test_get_aligned_pairs(self):
        a = self.build_read()
        a.query_sequence = "A" * 9
        a.cigarstring = "9M"
        a.set_tag("MD", "9")
        self.assertEqual(
            a.get_aligned_pairs(with_seq=True),
            [
                (0, 20, "A"),
                (1, 21, "A"),
                (2, 22, "A"),
                (3, 23, "A"),
                (4, 24, "A"),
                (5, 25, "A"),
                (6, 26, "A"),
                (7, 27, "A"),
                (8, 28, "A"),
            ],
        )

        a.set_tag("MD", "4C4")
        self.assertEqual(
            a.get_aligned_pairs(with_seq=True),
            [
                (0, 20, "A"),
                (1, 21, "A"),
                (2, 22, "A"),
                (3, 23, "A"),
                (4, 24, "c"),
                (5, 25, "A"),
                (6, 26, "A"),
                (7, 27, "A"),
                (8, 28, "A"),
            ],
        )

        a.cigarstring = "5M2D4M"
        a.set_tag("MD", "4C^TT4")
        self.assertEqual(
            a.get_aligned_pairs(with_seq=True),
            [
                (0, 20, "A"),
                (1, 21, "A"),
                (2, 22, "A"),
                (3, 23, "A"),
                (4, 24, "c"),
                (None, 25, "T"),
                (None, 26, "T"),
                (5, 27, "A"),
                (6, 28, "A"),
                (7, 29, "A"),
                (8, 30, "A"),
            ],
        )

        a.cigarstring = "5M2D2I2M"
        a.set_tag("MD", "4C^TT2")
        self.assertEqual(
            a.get_aligned_pairs(with_seq=True),
            [
                (0, 20, "A"),
                (1, 21, "A"),
                (2, 22, "A"),
                (3, 23, "A"),
                (4, 24, "c"),
                (None, 25, "T"),
                (None, 26, "T"),
                (5, None, None),
                (6, None, None),
                (7, 27, "A"),
                (8, 28, "A"),
            ],
        )


class TestAlignedPairs(unittest.TestCase):
    filename = os.path.join(BAM_DATADIR, "example_aligned_pairs.bam")

    def testReferenceBases(self):
        """reference bases should always be the same nucleotide
        """
        reference_bases = collections.defaultdict(list)
        with pysam.AlignmentFile(self.filename) as inf:
            for c in inf.pileup():
                for r in c.pileups:
                    for read, ref, base in r.alignment.get_aligned_pairs(with_seq=True):
                        if ref is None:
                            continue
                        reference_bases[ref].append(base.upper())

        for x, y in reference_bases.items():
            self.assertEqual(len(set(y)), 1)



    def testOrient(self):
        """reference bases should always be the same nucleotide
        """
        filename = os.path.join(BAM_DATADIR, "MM-orient.bam")
        expect = {
            "top-fwd": [
                {("C", 0, "m"): [(7, 128), (30, 153), (31, 179)]},
                {("C", 0, "m"): [(7, 128), (30, 153), (31, 179)]},
            ],
            "top-rev": [
                {("C", 1, "m"): [(4, 179), (5, 153), (28, 128)]},
                {("C", 0, "m"): [(31, 179), (30, 153), (7, 128)]},
            ],
            "bot-fwd": [
                {("G", 1, "m"): [(1, 115), (2, 141), (18, 166), (23, 192)]},
                {("G", 1, "m"): [(1, 115), (2, 141), (18, 166), (23, 192)]},
            ],
            "bot-rev": [
                {("G", 0, "m"): [(12, 192), (17, 166), (33, 141), (34, 115)]},
                {("G", 1, "m"): [(23, 192), (18, 166), (2, 141), (1, 115)]},
            ],
        }

        with pysam.AlignmentFile(filename, check_sq=False) as inf:
            for r in inf:
                self.assertDictEqual(r.modified_bases, expect[r.query_name][0])
                self.assertDictEqual(r.modified_bases_forward, expect[r.query_name][1])
                for (B, s, _), mods in r.modified_bases.items():
                    C = B.translate(maketrans("ACGTacgtNnXx", "TGCAtgcaNnXx"))
                    for pos, _ in mods:
                        if r.is_reverse:
                            if s == 1:
                                self.assertEqual(
                                    C, r.query_sequence[pos], r.to_string()
                                )
                            else:
                                self.assertEqual(
                                    C, r.query_sequence[pos], r.to_string()
                                )
                        else:
                            if s == 0:
                                self.assertEqual(
                                    B, r.query_sequence[pos], r.to_string()
                                )
                            else:
                                self.assertEqual(
                                    B, r.query_sequence[pos], r.to_string()
                                )


class TestTags(ReadTest):
    def testMissingTag(self):
        a = self.build_read()
        self.assertRaises(KeyError, a.get_tag, "XP")

    def testEmptyTag(self):
        a = self.build_read()
        self.assertRaises(KeyError, a.get_tag, "XT")

    def testSetTag(self):
        a = self.build_read()
        self.assertEqual(False, a.has_tag("NM"))
        a.set_tag("NM", 2)
        self.assertEqual(True, a.has_tag("NM"))
        self.assertEqual(a.get_tag("NM"), 2)
        a.set_tag("NM", 3)
        self.assertEqual(a.get_tag("NM"), 3)
        a.set_tag("NM", None)
        self.assertEqual(False, a.has_tag("NM"))
        # check if deleting a non-existing tag is fine
        a.set_tag("NM", None)
        a.set_tag("NM", None)

